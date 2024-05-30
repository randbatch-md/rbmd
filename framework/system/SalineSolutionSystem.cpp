#include "SalineSolutionSystem.h"
#include "Application.h"
#include <fstream>
#include <vtkm/Pair.h>
#include <vtkm/VectorAnalysis.h>
#include <vtkm/cont/Algorithm.h>
#include <vtkm/cont/ArrayHandleGroupVecVariable.h>
#include <vtkm/cont/ConvertNumComponentsToOffsets.h>
#include <vtkm/Math.h>
#include <cmath>  // erfc(x)
#include "math/Math.h"
#include "math/Utiles.h"
#include "RBEPSample.h"
#include "system/worklet/SystemWorklet.h"
#include "MeshFreeCondition.h"
#include "FieldName.h"

//RegisterObject(SalineSolutionSystem);

SalineSolutionSystem::SalineSolutionSystem(const Configuration& cfg)
  : MDSystem(cfg) 
  , _executioner((_app.GetExecutioner()))
  , _kbT(Get<IdComponent>("kbT"))
  , _Kmax(Get<IdComponent>("kmax"))
  , _alpha(Get<Real>("alpha"))
  , _farforce_type(Get<std::string>("farforce_type"))
  , _nearforce_type(Get<std::string>("nearforce_type"))
  , _temp_con_type(Get<std::string>("temp_con_type"))
{
  SetParameter(PARA_RBE_P, _RBE_P);
  SetParameter(PARA_ALPHA, _alpha);
  SetParameter(PARA_KMAX, _Kmax);
  SetParameter(PARA_TEMPT_SUM, Real{ 0.0 });
  SetParameter(PARA_TEMPT, Real{ 0.0 });
}

void SalineSolutionSystem::Init()
{
  MDSystem::Init();

  //init variable
  InitialCondition();

  ComputeForce(); // Presolve force

  _Elefartimer_counting = 0.0;
}

void SalineSolutionSystem::InitialCondition()
{
  MDSystem::InitialCondition();
  _rho = GetParameter<Real>(PARA_RHO);
  _nosehooverxi = 0.0;

  SetParameter(PARA_RDF_RHO, _rho / Real{ 2.0 });
  SetCharge();
  PreForce();
}

void SalineSolutionSystem::ComputeForce()
{
  ComputeAllForce();
  TempConTypeForce();
}

void SalineSolutionSystem::ComputeAllForce()
{
  //SystemWorklet::SumFarNearLJForce(EleNewForce(),  EleNearForce(), LJForce(), _all_force); //RBE + LJ
    SystemWorklet::SumFarNearForce(EleNewForce(), NearForce(), _all_force); //RBE + LJ
}

void SalineSolutionSystem::UpdateVelocity()
{
  try
  {
    SystemWorklet::UpdateVelocity(_dt, _unit_factor._fmt2v ,_all_force, _mass, _velocity);
    ComputeTempe();
    UpdateVelocityByTempConType();
  }
  catch (const std::exception& e)
  {
    std::cout << e.what() << std::endl;
  }
}

void SalineSolutionSystem::UpdatePosition()
{
  //SystemWorklet::UpdatePosition(_dt, _velocity, _locator, _position);
  auto&& position_flag = GetFieldAsArrayHandle<Id3>(field::position_flag);
  SystemWorklet::UpdatePositionFlag(_dt, _velocity, _locator, _position, position_flag);
  _locator.SetPosition(_position);
  SetCenterTargetPositions();
}

void SalineSolutionSystem::PreForce()
{
  _Vlength = GetParameter<Real>(PARA_VLENGTH);
  _dt = _executioner->Dt();
  // prepare for RBE force
  auto velocity_type = GetParameter<std::string>(gtest::velocity_type);
  auto random = (velocity_type != "TEST") ? true : false;
  RBEPSAMPLE rbe_presolve_psample = { _alpha, _Vlength, _RBE_P };
  rbe_presolve_psample._RBE_random = random;
  _psample = rbe_presolve_psample.Fetch_P_Sample(
    Real(0.0), (vtkm::Sqrt(_alpha / 2.0) * _Vlength / vtkm::Pi()));

  // prepare for Langevin dynamics
  RBEPSAMPLE sample_presolve_1d;
  sample_presolve_1d._RBE_random = random;
  Real gaussianx = sample_presolve_1d.FetchSample_1D(Real(0.0), Real(1.0));
  Real gaussiany = sample_presolve_1d.FetchSample_1D(Real(0.0), Real(1.0));
  Real gaussianz = sample_presolve_1d.FetchSample_1D(Real(0.0), Real(1.0));
  _gaussian = { gaussianx, gaussiany, gaussianz };
}

void SalineSolutionSystem::SetCharge()
{
  _charge = GetFieldAsArrayHandle<Real>(field::charge);
  auto n = _position.GetNumberOfValues();
  _charge.Allocate(n);
  _charge.Fill(-1.0, 0);
  _charge.Fill(1.0, n / 2);
  _topology.SetCharge(_charge);
}

vtkm::cont::ArrayHandle<Vec3f> SalineSolutionSystem::LJForce()
{
  if (_nearforce_type == "RBL")
  {
    ComputeRBLLJForce(_LJforce);
  }
  else if (_nearforce_type == "VERLETLIST")
  {
    ComputeVerletlistLJForce(_LJforce);
  }
  else if (_nearforce_type == "ORIGINAL")
  {
    ComputeOriginalLJForce(_LJforce);
  }
  return _LJforce;
}

vtkm::cont::ArrayHandle<Vec3f> SalineSolutionSystem::EleNearForce()
{
  SystemWorklet::ComputeNearElectrostatics(_atoms_id, _locator, _topology,_force_function, _ele_near_force);
  return _ele_near_force;
}

vtkm::cont::ArrayHandle<Vec3f> SalineSolutionSystem::EleNewForce()
{
  _EleFartimer.Start();
  if (_farforce_type == "RBE")
  {
    // New RBE force part
    ComputeRBEEleForce(_psample, _RBE_P, _ele_new_force);
    _Elefartimer_counting = _Elefartimer_counting + _EleFartimer.GetElapsedTime();
    std::cout << "RBE time: " << _Elefartimer_counting << std::endl;
  }
  if (_farforce_type == "EWALD")
  {
    // New Ewald far part
    ComputeEwaldEleForce(_Kmax, _ele_new_force);
  }
  return _ele_new_force;
}

vtkm::cont::ArrayHandle<Vec3f> SalineSolutionSystem::NearForce()
{
  if (_nearforce_type == "RBL")
  {
    ComputeRBLNearForce(_nearforce);
  }
  else if (_nearforce_type == "VERLETLIST")
  {
    ComputeVerletlistNearForce(_nearforce);
  }
  else if (_nearforce_type == "ORIGINAL")
  {
    SystemWorklet::SumFarNearForce(EleNearForce(), LJForce(), _nearforce);
  }
  return _nearforce;
}

void SalineSolutionSystem::TempConTypeForce()
{
  vtkm::cont::ArrayHandle<Real> mass;
  mass.Allocate(_all_force.GetNumberOfValues());
  auto writePortal = mass.WritePortal();
  for (size_t i = 0; i < _all_force.GetNumberOfValues(); i++)
  {
    writePortal.Set(i, 1);
  }
  if (_temp_con_type == "LANGEVIN")
  {
    // Underdamped Langevin
    Real kBT = 1.0;
    Real gamma = 100.0;
    //Maybe gamma = 50 is a good choice for dt = 2e-3. The choice of gamma depends on the value of dt.
    //In another words, gamma * dt can not be too small.
    //auto&& velocity = GetFieldAsArrayHandle<Vec3f>(field::velocity);
    SystemWorklet::UnderdampedLangevin(_gaussian, kBT, gamma, _dt, mass, _velocity, _all_force);
  }
}

void SalineSolutionSystem::ComputeTempe()
{
  auto n = _position.GetNumberOfValues();
  ArrayHandle<Real> sq_velocity;
  sq_velocity.Allocate(n);

  SystemWorklet::ComputerKineticEnergy(_velocity, _mass, sq_velocity);
  _tempT_sum =
    vtkm::cont::Algorithm::Reduce(sq_velocity, vtkm::TypeTraits<Real>::ZeroInitialization());
  _tempT = 0.5 * _tempT_sum / (3 * n / 2.0);
  SetParameter(PARA_TEMPT_SUM, _tempT_sum);
  SetParameter(PARA_TEMPT, _tempT);
}

void SalineSolutionSystem::UpdateVelocityByTempConType()
{
    if(_temp_con_type == "NOSE_HOOVER")
    {
      //Nose Hoover
      //Maybe tauT = 20.0 * dt is a good choice for dt = 2e-3 from the test.
      //Because the temperature curve is the smoothest of all test simulations.
      //In fact, 5.0 and 10.0 are also optional.
      //As long as the coefficent is not too large, such as larger than 100 * dt.
      SystemWorklet::UpdateVelocityNoseHoover(
        _dt, _unit_factor._fmt2v, _nosehooverxi, _all_force, _mass, _velocity);
      Real tauT = 20.0 * _dt;
      _nosehooverxi += 0.5 * _dt * (_tempT / _kbT - 1.0) / tauT;
    }
    else if(_temp_con_type == "TEMP_RESCALE")
    {
      Real coeff_rescale = vtkm::Sqrt(_kbT / _tempT);
      SystemWorklet::UpdateVelocityRescale(coeff_rescale, _velocity);
    }
    else if(_temp_con_type == "BERENDSEN")
    {
      //
      //Velocity Rescale: Berendsen
      //Maybe dt_divide_taut = 0.05 is a good choice for dt = 2e-3. 0.005, 0.01, 0.1 is optional.
      //The selection of dt_divide_taut determines the temperature equilibrium time.
      //
      Real dt_divide_taut = 0.005;
      Real coeff_Berendsen = vtkm::Sqrt(1.0 + dt_divide_taut * (_kbT / _tempT - 1.0));
      SystemWorklet::UpdateVelocityRescale(coeff_Berendsen, _velocity);
    }
}

void SalineSolutionSystem::PreSolve()
{
  _locator.SetPosition(_position);
  PreForce();
}
 
void SalineSolutionSystem::Solve() 
{
  // stage1:
  UpdateVelocity();

  // stage2:
  UpdatePosition();

  // stage3:
  ComputeForce();
  UpdateVelocity();
}

void SalineSolutionSystem::PostSolve()
{}

void SalineSolutionSystem::InitField()
{
  MDSystem::InitField();
  AddField(field::pts_type , ArrayHandle<Id>{});
  AddField(field::position_flag, ArrayHandle<Id3>{});
  AddField(field::center_position, ArrayHandle<Vec3f>{});
  AddField(field::target_position, ArrayHandle<Vec3f>{});
  AddField(field::epsilon, ArrayHandle<Real>{});
  AddField(field::sigma, ArrayHandle<Real>{});
}

void SalineSolutionSystem::SetCenterTargetPositions()
{
  auto num_pos = _position.GetNumberOfValues();
  vtkm::cont::ArrayHandle<vtkm::Vec3f> center_position_temp;
  center_position_temp.Allocate(num_pos / 2);
  auto&& write_prot_center = center_position_temp.WritePortal();
  auto&& read_prot_center = _position.ReadPortal();

  vtkm::cont::ArrayHandle<vtkm::Vec3f> target_position_temp;
  target_position_temp.Allocate(num_pos / 2);
  auto&& write_prot_target = target_position_temp.WritePortal();
  auto&& read_prot_target = _position.ReadPortal();

  for (int i = 0; i < num_pos / 2; i++)
  {
    write_prot_center.Set(i, read_prot_center.Get(i));
    write_prot_target.Set(i, read_prot_target.Get(i + (num_pos / 2)));
  }

  auto center_position = GetFieldAsArrayHandle<Vec3f>(field::center_position);
  vtkm::cont::ArrayCopy(center_position_temp, center_position);

  auto target_position = GetFieldAsArrayHandle<Vec3f>(field::target_position);
  vtkm::cont::ArrayCopy(target_position_temp, target_position);
}

void SalineSolutionSystem::TimeIntegration() {}

void SalineSolutionSystem::SetForceFunction()
{
  InitERF();
  auto cut_off = GetParameter<Real>(PARA_CUTOFF);
  auto alpha = GetParameter<Real>(PARA_ALPHA);
  auto volume = GetParameter<Real>(PARA_VOLUME);
  auto vlength = GetParameter<Real>(PARA_VLENGTH);
  auto Kmax = GetParameter<IdComponent>(PARA_KMAX);

  _force_function.SetParameters(cut_off, alpha, volume, vlength, Kmax);
}

void SalineSolutionSystem::SetTopology()
{
  auto pts_type = GetFieldAsArrayHandle<Id>(field::pts_type);
  auto molecule_id = GetFieldAsArrayHandle<Id>(field::molecule_id);
  auto epsilon = GetFieldAsArrayHandle<Real>(field::epsilon);
  auto sigma = GetFieldAsArrayHandle<Real>(field::sigma);
  _topology.SetAtomsType(pts_type);
  _topology.SetMolecularId(molecule_id);
  _topology.SetEpsAndSigma(epsilon, sigma);
}