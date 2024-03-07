#include "LennardJonesSystem.h"
#include "Application.h"
#include <vtkm/Pair.h>
#include <vtkm/VectorAnalysis.h>
#include <vtkm/cont/Algorithm.h>
#include <vtkm/cont/ArrayHandleGroupVecVariable.h>
#include <vtkm/cont/ConvertNumComponentsToOffsets.h>
#include <vtkm/worklet/WorkletMapField.h>
#include <vtkm/worklet/WorkletMapTopology.h>
#include <fstream>
#include "math/Utiles.h"
#include "system/worklet/SystemWorklet.h"
#include "MeshFreeCondition.h"
#include "FieldName.h"
#include <vtkm/VecVariable.h>
#include <vtkm/cont/ArrayHandleConstant.h>

//RegisterObject(LennardJonesSystem);

LennardJonesSystem::LennardJonesSystem(const Configuration& cfg)
  : MDSystem(cfg)
  , _nearforce_type(Get<std::string>("nearforce_type"))
{
  SetParameter(PARA_KMAX, Get<IdComponent>("kmax"));
  SetParameter(PARA_ALPHA, Get<Real>("alpha"));
  SetParameter(PARA_TEMPT_SUM, Real{ 0.0 });
  SetParameter(PARA_TEMPT, Real{ 0.0 });
}

void LennardJonesSystem::Init()
{
  MDSystem::Init();

  //init variable
  InitialCondition();

  ComputeForce();
}

void LennardJonesSystem::InitialCondition()
{
  MDSystem::InitialCondition();

  _nosehooverxi = 0.0;
  _rho = GetParameter<Real>(PARA_RHO);
  SetParameter(PARA_RDF_RHO, _rho);
  SetCharge();

  vtkm::Vec<vtkm::Range, 3> range = GetParameter<vtkm::Vec<vtkm::Range, 3>>(PARA_RANGE);
  vtkm::Vec3f left_bottom{
    { static_cast<vtkm::FloatDefault>(range[0].Min),},
    { static_cast<vtkm::FloatDefault>(range[1].Min),},
    { static_cast<vtkm::FloatDefault>(range[2].Min),}
  };
  vtkm::Vec3f right_top{
    { static_cast<vtkm::FloatDefault>(range[0].Max),},
    { static_cast<vtkm::FloatDefault>(range[1].Max),},
    { static_cast<vtkm::FloatDefault>(range[2].Max),}
  };
}

void LennardJonesSystem::ComputeForce()
{
  if (_nearforce_type == "RBL")
  {
    ComputeRBLLJForce(_force);
  }
  else if (_nearforce_type == "VERLETLIST")
  {
    ComputeVerletlistLJForce(_force);
  }
  else if (_nearforce_type == "ORIGINAL")
  {
    ComputeOriginalLJForce(_force);
  }
}

void LennardJonesSystem::UpdateVelocity()
{
  try
  {
    SystemWorklet::UpdateVelocity(_dt, _unit_factor._fmt2v, _force, _mass, _velocity);
    ComputeTempe();

    Real _kbT = 1.0;
    Real dt_divide_taut = 0.02;
    Real coeff_Berendsen = vtkm::Sqrt(1.0 + dt_divide_taut * (_kbT / _tempT - 1.0));
    SystemWorklet::UpdateVelocityRescale(coeff_Berendsen, _velocity);
  }
  catch (const std::exception& e)
  {
    std::cout << e.what() << std::endl;
  }
}

void LennardJonesSystem::UpdatePosition()
{
  //SystemWorklet::UpdatePosition(_dt,_velocity, _locator, _position);
  auto&& position_flag = GetFieldAsArrayHandle<Id3>(field::position_flag);
  SystemWorklet::UpdatePositionFlag(_dt, _velocity, _locator, _position, position_flag);
  _locator.SetPosition(_position);
  SetCenterTargetPositions();
}

void LennardJonesSystem::SetCharge()
{
  _charge = GetFieldAsArrayHandle<Real>(field::charge);
  auto n = _position.GetNumberOfValues();
  _charge.AllocateAndFill(n, 0);
  _topology.SetCharge(_charge);
}

void LennardJonesSystem::ComputeTempe()
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

void LennardJonesSystem::PreSolve()
{
  _locator.SetPosition(_position);
  std::shared_ptr<Executioner>& executioner = _app.GetExecutioner();
  _dt = executioner->Dt();
}

void LennardJonesSystem::InitField() 
{
  MDSystem::InitField();
  AddField(field::pts_type , ArrayHandle<Id>{});
  AddField(field::position_flag, ArrayHandle<Id3>{});
  AddField(field::center_position, ArrayHandle<Vec3f>{});
  AddField(field::target_position, ArrayHandle<Vec3f>{});
  AddField(field::epsilon, ArrayHandle<Real>{});
  AddField(field::sigma, ArrayHandle<Real>{});
}

void LennardJonesSystem::TimeIntegration() {}

void LennardJonesSystem::Solve() 
{
  // stage1:
  UpdateVelocity();

  // stage2:
  UpdatePosition();

  // stage3:
  ComputeForce();
  UpdateVelocity();
}

void LennardJonesSystem::ESS_RBL()
{
  ComputeTempe();
  auto N = _position.GetNumberOfValues();
  vtkm::Float64 kinetic_energy = _tempT_sum / (2.0 * N);
  std::shared_ptr<Executioner>& executioner = _app.GetExecutioner();
  _dt = executioner->Dt();
  vtkm::Float64 gamma_ess = 10.0 * _dt;
  vtkm::Float64 xi_ess = vtkm::Sqrt(
    vtkm::Abs(1.0 + _dt / (gamma_ess * kinetic_energy) * (-4.78 - _rbl_pe_ave - kinetic_energy)));
  //std::cout << "kinetic_energy = "<< kinetic_energy << "  rbl_pe_ave = " << _rbl_pe_ave << "  xi_ess = " << xi_ess << std::endl;
  SystemWorklet::UpdateVelocityRescale(xi_ess, _velocity);
}

void LennardJonesSystem::PostSolve()
{

}

void LennardJonesSystem::SetCenterTargetPositions()
{
  auto center_position = GetFieldAsArrayHandle<Vec3f>(field::center_position);
  vtkm::cont::ArrayCopy(_position, center_position);

  auto target_position = GetFieldAsArrayHandle<Vec3f>(field::target_position);
  vtkm::cont::ArrayCopy(_position, target_position);
}

void LennardJonesSystem::SetForceFunction()
{
  //auto cut_off = GetParameter<Real>(PARA_CUTOFF);
  //_force_function.SetCutOff(cut_off);
}

void LennardJonesSystem::SetTopology()
{
  auto pts_type = GetFieldAsArrayHandle<Id>(field::pts_type);
  auto molecule_id = GetFieldAsArrayHandle<Id>(field::molecule_id);
  auto epsilon = GetFieldAsArrayHandle<Real>(field::epsilon);
  auto sigma = GetFieldAsArrayHandle<Real>(field::sigma);
  _topology.SetAtomsType(pts_type);
  _topology.SetMolecularId(molecule_id);
  _topology.SetEpsAndSigma(epsilon, sigma);
}
