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
  SetParameter(PARA_PRESSURE, Real{ 0.0 });
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
  set_global_box();
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

  //
  fix_press_berendsen();
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
  _tempT = 0.5 * _tempT_sum / ((3 * n -3) / 2.0);
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

void LennardJonesSystem::fix_press_berendsen()
{
  bulkmodulus = 10.0;
  p_start[0] = p_start[1] = p_start[2] = 10.0;
  p_stop[0] = p_stop[1] = p_stop[2] = 10.0;
  p_period[0] = p_period[1] = p_period[2] = 20;

  // compute new T,P

  Compute_Pressure_Scalar();
  Couple();

  auto currentstep = _app.GetExecutioner()->CurrentStep();
  auto beginstep = 0;
  auto endstep = _app.GetExecutioner()->NumStep();

  auto delta = currentstep - beginstep;
  if (delta != 0.0)
  {
    delta = delta / (endstep - beginstep);
  }
  for (int i = 0; i < 3; i++)
  {
    p_target[i] = p_start[i] + delta * (p_stop[i] - p_start[i]);
    dilation[i] =
      pow(1.0 - _dt / p_period[i] * (p_target[i] - p_current[i]) / bulkmodulus, 1.0 / 3.0);
  }

  //for (int i = 0; i < 3; i++)
  //{
  //  std::cout << "i=" << i  << ",dilation=" << dilation[i] << std::endl;
  //}
  // remap simulation box and atoms
  // redo KSpace coeffs since volume has changed

  remap();
}

void LennardJonesSystem::ComputeVirial()
{
  auto cut_off = GetParameter<Real>(PARA_CUTOFF);

  SystemWorklet::LJVirial(
   cut_off, _atoms_id, _locator, _topology, _force_function, _virial_atom);

  //auto range = GetParameter<vtkm::Vec<vtkm::Range, 3>>(PARA_RANGE);
  //auto Vlength = range[0].Max - range[0].Min;
  //SystemWorklet::LJVirialPBC(
   // cut_off, Vlength,_atoms_id, _locator, _topology, _force_function, _virial_atom);


  //for (int i = 0; i <_virial_atom.GetNumberOfValues();++i)
  //{
  //  std::cout << "i=" << i << ",_virial_atom=" << _virial_atom.ReadPortal().Get(i)[0] << ","
  //            << _virial_atom.ReadPortal().Get(i)[1] << "," << _virial_atom.ReadPortal().Get(i)[2] << ","
  //            << _virial_atom.ReadPortal().Get(i)[3] << "," << _virial_atom.ReadPortal().Get(i)[4] << ","
  //            << _virial_atom.ReadPortal().Get(i)[5]
  //      << std::endl;
  //}

  //reduce virial_atom
  virial = { 0, 0, 0, 0, 0, 0 };
  for (int i = 0; i < _virial_atom.GetNumberOfValues(); ++i)
  {
    // 获取当前原子的virial
    Vec6f vatom = _virial_atom.ReadPortal().Get(i);
    for (int j = 0; j < 6; ++j)
    {
      virial[j] += vatom[j];
    }
  }
  //
  //for (int i = 0; i < 6; ++i)
  //{
  //  std::cout << "total_virial[" << i << "] = " << virial[i] << std::endl;
  //}
}

void LennardJonesSystem::Compute_Pressure_Scalar()
{
  //compute temperature
  auto temperature = GetParameter<Real>(PARA_TEMPT);

  // compute  virial
  vtkm::Vec<vtkm::Range, 3> range = GetParameter<vtkm::Vec<vtkm::Range, 3>>(PARA_RANGE);
  auto volume =
    (range[0].Max - range[0].Min) * (range[1].Max - range[1].Min) * (range[2].Max - range[2].Min);

  auto inv_volume = 1.0 / volume;
  ComputeVirial();

  //compute dof
  auto n = _position.GetNumberOfValues();
  auto extra_dof = 3; //dimension =3
  auto dof = 3 * n - extra_dof;

  //compute pressure_scalar
  _pressure_scalar = (dof * _unit_factor.boltz * temperature + virial[0] + virial[1] + virial[2]) /
    3.0 * inv_volume * _unit_factor.nktv2p;

  SetParameter(PARA_PRESSURE, _pressure_scalar);
  //std::cout << " pressure=" << _pressure_scalar << std::endl;
}

void LennardJonesSystem::Compute_Temp_Scalar() {}

void LennardJonesSystem::Couple()
{
  p_current[0] = p_current[1] = p_current[2] = _pressure_scalar;
}

void LennardJonesSystem::x2lamda(Id n)
{
  //
  Vec3f delta;
  //auto position = GetFieldAsArrayHandle<Vec3f>(field::position);
  //
  vtkm::Vec<vtkm::Range, 3> range = GetParameter<vtkm::Vec<vtkm::Range, 3>>(PARA_RANGE);
  for (int i = 0; i < n; i++)
  {
    delta[0] = _position.ReadPortal().Get(i)[0] - range[0].Min;
    delta[1] = _position.ReadPortal().Get(i)[1] - range[1].Min;
    delta[2] = _position.ReadPortal().Get(i)[2] - range[2].Min;

    Vec3f lamda_position = {h_inv[0] * delta[0] + h_inv[5] * delta[1] + h_inv[4] * delta[2],
                            h_inv[1] * delta[1] + h_inv[3] * delta[2],
                            h_inv[2] * delta[2] };

    _position.WritePortal().Set(i, lamda_position);        
  }
  //_locator.SetPosition(_position);
  //for (int i = 0; i < n; i++)
  //{
  //  std::cout << "i=" << i << ", lamda_position=" << _position.ReadPortal().Get(i)[0] << ","
  //            << _position.ReadPortal().Get(i)[1] << "," << _position.ReadPortal().Get(i)[2]
  //            << std::endl;
  //}
}

void LennardJonesSystem::lamda2x(Id n)
{
  //auto position = GetFieldAsArrayHandle<Vec3f>(field::position);
  vtkm::Vec<vtkm::Range, 3> range = GetParameter<vtkm::Vec<vtkm::Range, 3>>(PARA_RANGE);
  Vec3f position_base;

  for (int i = 0; i < n; i++)
  {
    position_base[0] = _position.ReadPortal().Get(i)[0];
    position_base[1] = _position.ReadPortal().Get(i)[1];
    position_base[2] = _position.ReadPortal().Get(i)[2];

    Vec3f x_position = { h[0] * position_base[0] + h[5] * position_base[1] +
                         h[4] * position_base[2] + static_cast<vtkm::FloatDefault>(range[0].Min),
                         h[1] * position_base[1] + h[3] * position_base[2] +
                           static_cast<vtkm::FloatDefault>(range[1].Min),
                         h[2] * position_base[2] + static_cast<vtkm::FloatDefault>(range[2].Min) };

    _position.WritePortal().Set(i, x_position);
  } 

  //for (int i = 0; i < n; i++)
  //{
  //  std::cout << "i=" << i << ", new_position=" <<
  //            _position.ReadPortal().Get(i)[0] << "," << 
  //            _position.ReadPortal().Get(i)[1] << ","  << _position.ReadPortal().Get(i)[2]
  //            << std::endl;
  //}
}

void LennardJonesSystem::set_global_box()
{
  vtkm::Vec<vtkm::Range, 3> range = GetParameter<vtkm::Vec<vtkm::Range, 3>>(PARA_RANGE);
  prd[0] = xprd = range[0].Max - range[0].Min;
  prd[1] = yprd = range[1].Max - range[1].Min;
  prd[2] = zprd = range[2].Max - range[2].Min;

  h[0] = xprd;
  h[1] = yprd;
  h[2] = zprd;
  h[3] = 0;
  h[4] = 0;
  h[5] = 0;

  //
  auto orthogonal = 1;
  if (orthogonal)
  {
    h_inv[0] = 1.0 / h[0];
    h_inv[1] = 1.0 / h[1];
    h_inv[2] = 1.0 / h[2];
    h_inv[3] = 0;
    h_inv[4] = 0;
    h_inv[5] = 0;
  }
}

void LennardJonesSystem::remap()
{
  Real oldlo, oldhi, ctr;
  auto n = _position.GetNumberOfValues();

  // convert pertinent atoms and rigid bodies to lamda coords
  x2lamda(n);

  // reset global and local box to new size/shape
  vtkm::Vec<vtkm::Range, 3> range = GetParameter<vtkm::Vec<vtkm::Range, 3>>(PARA_RANGE);
  for (int i = 0; i < 3; i++)
  {
    oldlo = range[i].Min;
    oldhi = range[i].Max;
    ctr = 0.5 * (oldlo + oldhi);
    range[i].Min = (oldlo - ctr) * dilation[i] + ctr;
    range[i].Max = (oldhi - ctr) * dilation[i] + ctr;
  }
  //
  SetParameter(PARA_RANGE, range);

  set_global_box();

  // convert pertinent atoms and rigid bodies back to box coords

  lamda2x(n);
}
