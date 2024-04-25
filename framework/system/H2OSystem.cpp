#include "H2OSystem.h"
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
#include <vtkm/worklet/WorkletMapField.h>
#include "system/worklet/MolecularWorklet.h"
#include <vtkm/cont/ArrayHandleGroupVec.h>
#include <vtkm/worklet/WorkletReduceByKey.h>
#include <vtkm/worklet/Keys.h>
#include "locator/ContPointLocator.h"
#include "forceFunction/ContForceFunction.h"
//RegisterObject(H2OSystem);

template<typename T>
void PrintArrayhandle(const vtkm::cont::ArrayHandle<T> arrayHandle)
{
  auto num = arrayHandle.GetNumberOfValues();
  auto read_protol = arrayHandle.ReadPortal();
  for (int i = 0; i < num; ++i)
  {
    std::cout << read_protol.Get(i) << std::endl;
  }
}

H2OSystem::H2OSystem(const Configuration& cfg)
  : MDSystem(cfg) 
  , _RBE_P(Get<IdComponent>("rbeP"))
  , _executioner((_app.GetExecutioner()))
  , _kbT(Get<IdComponent>("kbT"))
  , _farforce_type(Get<std::string>("farforce_type"))
  , _temp_con_type(Get<std::string>("temp_con_type"))
  , _nearforce_type(Get<std::string>("nearforce_type"))
  , _use_shake(Get<bool>("use_shake"))
  , _Kmax(Get<IdComponent>("kmax"))
  , _alpha(Get<Real>("alpha"))
  , _cut_off(GetParameter<Real>(PARA_CUTOFF))
{
  SetParameter(PARA_RBE_P, _RBE_P);
  SetParameter(PARA_ALPHA, _alpha);
  SetParameter(PARA_KMAX, _Kmax);
  SetParameter(PARA_TEMPT_SUM, Real{ 0.0 });
  SetParameter(PARA_TEMPT, Real{ 0.0 });
}

void H2OSystem::Init()
{
  MDSystem::Init();

  InitialCondition();

  ComputeForce(); // Presolve force
}

void H2OSystem::PreSolve()
{
  _locator.SetPosition(_position);
  PreForce();
}

void H2OSystem::Solve()
{
  // stage1:
  //ComputeForce();   // Only compute once during evaluation
  UpdateVelocity();

  // stage2:
  UpdatePosition();

  //New added
  if (_use_shake)
  {
    ConstraintA();
  }

  // stage3:
  ComputeForce();
  UpdateVelocity();

  //New added
  if (_use_shake)
  {
    ConstraintB();
  }

  ComputeTempe();
  UpdateVelocityByTempConType();
}

void H2OSystem::PostSolve() {}

void H2OSystem::InitialCondition()
{
  MDSystem::InitialCondition();
  //_rho
  _Volume = GetParameter<Real>(PARA_VOLUME);
  SetCenterTargetPositions();
  auto target_position = GetFieldAsArrayHandle<Vec3f>(field::target_position);
  SetParameter(PARA_RDF_RHO, target_position.GetNumberOfValues() / _Volume);

  PreForce();
  _nosehooverxi = 0.0;

  set_global_box();
}

void H2OSystem::ComputeForce()
{
  ComputeAllForce();
  TempConTypeForce();
}

void H2OSystem::ComputeAllForce()
{
    // FarNearLJforce
    //SystemWorklet::SumFarNearLJForce(EleNewForce(), EleNearForce(), LJForce(), _all_force); //RBE + LJ
    SystemWorklet::SumFarNearForce(EleNewForce(), NearForce(), _all_force); //RBE + LJ

    // special coul
    Invoker{}(MolecularWorklet::AddForceWorklet{}, SpecialCoulForce(), _all_force);

    //all force with bondforce
    Invoker{}(MolecularWorklet::AddForceWorklet{}, BondForce(), _all_force);
    
    //all force with angleforce
    Invoker{}(MolecularWorklet::AddForceWorklet{}, AngleForce(), _all_force);  
}

void H2OSystem::UpdateVelocity()
{
  try
  {
    auto n = _position.GetNumberOfValues();

    _old_velocity.Allocate(n);
    for (int i = 0; i < n; i++)
    {
      _old_velocity.WritePortal().Set(i, _velocity.ReadPortal().Get(i));
    }

    SystemWorklet::UpdateVelocity(_dt, _unit_factor._fmt2v, _all_force, _mass, _velocity);
  }
  catch (const std::exception& e)
  {
    std::cout << e.what() << std::endl;
  }
}

void H2OSystem::UpdatePosition()
{
  auto n = _position.GetNumberOfValues();
  // store old position
  _old_position.Allocate(n);
  for (int i = 0; i < n; i++)
  {
    _old_position.WritePortal().Set(i, _position.ReadPortal().Get(i));
  }

  //if (_use_shake == 0)
  //{
    //std::cout << "Before UpdatePosition: position flag[46][2] = "
    //          << position_flag.ReadPortal().Get(46)[2] << std::endl;
    //std::cout << "Before UpdatePosition:_position[46][2] = " << _position.ReadPortal().Get(46)[2]
    //          << std::endl;
    //std::cout << "\n";
    SystemWorklet::UpdatePosition(_dt, _velocity, _locator, _position);

    //std::cout << "After UpdatePosition: position flag[46][2] = "
    //          << position_flag.ReadPortal().Get(46)[2] << std::endl;
    //std::cout << "After UpdatePosition:_position[46][2] = " << _position.ReadPortal().Get(46)[2]
    //          << std::endl;
    //
   
    //
    fix_press_berendsen();

    _locator.SetPosition(_position);

    SetCenterTargetPositions();
  //}
  
}

void H2OSystem::UpdateVelocityByTempConType()
{
    if (_temp_con_type == "NOSE_HOOVER")
    {
      //Nose Hoover
      //Maybe tauT = 20.0 * dt is a good choice for dt = 2e-3 from the test.
      //Because the temperature curve is the smoothest of all test simulations.
      //In fact, 5.0 and 10.0 are also optional.
      //As long as the coefficent is not too large, such as larger than 100 * dt.
      SystemWorklet::UpdateVelocityNoseHoover(
        _dt, _unit_factor._fmt2v, _nosehooverxi, _all_force, _mass, _velocity);
      Real tauT = vtkm::Pow(10.0,-1) * _dt;
      _nosehooverxi += 0.5 * _dt * (_tempT / _kbT - 1.0) / tauT;
    }
    else if(_temp_con_type ==  "TEMP_RESCALE")
    {
      Real coeff_rescale = vtkm::Sqrt(_kbT / _tempT);
      SystemWorklet::UpdateVelocityRescale(coeff_rescale, _velocity);
    }
    else if (_temp_con_type == "BERENDSEN")
    {
      //
      //Velocity Rescale: Berendsen
      //Maybe dt_divide_taut = 0.05 is a good choice for dt = 2e-3. 0.005, 0.01, 0.1 is optional.
      //The selection of dt_divide_taut determines the temperature equilibrium time.
      //
      Real dt_divide_taut = 0.1;
      Real coeff_Berendsen = vtkm::Sqrt(1.0 + dt_divide_taut * (_kbT / _tempT - 1.0));
      SystemWorklet::UpdateVelocityRescale(coeff_Berendsen, _velocity);
    }
}

void H2OSystem::SetCenterTargetPositions() 
{
  auto atom_id_center = GetFieldAsArrayHandle<Id>(field::atom_id_center);
  auto atom_id_target = GetFieldAsArrayHandle<Id>(field::atom_id_target);
  auto center_position = GetFieldAsArrayHandle<Vec3f>(field::center_position);
  auto target_position = GetFieldAsArrayHandle<Vec3f>(field::target_position);

  Invoker{}(
    MolecularWorklet::GetPositionByTypeWorklet{}, atom_id_center, _position, center_position);
  Invoker{}(
    MolecularWorklet::GetPositionByTypeWorklet{}, atom_id_target, _position, target_position);
}

void H2OSystem::PreForce()
{
  _Vlength = GetParameter<Real>(PARA_VLENGTH);
  _dt = _executioner->Dt();
  // prepare for RBE force
  auto velocity_type = GetParameter<std::string>(gtest::velocity_type);
  auto random = (velocity_type != "TEST") ? true : false;
  RBEPSAMPLE rbe_presolve_psample = { _alpha, _Vlength, _RBE_P};
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

vtkm::cont::ArrayHandle<Vec3f> H2OSystem::LJForce()
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

vtkm::cont::ArrayHandle<Vec3f> H2OSystem::EleNearForce()
{
  if (_use_erf == true)
  {
    SystemWorklet::ComputeNearElectrostaticsERF(_atoms_id, _static_table, _locator, _topology, _force_function, _ele_near_force);
  }
  else
  {
    SystemWorklet::ComputeNearElectrostatics(_atoms_id, _locator, _topology, _force_function, _ele_near_force);
  }
  Invoker{}(MolecularWorklet::UnitRescaleWorklet{ _unit_factor._qqr2e }, _ele_near_force);
  return _ele_near_force;
}

vtkm::cont::ArrayHandle<Vec3f> H2OSystem::NearForce()
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

vtkm::cont::ArrayHandle<Vec3f> H2OSystem::BondForce()
{
  // original
  auto bondlist = GetFieldAsArrayHandle<Id>(field::bond_atom_id);
  auto bond_type = GetFieldAsArrayHandle<Id>(field::bond_type);
  vtkm::IdComponent bondlist_num = bondlist.GetNumberOfValues();
  auto&& bondlist_group = vtkm::cont::make_ArrayHandleGroupVec<2>(bondlist);
  auto bond_num = bondlist_group.GetNumberOfValues();


  // bond_coeffs_k   bond_coeffs_equilibrium
  auto bond_coeffs_k = GetFieldAsArrayHandle<Real>(field::bond_coeffs_k);
  auto bond_coeffs_equilibrium = GetFieldAsArrayHandle<Real>(field::bond_coeffs_equilibrium);

  // forcebond
  vtkm::cont::ArrayHandle<Vec3f> forcebond;
  vtkm::cont::ArrayHandle<Real> bond_energy;
  auto&& forcebond_group = vtkm::cont::make_ArrayHandleGroupVec<2>(forcebond);

  Invoker{}(MolecularWorklet::ComputeBondHarmonicWorklet{ _Vlength },
            bond_type,
            bondlist_group,
            bond_coeffs_k,
            bond_coeffs_equilibrium,
            _position,
            forcebond_group,
            bond_energy,
            _locator);

  auto bond_energy_avr =
    vtkm::cont::Algorithm::Reduce(bond_energy, vtkm::TypeTraits<Real>::ZeroInitialization()) /
    bond_num;
  SetParameter(PARA_BOND_ENERGY, bond_energy_avr);

  //reduce bond force
  vtkm::worklet::Keys<vtkm::Id> keys_bond(bondlist);
  vtkm::cont::ArrayHandle<Vec3f> reduce_force_bond;
  Invoker{}(MolecularWorklet::ReduceForceWorklet{}, keys_bond, forcebond, reduce_force_bond);
  return reduce_force_bond;
}

vtkm::cont::ArrayHandle<Vec3f> H2OSystem::AngleForce()
{
  //angle
  auto angle_list = GetFieldAsArrayHandle<Id>(field::angle_atom_id);
  auto angle_type = GetFieldAsArrayHandle<Id>(field::angle_type);
  vtkm::IdComponent anglelist_num = angle_list.GetNumberOfValues();
  auto&& anglelist_group = vtkm::cont::make_ArrayHandleGroupVec<3>(angle_list);
  auto angle_num = anglelist_group.GetNumberOfValues();

  // angle_coeffs_k   angle_coeffs_equilibrium
  auto angle_coeffs_k = GetFieldAsArrayHandle<Real>(field::angle_coeffs_k);
  auto angle_coeffs_equilibrium = GetFieldAsArrayHandle<Real>(field::angle_coeffs_equilibrium);

  // force_angle
  vtkm::cont::ArrayHandle<Vec3f> force_angle;
  vtkm::cont::ArrayHandle<Real> angle_energy;
  auto&& forceangle_group = vtkm::cont::make_ArrayHandleGroupVec<3>(force_angle);

  Invoker{}(MolecularWorklet::ComputeAngleHarmonicWorklet{ _Vlength },
            angle_type,
            anglelist_group,
            angle_coeffs_k,
            angle_coeffs_equilibrium,
            _position,
            forceangle_group,
            angle_energy,
            _locator);

  auto angle_energy_avr =
    vtkm::cont::Algorithm::Reduce(angle_energy, vtkm::TypeTraits<Real>::ZeroInitialization()) /
    angle_num;
  SetParameter(PARA_ANGLE_ENERGY, angle_energy_avr);

  //reduce angle force
  vtkm::worklet::Keys<vtkm::Id> keys_angle(angle_list);
  vtkm::cont::ArrayHandle<Vec3f> reduce_force_angle;
  Invoker{}(MolecularWorklet::ReduceForceWorklet{}, keys_angle, force_angle, reduce_force_angle);
  return reduce_force_angle;
}

vtkm::cont::ArrayHandle<Vec3f> H2OSystem::SpecialCoulForce()
{ 
  auto vlength = GetParameter<Real>(PARA_VLENGTH);
  auto source_array = GetFieldAsArrayHandle<Id>(field::special_source_array);
  auto offsets_array = GetFieldAsArrayHandle<Id>(field::special_offsets_array);
  auto groupVecArray = vtkm::cont::make_ArrayHandleGroupVecVariable(source_array, offsets_array);
  SystemWorklet::ComputeSpecialCoul(vlength, _atoms_id, groupVecArray, _force_function, _topology, _locator, _spec_coul_force);

  return _spec_coul_force;
}

vtkm::cont::ArrayHandle<Vec3f> H2OSystem::EleNewForce()
{
  //_EleFartimer.Start();
  if(_farforce_type == "RBE")
  {
    // New RBE force part
    ComputeRBEEleForce(_psample, _RBE_P, _ele_new_force);
    //_Elefartimer_counting = _Elefartimer_counting + _EleFartimer.GetElapsedTime();
    //std::cout << "RBE time: " << _Elefartimer_counting << std::endl;
  }
  else if(_farforce_type == "EWALD")
  {
    // New Ewald far part
    ComputeEwaldEleForce(_Kmax,_ele_new_force);
    //_Elefartimer_counting = _Elefartimer_counting + _EleFartimer.GetElapsedTime();
    //std::cout << "Ewald time: " << _Elefartimer_counting << std::endl;
  }
  //_EleFartimer.Stop();

  Invoker{}(MolecularWorklet::UnitRescaleWorklet{ _unit_factor._qqr2e }, _ele_new_force);
  return _ele_new_force;
}

void H2OSystem::TempConTypeForce()
{
  vtkm::cont::ArrayHandle<Real> mass;
  mass.Allocate(_all_force.GetNumberOfValues());
  for (size_t i = 0; i < _all_force.GetNumberOfValues(); i++)
  {
    mass.WritePortal().Set(i, 1);
  }
  if (_temp_con_type == "LANGEVIN")
  {
    // Underdamped Langevin
    Real kBT = 1.0;
    Real gamma = 100.0;
    SystemWorklet::UnderdampedLangevin(_gaussian, kBT, gamma, _dt, mass, _velocity, _all_force);
  }
}

void H2OSystem::ComputeTempe()
{
  auto n = _position.GetNumberOfValues();
  ArrayHandle<Real> sq_velocity;
  sq_velocity.Allocate(n);

  SystemWorklet::ComputerKineticEnergy(_velocity, _mass, sq_velocity);
  Invoker{}(MolecularWorklet::UnitRescaleWorklet{ _unit_factor._mvv2e }, sq_velocity);
  _tempT_sum =
    vtkm::cont::Algorithm::Reduce(sq_velocity, vtkm::TypeTraits<Real>::ZeroInitialization());

  //////////////////////////////////////////////////////////////////////
  //dof_shake is important!!!
  //It is used only for shake option, which may be modified in the future.
  //When shake is called, dof_shake = n(number of atoms of water); otherwise, dof_shake = 0
  // 3 is the extra dof
  ///////////////////////////////////////////////////////////////////////
  IdComponent dof_shake = _use_shake ? n : 0;
  auto field = _use_shake ? 3 : 0;
  //int field = 0;
  Real temperature_kB = _unit_factor._kB;
  _tempT = 0.5 * _tempT_sum / ((3 * n - dof_shake - field) * temperature_kB / 2.0);

  SetParameter(PARA_TEMPT_SUM, _tempT_sum);
  SetParameter(PARA_TEMPT, _tempT);
}

void H2OSystem::SetForceFunction()
{
  InitERF();
  auto cut_off = GetParameter<Real>(PARA_CUTOFF);
  auto alpha = GetParameter<Real>(PARA_ALPHA);
  auto volume = GetParameter<Real>(PARA_VOLUME);
  auto vlength = GetParameter<Real>(PARA_VLENGTH);
  auto Kmax = GetParameter<IdComponent>(PARA_KMAX);
  auto rbe_p = GetParameter<IdComponent>(PARA_RBE_P);

  _force_function.SetParameters(cut_off, alpha, volume, vlength, Kmax);
}

void H2OSystem::SetTopology()
{
  auto pts_type = GetFieldAsArrayHandle<Id>(field::pts_type);
  auto molecule_id = GetFieldAsArrayHandle<Id>(field::molecule_id);
  auto epsilon = GetFieldAsArrayHandle<Real>(field::epsilon);
  auto sigma = GetFieldAsArrayHandle<Real>(field::sigma);
  _topology.SetAtomsType(pts_type);
  _topology.SetMolecularId(molecule_id);
  _topology.SetEpsAndSigma(epsilon, sigma);

  auto source_array = GetFieldAsArrayHandle<Id>(field::special_source_array);
  auto offsets_array = GetFieldAsArrayHandle<Id>(field::special_offsets_array);
  _topology.SetSourceAndOffsets(source_array, offsets_array);
}

void H2OSystem::InitField() 
{
  MDSystem::InitField();
  AddField(field::position_flag, ArrayHandle<Id3>{});
  AddField(field::bond_atom_id, ArrayHandle<Id>{});
  AddField(field::bond_type, ArrayHandle<Id>{});
  AddField(field::bond_coeffs_k, ArrayHandle<Real>{});
  AddField(field::bond_coeffs_equilibrium, ArrayHandle<Real>{});
  AddField(field::angle_atom_id, ArrayHandle<Id>{});
  AddField(field::angle_type, ArrayHandle<Id>{});
  AddField(field::angle_coeffs_k, ArrayHandle<Real>{});
  AddField(field::angle_coeffs_equilibrium, ArrayHandle<Real>{});
  AddField(field::atom_id_center, ArrayHandle<Id>{});
  AddField(field::atom_id_target, ArrayHandle<Id>{});
  AddField(field::pts_type, ArrayHandle<Id>{});
  AddField(field::center_position, ArrayHandle<Vec3f>{});
  AddField(field::target_position, ArrayHandle<Vec3f>{});
  AddField(field::epsilon, ArrayHandle<Real>{});
  AddField(field::sigma, ArrayHandle<Real>{});
  AddField(field::dihedrals_atom_id, ArrayHandle<Id>{});
  AddField(field::dihedrals_type, ArrayHandle<Id>{});
  AddField(field::signal_atoms_id, ArrayHandle<Id>{});
  AddField(field::special_source_array, ArrayHandle<Id>{});
  AddField(field::special_offsets_array, ArrayHandle<Id>{});
}

void H2OSystem::TimeIntegration() {}

void H2OSystem::Rattle()
{
  auto angle_list = GetFieldAsArrayHandle<Id>(field::angle_atom_id);
  auto&& anglelist_group = vtkm::cont::make_ArrayHandleGroupVec<3>(angle_list);

  auto data_range = GetParameter<vtkm::Vec<vtkm::Range, 3>>(PARA_RANGE);
  Vec<Vec2f, 3> range{
    { static_cast<Real>(data_range[0].Min), static_cast<Real>(data_range[0].Max) },
    { static_cast<Real>(data_range[1].Min), static_cast<Real>(data_range[1].Max) },
    { static_cast<Real>(data_range[2].Min), static_cast<Real>(data_range[2].Max) }
  };

   Invoker{}(
    MolecularWorklet::ConstraintWaterVelocityBondAngleWorklet{ _Vlength, _dt, _unit_factor._fmt2v, range },
    anglelist_group,
    _position,
    _velocity,
    _all_force,
    _mass,
    _locator);
}

void H2OSystem::ConstraintA()
{
  auto angle_list = GetFieldAsArrayHandle<Id>(field::angle_atom_id);
  auto&& anglelist_group = vtkm::cont::make_ArrayHandleGroupVec<3>(angle_list);

  auto data_range = GetParameter<vtkm::Vec<vtkm::Range, 3>>(PARA_RANGE);
  Vec<Vec2f, 3> range{
    { static_cast<Real>(data_range[0].Min), static_cast<Real>(data_range[0].Max) },
    { static_cast<Real>(data_range[1].Min), static_cast<Real>(data_range[1].Max) },
    { static_cast<Real>(data_range[2].Min), static_cast<Real>(data_range[2].Max) }
  };
  
  auto&& position_flag = GetFieldAsArrayHandle<Id3>(field::position_flag);

  Invoker{}(
    MolecularWorklet::NewConstraintAWaterBondAngleWorklet{ _Vlength, _dt, _unit_factor._fmt2v, range },
    anglelist_group,
    _old_position,
    _old_velocity,
    _all_force,
    _mass,
    position_flag,
    _locator);

  for(int i = 0; i < _position.GetNumberOfValues(); i++)
  {
    _position.WritePortal().Set(i, _old_position.ReadPortal().Get(i));
    _velocity.WritePortal().Set(i, _old_velocity.ReadPortal().Get(i));
  }

  /* auto&& position_flag = GetFieldAsArrayHandle<Id3>(field::position_flag);
  auto&& output_flag = _app.GetExecutioner()->CurrentStep() >= _flag_steps ? 1 : 0;

  std::cout << "Before UpdatePosition: position flag[272][2] = "
            << position_flag.ReadPortal().Get(272)[2] << std::endl;
  std::cout << "Before UpdatePosition:_position[272][2] = " << _position.ReadPortal().Get(272)[2]
            << std::endl;
  std::cout << "\n";
  SystemWorklet::UpdatePosition(_dt, output_flag, _velocity, _locator, _position, position_flag);

  std::cout << "After UpdatePosition: position flag[272][2] = "
            << position_flag.ReadPortal().Get(272)[2] << std::endl;
  std::cout << "After UpdatePosition:_position[272][2] = " << _position.ReadPortal().Get(272)[2]
            << std::endl;

  _locator.UpdatePosition(_position);
  _force_function.SetPtsPosition(_position);

  SetCenterTargetPositions();*/
}

void H2OSystem::ConstraintB()
{
  auto angle_list = GetFieldAsArrayHandle<Id>(field::angle_atom_id);
  auto&& anglelist_group = vtkm::cont::make_ArrayHandleGroupVec<3>(angle_list);

  auto data_range = GetParameter<vtkm::Vec<vtkm::Range, 3>>(PARA_RANGE);
  Vec<Vec2f, 3> range{
    { static_cast<Real>(data_range[0].Min), static_cast<Real>(data_range[0].Max) },
    { static_cast<Real>(data_range[1].Min), static_cast<Real>(data_range[1].Max) },
    { static_cast<Real>(data_range[2].Min), static_cast<Real>(data_range[2].Max) }
  };

  Invoker{}(
    MolecularWorklet::NewConstraintBWaterBondAngleWorklet{ _Vlength, _dt, _unit_factor._fmt2v, range },
    anglelist_group,
    _position,
    _old_velocity,
    _all_force,
    _mass,
    _locator);


  for(int i = 0; i < _velocity.GetNumberOfValues(); i++)
  {
    _velocity.WritePortal().Set(i , _old_velocity.ReadPortal().Get(i));
  }
}

void H2OSystem::fix_press_berendsen()
{
  bulkmodulus = 10.0;
  p_start[0] = p_start[1] = p_start[2] = 1.0;
  p_stop[0] = p_stop[1] = p_stop[2] = 1.0;
  p_period[0] = p_period[1] = p_period[2] = 500;

  // compute new T,P

  Compute_Pressure_Scalar();
  Couple();

  auto currentstep = _app.GetExecutioner()->CurrentStep();
  auto beginstep = 0;
  auto endstep = _app.GetExecutioner()->NumStep();

  auto delta = currentstep - beginstep;
  std::cout << delta << std::endl;
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

void H2OSystem::ComputeVirial()
{
  auto cut_off = GetParameter<Real>(PARA_CUTOFF);
  SystemWorklet::LJVirial(
    cut_off, _atoms_id, _locator, _topology, _force_function, _virial_atom);


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

void H2OSystem::Compute_Pressure_Scalar()
{
  //compute temperature
  auto temperature = GetParameter<Real>(PARA_TEMPT);

  // compute  virial
  vtkm::Vec<vtkm::Range, 3> range = GetParameter<vtkm::Vec<vtkm::Range, 3>>(PARA_RANGE);
  auto volume =  (range[0].Max - range[0].Min) * 
                 (range[1].Max - range[1].Min) * 
                 (range[2].Max - range[2].Min);

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
  std::cout << " pressure=" << _pressure_scalar << std::endl;
}

void H2OSystem::Compute_Temp_Scalar() {}

void H2OSystem::Couple()
{
  p_current[0] = p_current[1] = p_current[2] = _pressure_scalar;
}

void H2OSystem::x2lamda(Id n)
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

    _position.WritePortal().Get(i)[0] =
      h_inv[0] * delta[0] + h_inv[5] * delta[1] + h_inv[4] * delta[2];
    _position.WritePortal().Get(i)[1] = h_inv[1] * delta[1] + h_inv[3] * delta[2];
    _position.WritePortal().Get(i)[2] = h_inv[2] * delta[2];
  }
  _locator.SetPosition(_position);
  //for (int i = 0; i < n; i++)
  //{
  //  for (int j = 0; j < 3; j++)
  //  {
  //    std::cout << "i=" << i << ", lamda_position=" << position.ReadPortal().Get(i)[j] << std::endl;
  //  }
  //}
}

void H2OSystem::lamda2x(Id n)
{
  //auto position = GetFieldAsArrayHandle<Vec3f>(field::position);
  vtkm::Vec<vtkm::Range, 3> range = GetParameter<vtkm::Vec<vtkm::Range, 3>>(PARA_RANGE);

  for (int i = 0; i < n; i++)
  {
    _position.WritePortal().Get(i)[0] = h[0] * _position.ReadPortal().Get(i)[0] +
      h[5] * _position.ReadPortal().Get(i)[1] + h[4] * _position.ReadPortal().Get(i)[2] +
      range[0].Min;

    _position.WritePortal().Get(i)[1] = h[1] * _position.ReadPortal().Get(i)[1] +
      h[3] * _position.ReadPortal().Get(i)[2] + range[1].Min;

    _position.WritePortal().Get(i)[2] = h[2] * _position.ReadPortal().Get(i)[2] + range[2].Min;
  }
  _locator.SetPosition(_position);

  //for (int i = 0; i < n; i++)
  //{
  //  for (int j = 0; j < 3; j++)
  //  {
  //    std::cout << "i=" << i << ", position=" << position.ReadPortal().Get(i)[j] << std::endl;
  //  }
  //}
}

void H2OSystem::set_global_box()
{
  vtkm::Vec<vtkm::Range, 3> range = GetParameter<vtkm::Vec<vtkm::Range, 3>>(PARA_RANGE);
  prd[0] = xprd = range[0].Max - range[0].Min;
  prd[1] = yprd = range[1].Max - range[1].Min;
  prd[2] = zprd = range[2].Max - range[2].Min;

  h[0] = xprd;
  h[1] = yprd;
  h[2] = zprd;

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

void H2OSystem::remap()
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