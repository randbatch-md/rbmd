#include "PCFFSystem.h"
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
#include <iomanip>

template<typename T>
void PrintArrayhandle(const vtkm::cont::ArrayHandle<T> arrayHandle)
{
  auto num = arrayHandle.GetNumberOfValues();
  auto read_protol = arrayHandle.ReadPortal();
  for (int i = 0; i < num; ++i)
  {
    std::cout << std::right << std::setw(10) << read_protol.Get(i)[0] << " " << std::right 
              << std::setw(10) << read_protol.Get(i)[1] << " " << std::right << std::setw(10)
              << read_protol.Get(i)[2]
              << std::endl;
  }
}

PCFFSystem::PCFFSystem(const Configuration& cfg)
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

void PCFFSystem::Init()
{
  MDSystem::Init();

  InitialCondition();

  ComputeForce(); // Presolve force
}

void PCFFSystem::PreSolve()
{
  _locator.SetPosition(_position);
  PreForce();
}

void PCFFSystem::Solve()
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

void PCFFSystem::PostSolve() {}

void PCFFSystem::InitialCondition()
{
  MDSystem::InitialCondition();
  //_rho
  _Volume = GetParameter<Real>(PARA_VOLUME);
  SetCenterTargetPositions();
  auto target_position = GetFieldAsArrayHandle<Vec3f>(field::target_position);
  SetParameter(PARA_RDF_RHO, target_position.GetNumberOfValues() / _Volume);

  PreForce();
  _nosehooverxi = 0.0;
}

void PCFFSystem::ComputeForce()
{
  ComputeAllForce();
  TempConTypeForce();
  //PrintArrayhandle(_all_force);
  //std::cout << "------------------------------------------------------------" << std::endl;
  //PrintArrayhandle(_position);
}

void PCFFSystem::ComputeAllForce()
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

void PCFFSystem::UpdateVelocity()
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

void PCFFSystem::UpdatePosition()
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

    _locator.SetPosition(_position);

    SetCenterTargetPositions();
  //}
  
}

void PCFFSystem::UpdateVelocityByTempConType()
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
      Real dt_divide_taut = 0.5;
      Real coeff_Berendsen = vtkm::Sqrt(1.0 + dt_divide_taut * (_kbT / _tempT - 1.0));
      SystemWorklet::UpdateVelocityRescale(coeff_Berendsen, _velocity);
    }
}

void PCFFSystem::SetCenterTargetPositions() 
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

void PCFFSystem::PreForce()
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

vtkm::cont::ArrayHandle<Vec3f> PCFFSystem::LJForce()
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

vtkm::cont::ArrayHandle<Vec3f> PCFFSystem::EleNearForce()
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

vtkm::cont::ArrayHandle<Vec3f> PCFFSystem::NearForce()
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

vtkm::cont::ArrayHandle<Vec3f> PCFFSystem::BondForce()
{
  // original
  auto bondlist = GetFieldAsArrayHandle<Id>(field::bond_atom_id);
  auto bond_type = GetFieldAsArrayHandle<Id>(field::bond_type);
  vtkm::IdComponent bondlist_num = bondlist.GetNumberOfValues();
  auto&& bondlist_group = vtkm::cont::make_ArrayHandleGroupVec<2>(bondlist);
  auto bond_num = bondlist_group.GetNumberOfValues();


  // bond_coeffs_k   bond_coeffs_equilibrium
  auto bond_coeffs_k = GetFieldAsArrayHandle<Real>(field::bond_coeffs_k);
  auto bond_coeffs_k3 = GetFieldAsArrayHandle<Real>(field::bond_coeffs_k3);
  auto bond_coeffs_k4 = GetFieldAsArrayHandle<Real>(field::bond_coeffs_k4);
  auto bond_coeffs_equilibrium = GetFieldAsArrayHandle<Real>(field::bond_coeffs_equilibrium);

  // forcebond
  vtkm::cont::ArrayHandle<Vec3f> forcebond;
  vtkm::cont::ArrayHandle<Real> bond_energy;
  auto&& forcebond_group = vtkm::cont::make_ArrayHandleGroupVec<2>(forcebond);
  Invoker{}(MolecularWorklet::ComputeBondClass2Worklet{ _Vlength },
            bond_type,
            bondlist_group,
            bond_coeffs_k,
            bond_coeffs_k3,
            bond_coeffs_k4,
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

vtkm::cont::ArrayHandle<Vec3f> PCFFSystem::AngleForce()
{
  //angle
  auto angle_list = GetFieldAsArrayHandle<Id>(field::angle_atom_id);
  auto angle_type = GetFieldAsArrayHandle<Id>(field::angle_type);
  vtkm::IdComponent anglelist_num = angle_list.GetNumberOfValues();
  auto&& anglelist_group = vtkm::cont::make_ArrayHandleGroupVec<3>(angle_list);
  auto angle_num = anglelist_group.GetNumberOfValues();

  // angle_coeffs_k   angle_coeffs_equilibrium
  auto angle_coeffs_k2 = GetFieldAsArrayHandle<Real>(field::angle_coeffs_k);
  auto angle_coeffs_k3 = GetFieldAsArrayHandle<Real>(field::angle_coeffs_k3);
  auto angle_coeffs_k4 = GetFieldAsArrayHandle<Real>(field::angle_coeffs_k4);
  auto angle_coeffs_equilibrium = GetFieldAsArrayHandle<Real>(field::angle_coeffs_equilibrium);
  auto bb_m  = GetFieldAsArrayHandle<Real>(field::bb_m);
  auto bb_r1 = GetFieldAsArrayHandle<Real>(field::bb_r1);
  auto bb_r2 = GetFieldAsArrayHandle<Real>(field::bb_r2);
  auto ba_n1 = GetFieldAsArrayHandle<Real>(field::ba_n1);
  auto ba_n2 = GetFieldAsArrayHandle<Real>(field::ba_n2);
  auto ba_r1 = GetFieldAsArrayHandle<Real>(field::ba_r1);
  auto ba_r2 = GetFieldAsArrayHandle<Real>(field::ba_r2);

  // force_angle
  vtkm::cont::ArrayHandle<Vec3f> force_angle;
  vtkm::cont::ArrayHandle<Real> angle_energy;
  auto&& forceangle_group = vtkm::cont::make_ArrayHandleGroupVec<3>(force_angle);

  Invoker{}(MolecularWorklet::ComputeAngleClass2Worklet{ _Vlength },
            angle_type,
            anglelist_group,
            angle_coeffs_k2,
            angle_coeffs_k3,
            angle_coeffs_k4,
            angle_coeffs_equilibrium,
            bb_m,
            bb_r1,
            bb_r2,
            ba_n1,
            ba_n2,
            ba_r1,
            ba_r2,
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

vtkm::cont::ArrayHandle<Vec3f> PCFFSystem::SpecialCoulForce()
{ 
  auto vlength = GetParameter<Real>(PARA_VLENGTH);
  auto source_array = GetFieldAsArrayHandle<Id>(field::special_source_array);
  auto offsets_array = GetFieldAsArrayHandle<Id>(field::special_offsets_array);
  auto groupVecArray = vtkm::cont::make_ArrayHandleGroupVecVariable(source_array, offsets_array);
  SystemWorklet::ComputeSpecialCoul(vlength, _atoms_id, groupVecArray, _force_function, _topology, _locator, _spec_coul_force);

  return _spec_coul_force;
}

vtkm::cont::ArrayHandle<Vec3f> PCFFSystem::EleNewForce()
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

  //Invoker{}(MolecularWorklet::UnitRescaleWorklet{ _unit_factor._qqr2e }, _ele_new_force);
  //vtkm::cont::ArrayHandle<Vec3f> arrayHandle;
  //arrayHandle.Allocate(_ele_new_force.GetNumberOfValues());

  //// 初始化数组中的每个元素为 (0, 0, 0)
  //Vec3f initValue(0.0f, 0.0f, 0.0f);
  //auto portal = arrayHandle.WritePortal();
  //for (vtkm::Id i = 0; i < _ele_new_force.GetNumberOfValues(); ++i)
  //{
  //  portal.Set(i, initValue);
  //}

  return _ele_new_force;
}

void PCFFSystem::TempConTypeForce()
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

void PCFFSystem::ComputeTempe()
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

void PCFFSystem::SetForceFunction()
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

void PCFFSystem::SetTopology()
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

void PCFFSystem::InitField() 
{
  MDSystem::InitField();
  AddField(field::position_flag, ArrayHandle<Id3>{});
  AddField(field::bond_atom_id, ArrayHandle<Id>{});
  AddField(field::bond_type, ArrayHandle<Id>{});
  AddField(field::bond_coeffs_k, ArrayHandle<Real>{});
  AddField(field::bond_coeffs_k3, ArrayHandle<Real>{});
  AddField(field::bond_coeffs_k4, ArrayHandle<Real>{});
  AddField(field::bond_coeffs_equilibrium, ArrayHandle<Real>{});
  AddField(field::angle_atom_id, ArrayHandle<Id>{});
  AddField(field::angle_type, ArrayHandle<Id>{});
  AddField(field::angle_coeffs_k, ArrayHandle<Real>{});
  AddField(field::angle_coeffs_k3, ArrayHandle<Real>{});
  AddField(field::angle_coeffs_k4, ArrayHandle<Real>{});
  AddField(field::angle_coeffs_equilibrium, ArrayHandle<Real>{});
  AddField(field::bb_m, ArrayHandle<Real>{});
  AddField(field::bb_r1, ArrayHandle<Real>{});
  AddField(field::bb_r2, ArrayHandle<Real>{});
  AddField(field::ba_n1, ArrayHandle<Real>{});
  AddField(field::ba_n2, ArrayHandle<Real>{});
  AddField(field::ba_r1, ArrayHandle<Real>{});
  AddField(field::ba_r2, ArrayHandle<Real>{});
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

void PCFFSystem::TimeIntegration() {}

void PCFFSystem::Rattle()
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

void PCFFSystem::ConstraintA()
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

void PCFFSystem::ConstraintB()
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