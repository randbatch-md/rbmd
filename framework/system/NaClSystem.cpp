#include "NaClSystem.h"
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

//RegisterObject(NaClSystem);

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

NaClSystem::NaClSystem(const Configuration& cfg)
  : MDSystem(cfg) 
  , _RBE_P(Get<IdComponent>("rbeP"))
  , _executioner((_app.GetExecutioner()))
  , _kbT(Get<IdComponent>("kbT"))
  , _farforce_type(Get<std::string>("farforce_type"))
  , _nearforce_type(Get<std::string>("nearforce_type"))
  , _temp_con_type(Get<std::string>("temp_con_type"))
  , _use_shake(Get<bool>("use_shake"))
  , _Kmax(Get<IdComponent>("kmax"))
  , _alpha(Get<Real>("alpha"))
  , _cut_off(_para.GetParameter<Real>(PARA_CUTOFF))
{
  _para.SetParameter(PARA_RBE_P, _RBE_P);
  _para.SetParameter(PARA_ALPHA, _alpha);
  _para.SetParameter(PARA_KMAX, _Kmax);
  _para.SetParameter(PARA_TEMPT_SUM, Real{ 0.0 });
  _para.SetParameter(PARA_TEMPT, Real{ 0.0 });
}

void NaClSystem::Init()
{
  MDSystem::Init();

  InitialCondition();

  ComputeForce(); // Presolve force
}

void NaClSystem::PreSolve()
{
  _locator.SetPosition(_position);
  PreForce();
}

void NaClSystem::Solve()
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
  //auto transient = static_cast<Transient&>(*(_app.GetExecutioner()));
  //if(transient.CurrentStep() > 4000)
  //{
  //  Rattle();
  //}

  ComputeTempe();
  UpdateVelocityByTempConType();
}

void NaClSystem::PostSolve() {}

void NaClSystem::InitialCondition()
{
  MDSystem::InitialCondition();
  //_rho
  _Volume = _para.GetParameter<Real>(PARA_VOLUME);
  SetCenterTargetPositions();
  auto target_position = _para.GetFieldAsArrayHandle<Vec3f>(field::target_position);
  _para.SetParameter(PARA_RDF_RHO, target_position.GetNumberOfValues() / _Volume);

  PreForce();
  _nosehooverxi = 0.0;
}

void NaClSystem::ComputeForce()
{
  ComputeAllForce();
  TempConTypeForce();
}

void NaClSystem::ComputeAllForce()
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

void NaClSystem::UpdateVelocity()
{
  try
  {
    auto n = _position.GetNumberOfValues();

    //_old_velocity.Allocate(n);
    //for (int i = 0; i < n; i++)
    //{
    //  _old_velocity.WritePortal().Set(i, _velocity.ReadPortal().Get(i));
    //}
    vtkm::cont::ArrayCopy(_velocity, _old_velocity);
    SystemWorklet::UpdateVelocity(_dt, _unit_factor._fmt2v, _all_force, _mass, _velocity);
  }
  catch (const std::exception& e)
  {
    std::cout << e.what() << std::endl;
  }
}

void NaClSystem::UpdatePosition()
{
  auto n = _position.GetNumberOfValues();
  // store old position
  //_old_position.Allocate(n);
  //for (int i = 0; i < n; i++)
  //{
  //  _old_position.WritePortal().Set(i, _position.ReadPortal().Get(i));
  //}
  vtkm::cont::ArrayCopy(_position, _old_position);
  SystemWorklet::UpdatePosition(_dt, _velocity, _locator, _position);

  _locator.SetPosition(_position);
  SetCenterTargetPositions();
}

void NaClSystem::UpdateVelocityByTempConType()
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
    else if(_temp_con_type ==  "TEMP_RESCALE")
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
      Real dt_divide_taut = 0.1;
      Real coeff_Berendsen = vtkm::Sqrt(1.0 + dt_divide_taut * (_kbT / _tempT - 1.0));
      SystemWorklet::UpdateVelocityRescale(coeff_Berendsen, _velocity);
    }
}

void NaClSystem::SetCenterTargetPositions() 
{
  auto atom_id_center = _para.GetFieldAsArrayHandle<Id>(field::atom_id_center );
  auto atom_id_target = _para.GetFieldAsArrayHandle<Id>(field::atom_id_target);
  auto center_position = _para.GetFieldAsArrayHandle<Vec3f>(field::center_position);
  auto target_position = _para.GetFieldAsArrayHandle<Vec3f>(field::target_position);

  Invoker{}(
    MolecularWorklet::GetPositionByTypeWorklet{}, atom_id_center, _position, center_position);
  Invoker{}(
    MolecularWorklet::GetPositionByTypeWorklet{}, atom_id_target, _position, target_position);
}

void NaClSystem::PreForce()
{
  _Vlength = _para.GetParameter<Real>(PARA_VLENGTH);
  _dt = _executioner->Dt();
  // prepare for RBE force
  auto velocity_type = _para.GetParameter<std::string>(gtest::velocity_type);
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

vtkm::cont::ArrayHandle<Vec3f> NaClSystem::LJForce()
{

  if (_nearforce_type == "RBL")
  {
    ComputeRBLLJForce(_LJforce);
  }
  if (_nearforce_type == "VERLETLIST")
  {
    ComputeVerletlistLJForce(_LJforce);
  }
  else if (_nearforce_type == "ORIGINAL")
  {
    ComputeOriginalLJForce(_LJforce);
  }
  return _LJforce;
}

vtkm::cont::ArrayHandle<Vec3f> NaClSystem::EleNearForce()
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

vtkm::cont::ArrayHandle<Vec3f> NaClSystem::NearForce()
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

vtkm::cont::ArrayHandle<Vec3f> NaClSystem::BondForce()
{
  // original
  auto bondlist = _para.GetFieldAsArrayHandle<Id>(field::bond_atom_id);
  auto bond_type = _para.GetFieldAsArrayHandle<Id>(field::bond_type);
  vtkm::IdComponent bondlist_num = bondlist.GetNumberOfValues();
  auto&& bondlist_group = vtkm::cont::make_ArrayHandleGroupVec<2>(bondlist);
  auto bond_num = bondlist_group.GetNumberOfValues();


  // bond_coeffs_k   bond_coeffs_equilibrium
  auto bond_coeffs_k = _para.GetFieldAsArrayHandle<Real>(field::bond_coeffs_k);
  auto bond_coeffs_equilibrium = _para.GetFieldAsArrayHandle<Real>(field::bond_coeffs_equilibrium);

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
  _para.SetParameter(PARA_BOND_ENERGY, bond_energy_avr);

  vtkm::cont::ArrayHandle<Vec3f> reduce_force_bond;
  //reduce bond force
  auto atom_id = _para.GetFieldAsArrayHandle<Id>(field::atom_id);
  auto atom_id_number = atom_id.GetNumberOfValues();
  auto original_number = bondlist.GetNumberOfValues();

  std::vector<Id> new_bondlist(original_number);
  memcpy(&new_bondlist[0], bondlist.ReadPortal().GetArray(), original_number * sizeof(vtkm::Id));

  std::vector<vtkm::Id> append_list(atom_id_number);
  memcpy(&append_list[0], atom_id.ReadPortal().GetArray(), atom_id_number * sizeof(vtkm::Id));

  //append key
  new_bondlist.insert(new_bondlist.end(), append_list.begin(), append_list.end());

  //append force bond
  std::vector<vtkm::Vec3f> new_forcebond(original_number);
  memcpy(&new_forcebond[0], forcebond.ReadPortal().GetArray(), original_number * sizeof(vtkm::Vec3f));

  std::vector<Vec3f> value(atom_id_number, vtkm::Vec3f{ 0.0, 0.0, 0.0 });
  new_forcebond.insert(new_forcebond.end(), value.begin(), value.end());

  vtkm::worklet::Keys<vtkm::Id> keys_bond(vtkm::cont::make_ArrayHandle(new_bondlist));
  Invoker{}(MolecularWorklet::ReduceForceWorklet{},
            keys_bond,
            vtkm::cont::make_ArrayHandle(new_forcebond),
            reduce_force_bond);

  return reduce_force_bond;
}

vtkm::cont::ArrayHandle<Vec3f> NaClSystem::AngleForce()
{
  //angle
  auto angle_list = _para.GetFieldAsArrayHandle<Id>(field::angle_atom_id);
  auto angle_type = _para.GetFieldAsArrayHandle<Id>(field::angle_type);
  vtkm::IdComponent anglelist_num = angle_list.GetNumberOfValues();
  auto&& anglelist_group = vtkm::cont::make_ArrayHandleGroupVec<3>(angle_list);
  auto angle_num = anglelist_group.GetNumberOfValues();

  // angle_coeffs_k   angle_coeffs_equilibrium
  auto angle_coeffs_k = _para.GetFieldAsArrayHandle<Real>(field::angle_coeffs_k);
  auto angle_coeffs_equilibrium = _para.GetFieldAsArrayHandle<Real>(field::angle_coeffs_equilibrium);

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
  _para.SetParameter(PARA_ANGLE_ENERGY, angle_energy_avr);

  vtkm::cont::ArrayHandle<Vec3f> reduce_force_angle;
  //reduce angle force
  auto atom_id = _para.GetFieldAsArrayHandle<Id>(field::atom_id);
  auto atom_id_number = atom_id.GetNumberOfValues();
  auto original_number = angle_list.GetNumberOfValues();

  std::vector<Id> new_anglelist(original_number);
  memcpy(&new_anglelist[0], angle_list.ReadPortal().GetArray(), original_number * sizeof(vtkm::Id));

  std::vector<vtkm::Id> append_list(atom_id_number);
  memcpy(&append_list[0], atom_id.ReadPortal().GetArray(), atom_id_number * sizeof(vtkm::Id));

  //append key
  new_anglelist.insert(new_anglelist.end(), append_list.begin(), append_list.end());

  //append force bond
  std::vector<vtkm::Vec3f> new_force_angle(original_number);
  memcpy(&new_force_angle[0], force_angle.ReadPortal().GetArray(),original_number * sizeof(vtkm::Vec3f));
  std::vector<Vec3f> value(atom_id_number, vtkm::Vec3f{ 0.0, 0.0, 0.0 });
  new_force_angle.insert(new_force_angle.end(), value.begin(), value.end());

  vtkm::worklet::Keys<vtkm::Id> keys_angle(vtkm::cont::make_ArrayHandle(new_anglelist));
  Invoker{}(MolecularWorklet::ReduceForceWorklet{},
            keys_angle,
            vtkm::cont::make_ArrayHandle(new_force_angle),
            reduce_force_angle);

  return reduce_force_angle;
}

vtkm::cont::ArrayHandle<Vec3f> NaClSystem::SpecialCoulForce()
{
  auto vLength = _para.GetParameter<Real>(PARA_VLENGTH);
  auto source_array = _para.GetFieldAsArrayHandle<Id>(field::special_source_array);
  auto offsets_array = _para.GetFieldAsArrayHandle<Id>(field::special_offsets_array);
  auto groupVecArray = vtkm::cont::make_ArrayHandleGroupVecVariable(source_array, offsets_array);

  SystemWorklet::ComputeSpecialCoul(vLength, _atoms_id, groupVecArray, _force_function, _topology, _locator, _spec_coul_force);
  return _spec_coul_force;
}

vtkm::cont::ArrayHandle<Vec3f> NaClSystem::EleNewForce()
{
  if (_farforce_type == "RBE")
  {
    // New RBE force part
    ComputeRBEEleForce(_psample, _RBE_P, _ele_new_force);
  }
  else if (_farforce_type == "EWALD")
  {
    // New Ewald far part
    ComputeEwaldEleForce(_Kmax, _ele_new_force);
  }

  Invoker{}(MolecularWorklet::UnitRescaleWorklet{ _unit_factor._qqr2e }, _ele_new_force);
  return _ele_new_force;
}

void NaClSystem::TempConTypeForce()
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
    SystemWorklet::UnderdampedLangevin(_gaussian, kBT, gamma, _dt, mass, _velocity, _all_force);
  }
}

void NaClSystem::ComputeTempe()
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
  Real temperature_kB = _unit_factor._kB;
  _tempT = 0.5 * _tempT_sum / ((3 * n - dof_shake - field) * temperature_kB / 2.0);

  _para.SetParameter(PARA_TEMPT_SUM, _tempT_sum);
  _para.SetParameter(PARA_TEMPT, _tempT);
}

void NaClSystem::SetForceFunction()
{
  InitERF();
  auto cut_off = _para.GetParameter<Real>(PARA_CUTOFF);
  auto alpha = _para.GetParameter<Real>(PARA_ALPHA);
  auto volume = _para.GetParameter<Real>(PARA_VOLUME);
  auto vlength = _para.GetParameter<Real>(PARA_VLENGTH);
  auto Kmax = _para.GetParameter<IdComponent>(PARA_KMAX);

  _force_function.SetParameters(cut_off, alpha, volume, vlength, Kmax);
}

void NaClSystem::SetTopology()
{
  auto pts_type = _para.GetFieldAsArrayHandle<Id>(field::pts_type);
  auto molecule_id = _para.GetFieldAsArrayHandle<Id>(field::molecule_id);
  auto epsilon = _para.GetFieldAsArrayHandle<Real>(field::epsilon);
  auto sigma = _para.GetFieldAsArrayHandle<Real>(field::sigma);

  _topology.SetAtomsType(pts_type);
  _topology.SetMolecularId(molecule_id);
  _topology.SetEpsAndSigma(epsilon, sigma);

  auto source_array = _para.GetFieldAsArrayHandle<Id>(field::special_source_array);
  auto offsets_array = _para.GetFieldAsArrayHandle<Id>(field::special_offsets_array);
  _topology.SetSourceAndOffsets(source_array, offsets_array);
}

void NaClSystem::InitField() 
{
  MDSystem::InitField();
  _para.AddField(field::position_flag, ArrayHandle<Id3>{});
  _para.AddField(field::bond_atom_id, ArrayHandle<Id>{});
  _para.AddField(field::bond_type, ArrayHandle<Id>{});
  _para.AddField(field::bond_coeffs_k, ArrayHandle<Real>{});
  _para.AddField(field::bond_coeffs_equilibrium, ArrayHandle<Real>{});
  _para.AddField(field::angle_atom_id, ArrayHandle<Id>{});
  _para.AddField(field::angle_type, ArrayHandle<Id>{});
  _para.AddField(field::angle_coeffs_k, ArrayHandle<Real>{});
  _para.AddField(field::angle_coeffs_equilibrium, ArrayHandle<Real>{});
  _para.AddField(field::atom_id_center , ArrayHandle<Id>{});
  _para.AddField(field::atom_id_target, ArrayHandle<Id>{});
  _para.AddField(field::pts_type , ArrayHandle<Id>{});
  _para.AddField(field::center_position, ArrayHandle<Vec3f>{});
  _para.AddField(field::target_position, ArrayHandle<Vec3f>{});
  _para.AddField(field::epsilon, ArrayHandle<Real>{});
  _para.AddField(field::sigma, ArrayHandle<Real>{});
  _para.AddField(field::signal_atoms_id, ArrayHandle<Id>{});
  _para.AddField(field::special_source_array, ArrayHandle<Id>{});
  _para.AddField(field::special_offsets_array, ArrayHandle<Id>{});
}

void NaClSystem::TimeIntegration() {}

void NaClSystem::Rattle()
{
  auto angle_list = _para.GetFieldAsArrayHandle<Id>(field::angle_atom_id);
  auto&& anglelist_group = vtkm::cont::make_ArrayHandleGroupVec<3>(angle_list);

  auto data_range = _para.GetParameter<vtkm::Vec<vtkm::Range, 3>>(PARA_RANGE);
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

void NaClSystem::ConstraintA()
{
  auto angle_list = _para.GetFieldAsArrayHandle<Id>(field::angle_atom_id);
  auto&& anglelist_group = vtkm::cont::make_ArrayHandleGroupVec<3>(angle_list);

  auto data_range = _para.GetParameter<vtkm::Vec<vtkm::Range, 3>>(PARA_RANGE);
  Vec<Vec2f, 3> range{
    { static_cast<Real>(data_range[0].Min), static_cast<Real>(data_range[0].Max) },
    { static_cast<Real>(data_range[1].Min), static_cast<Real>(data_range[1].Max) },
    { static_cast<Real>(data_range[2].Min), static_cast<Real>(data_range[2].Max) }
  };

  auto&& position_flag = _para.GetFieldAsArrayHandle<Id3>(field::position_flag);
  Invoker{}(
    MolecularWorklet::NewConstraintAWaterBondAngleWorklet{ _Vlength, _dt, _unit_factor._fmt2v, range },
    anglelist_group,
    _old_position,
    _old_velocity,
    _all_force,
    _mass,
    position_flag,
    _locator);

  vtkm::cont::ArrayCopy(_old_position, _position);
  vtkm::cont::ArrayCopy(_old_velocity, _velocity);

  //for(int i = 0; i < _position.GetNumberOfValues(); i++)
  //{
  //  _position.WritePortal().Set(i, _old_position.ReadPortal().Get(i));
  //  _velocity.WritePortal().Set(i, _old_velocity.ReadPortal().Get(i));
  //}
}

void NaClSystem::ConstraintB()
{
  auto angle_list = _para.GetFieldAsArrayHandle<Id>(field::angle_atom_id);
  auto&& anglelist_group = vtkm::cont::make_ArrayHandleGroupVec<3>(angle_list);

  auto data_range = _para.GetParameter<vtkm::Vec<vtkm::Range, 3>>(PARA_RANGE);
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

  //for(int i = 0; i < _velocity.GetNumberOfValues(); i++)
  //{
  //  _velocity.WritePortal().Set(i , _old_velocity.ReadPortal().Get(i));
  //}
  vtkm::cont::ArrayCopy(_old_velocity, _velocity);
}