#include "ExecutionNVT.h"
#include "Application.h"
#include "FieldName.h"
#include "MeshFreeCondition.h"
#include "math/Math.h"
#include "math/Utiles.h"
#include "RBEPSample.h"
#include "run/worklet/MolecularWorklet.h"
#include "run/worklet/RunWorklet.h"
#include <cmath> // erfc(x)
#include <fstream>
#include <vtkm/Math.h>
#include <vtkm/Pair.h>
#include <vtkm/VectorAnalysis.h>
#include <vtkm/cont/Algorithm.h>
#include <vtkm/cont/ArrayHandleGroupVec.h>
#include <vtkm/cont/ArrayHandleGroupVecVariable.h>
#include <vtkm/cont/ConvertNumComponentsToOffsets.h>
#include <vtkm/worklet/Keys.h>
#include <vtkm/worklet/WorkletMapField.h>
#include <vtkm/worklet/WorkletReduceByKey.h>

//RegisterObject(ExecutionNVT);

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

ExecutionNVT::ExecutionNVT(const Configuration& cfg)
  : ExecutionMD(cfg)
  , _executioner((_app.GetExecutioner()))
{
  ExecutionMD::SetParameters();
  InitParameters();
}

void ExecutionNVT::Init()
{
  ExecutionMD::Init();

  InitialCondition();

  ComputeForce(); // Presolve force
}

void ExecutionNVT::PreSolve()
{
  _locator.SetPosition(_position);
  PreForce();
}

void ExecutionNVT::Solve()
{
  // stage1:
  //ComputeForce();   // Only compute once during evaluation
  UpdateVelocity();

  // stage2:
  UpdatePosition();

  //New added
  if (_para.GetParameter<std::string>(PARA_FIX_SHAKE) == "true")
  {
    ConstraintA();
  }

  // stage3:
  ComputeForce();
  UpdateVelocity();

  //New added
  if (_para.GetParameter<std::string>(PARA_FIX_SHAKE) == "true")
  {
    ConstraintB();
  }

  ComputeTempe();
  UpdateVelocityByTempConType();
}

void ExecutionNVT::PostSolve() {}

void ExecutionNVT::InitialCondition()
{
  ExecutionMD::InitialCondition();
  //_rho
  _Volume = _para.GetParameter<Real>(PARA_VOLUME);
  SetCenterTargetPositions();
  auto target_position = _para.GetFieldAsArrayHandle<Vec3f>(field::target_position);
  _para.SetParameter(PARA_RDF_RHO, target_position.GetNumberOfValues() / _Volume);

  PreForce();
  _nosehooverxi = 0.0;
}

void ExecutionNVT::ComputeForce()
{
  ComputeAllForce();
  TempConTypeForce();
}

void ExecutionNVT::ComputeAllForce()
{
  // FarNearLJforce
  //RunWorklet::SumFarNearLJForce(EleNewForce(), EleNearForce(), LJForce(), _all_force); //RBE + LJ
  RunWorklet::SumFarNearForce(EleNewForce(), NearForce(), _all_force); //RBE + LJ
  // special coul
  Invoker{}(MolecularWorklet::AddForceWorklet{}, SpecialCoulForce(), _all_force);

  //all force with bondforce
  Invoker{}(MolecularWorklet::AddForceWorklet{}, BondForce(), _all_force);

  //all force with angleforce
  Invoker{}(MolecularWorklet::AddForceWorklet{}, AngleForce(), _all_force);

  if (_para.GetParameter<bool>(PARA_DIHEDRALS_FORCE))
  {
    //all force with dihedral+force
    Invoker{}(MolecularWorklet::AddForceWorklet{}, DihedralsForce(), _all_force);
  }
}

void ExecutionNVT::UpdateVelocity()
{
  try
  {
    //auto n = _position.GetNumberOfValues();
    //
    //_old_velocity.Allocate(n);
    //for (int i = 0; i < n; i++)
    //{
    //  _old_velocity.WritePortal().Set(i, _velocity.ReadPortal().Get(i));
    //}
    vtkm::cont::ArrayCopy(_velocity, _old_velocity);

    RunWorklet::UpdateVelocity(_dt, _unit_factor._fmt2v, _all_force, _mass, _velocity);
  }
  catch (const std::exception& e)
  {
    std::cout << e.what() << std::endl;
  }
}

void ExecutionNVT::UpdatePosition()
{
  //auto n = _position.GetNumberOfValues();
  // store old position
  //_old_position.Allocate(n);
  //for (int i = 0; i < n; i++)
  //{
  //  _old_position.WritePortal().Set(i, _position.ReadPortal().Get(i));
  //}

  vtkm::cont::ArrayCopy(_position, _old_position);

  if (_para.GetParameter<std::string>(PARA_FIX_SHAKE) == "null")
  {
    auto&& position_flag = _para.GetFieldAsArrayHandle<Id3>(field::position_flag);
    RunWorklet::UpdatePositionFlag(_dt, _velocity, _locator, _position, position_flag);
  }
  else
  {
    RunWorklet::UpdatePosition(_dt, _velocity, _locator, _position);
  }

  _locator.SetPosition(_position);
  SetCenterTargetPositions();
}

void ExecutionNVT::UpdateVelocityByTempConType()
{
    // tauT 和 dt_divide_taut 与 temperature[3] 相关 目前暂定统一为一个常量
  _temp_con_type = _para.GetParameter<std::string>(PARA_TEMP_CTRL_TYPE);
  if (_temp_con_type == "NOSE_HOOVER")
  {
    //Nose Hoover
    //Maybe tauT = 20.0 * dt is a good choice for dt = 2e-3 from the test.
    //Because the temperature curve is the smoothest of all test simulations.
    //In fact, 5.0 and 10.0 are also optional.
    //As long as the coefficent is not too large, such as larger than 100 * dt.
    RunWorklet::UpdateVelocityNoseHoover(
      _dt, _unit_factor._fmt2v, _nosehooverxi, _all_force, _mass, _velocity);
    //Real tauT = 20.0 * _dt;
    Real tauT = vtkm::Pow(10.0, -1) * _dt;
    _nosehooverxi += 0.5 * _dt * (_tempT / _kbT - 1.0) / tauT;
  }
  else if (_temp_con_type == "TEMP_RESCALE")
  {
    Real coeff_rescale = vtkm::Sqrt(_kbT / _tempT);
    RunWorklet::UpdateVelocityRescale(coeff_rescale, _velocity);
  }
  else if (_temp_con_type == "BERENDSEN")
  {
    //
    //Velocity Rescale: Berendsen
    //Maybe dt_divide_taut = 0.05 is a good choice for dt = 2e-3. 0.005, 0.01, 0.1 is optional.
    //The selection of dt_divide_taut determines the temperature equilibrium time.
    
    //Real dt_divide_taut = 0.02;
    Real dt_divide_taut = 0.1;
    Real coeff_Berendsen = vtkm::Sqrt(1.0 + dt_divide_taut * (_kbT / _tempT - 1.0));
    RunWorklet::UpdateVelocityRescale(coeff_Berendsen, _velocity);
  }
}

void ExecutionNVT::SetCenterTargetPositions()
{
  auto atom_id_center = _para.GetFieldAsArrayHandle<Id>(field::atom_id_center);
  auto atom_id_target = _para.GetFieldAsArrayHandle<Id>(field::atom_id_target);
  auto center_position = _para.GetFieldAsArrayHandle<Vec3f>(field::center_position);
  auto target_position = _para.GetFieldAsArrayHandle<Vec3f>(field::target_position);

  Invoker{}(
    MolecularWorklet::GetPositionByTypeWorklet{}, atom_id_center, _position, center_position);
  Invoker{}(
    MolecularWorklet::GetPositionByTypeWorklet{}, atom_id_target, _position, target_position);
}

void ExecutionNVT::PreForce()
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

vtkm::cont::ArrayHandle<Vec3f> ExecutionNVT::LJForce()
{
  _nearforce_type = _para.GetParameter<std::string>(PARA_NEIGHBOR_TYPE);
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
    //ComputeOriginalLJForce(_LJforce);
    ComputeSpecialBondsLJForce(_LJforce);
  }
  return _LJforce;
}

vtkm::cont::ArrayHandle<Vec3f> ExecutionNVT::EleNearForce()
{
  if (_use_erf == true)
  {
    RunWorklet::ComputeNearElectrostaticsERF(
      _atoms_id, _static_table, _locator, _topology, _force_function, _ele_near_force);
  }
  else
  {
    RunWorklet::ComputeNearElectrostatics(
      _atoms_id, _locator, _topology, _force_function, _ele_near_force);
  }
  Invoker{}(MolecularWorklet::UnitRescaleWorklet{ _unit_factor._qqr2e }, _ele_near_force);
  return _ele_near_force;
}

vtkm::cont::ArrayHandle<Vec3f> ExecutionNVT::NearForce()
{
  _nearforce_type = _para.GetParameter<std::string>(PARA_NEIGHBOR_TYPE);

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
    RunWorklet::SumFarNearForce(EleNearForce(), LJForce(), _nearforce);
  }
  return _nearforce;
}

vtkm::cont::ArrayHandle<Vec3f> ExecutionNVT::BondForce()
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
  if (_para.GetParameter<bool>(PARA_DIHEDRALS_FORCE))
  {
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
    memcpy(
      &new_forcebond[0], forcebond.ReadPortal().GetArray(), original_number * sizeof(vtkm::Vec3f));

    std::vector<Vec3f> value(atom_id_number, vtkm::Vec3f{ 0.0, 0.0, 0.0 });
    new_forcebond.insert(new_forcebond.end(), value.begin(), value.end());

    vtkm::worklet::Keys<vtkm::Id> keys_bond(vtkm::cont::make_ArrayHandle(new_bondlist));
    Invoker{}(MolecularWorklet::ReduceForceWorklet{},
              keys_bond,
              vtkm::cont::make_ArrayHandle(new_forcebond),
              reduce_force_bond);
  }
  else
  {
    vtkm::worklet::Keys<vtkm::Id> keys_bond(bondlist);
    Invoker{}(MolecularWorklet::ReduceForceWorklet{}, keys_bond, forcebond, reduce_force_bond);
  }
  return reduce_force_bond;
}

vtkm::cont::ArrayHandle<Vec3f> ExecutionNVT::AngleForce()
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
  if (_para.GetParameter<bool>(PARA_DIHEDRALS_FORCE))
  {
    //reduce angle force
    auto atom_id = _para.GetFieldAsArrayHandle<Id>(field::atom_id);
    auto atom_id_number = atom_id.GetNumberOfValues();
    auto original_number = angle_list.GetNumberOfValues();

    std::vector<Id> new_anglelist(original_number);
    memcpy(
      &new_anglelist[0], angle_list.ReadPortal().GetArray(), original_number * sizeof(vtkm::Id));

    std::vector<vtkm::Id> append_list(atom_id_number);
    memcpy(&append_list[0], atom_id.ReadPortal().GetArray(), atom_id_number * sizeof(vtkm::Id));

    //append key
    new_anglelist.insert(new_anglelist.end(), append_list.begin(), append_list.end());

    //append force bond
    std::vector<vtkm::Vec3f> new_force_angle(original_number);
    memcpy(&new_force_angle[0],
           force_angle.ReadPortal().GetArray(),
           original_number * sizeof(vtkm::Vec3f));
    std::vector<Vec3f> value(atom_id_number, vtkm::Vec3f{ 0.0, 0.0, 0.0 });
    new_force_angle.insert(new_force_angle.end(), value.begin(), value.end());

    vtkm::worklet::Keys<vtkm::Id> keys_angle(vtkm::cont::make_ArrayHandle(new_anglelist));
    Invoker{}(MolecularWorklet::ReduceForceWorklet{},
              keys_angle,
              vtkm::cont::make_ArrayHandle(new_force_angle),
              reduce_force_angle);
  }
  else
  {
    vtkm::worklet::Keys<vtkm::Id> keys_angle(angle_list);
    Invoker{}(MolecularWorklet::ReduceForceWorklet{}, keys_angle, force_angle, reduce_force_angle);
  }

  return reduce_force_angle;
}

vtkm::cont::ArrayHandle<Vec3f> ExecutionNVT::DihedralsForce()
{
  //dihedrals
  auto dihedrals_list = _para.GetFieldAsArrayHandle<Id>(field::dihedrals_atom_id);
  auto dihedrals_type = _para.GetFieldAsArrayHandle<Id>(field::dihedrals_type);
  vtkm::IdComponent dihedralslist_num = dihedrals_list.GetNumberOfValues();
  auto&& dihedralslist_group = vtkm::cont::make_ArrayHandleGroupVec<4>(dihedrals_list);
  auto dihedrals_num = dihedralslist_group.GetNumberOfValues();

  // dihedrals_coeffs_k   dihedrals_coeffs_sign     dihedrals_coeffs_multiplicity
  auto dihedrals_coeffs_k = _para.GetFieldAsArrayHandle<Real>(field::dihedrals_coeffs_k);
  auto dihedrals_coeffs_sign =
    _para.GetFieldAsArrayHandle<vtkm::IdComponent>(field::dihedrals_coeffs_sign);
  auto dihedrals_coeffs_multiplicity =
    _para.GetFieldAsArrayHandle<vtkm::IdComponent>(field::dihedrals_coeffs_multiplicity);

  // force_dihedrals
  vtkm::cont::ArrayHandle<Vec3f> force_dihedrals;
  vtkm::cont::ArrayHandle<Real> dihedrals_energy;
  auto&& forcedihedrals_group = vtkm::cont::make_ArrayHandleGroupVec<4>(force_dihedrals);

  //auto a1 = dihedrals_type.GetNumberOfValues();
  //auto a2 = dihedralslist_group.GetNumberOfValues();
  //auto a3 = dihedrals_coeffs_k.GetNumberOfValues();
  //auto a4 = dihedrals_coeffs_sign.GetNumberOfValues();
  //auto a5 = dihedrals_coeffs_multiplicity.GetNumberOfValues();
  Invoker{}(MolecularWorklet::ComputeDihedralHarmonicWorklet{ _Vlength },
            dihedrals_type,
            dihedralslist_group,
            dihedrals_coeffs_k,
            dihedrals_coeffs_sign,
            dihedrals_coeffs_multiplicity,
            _position,
            forcedihedrals_group,
            dihedrals_energy,
            _locator);

  auto dihedrals_energy_avr =
    vtkm::cont::Algorithm::Reduce(dihedrals_energy, vtkm::TypeTraits<Real>::ZeroInitialization()) /
    dihedrals_num;
  _para.SetParameter(PARA_DIHEDRAL_ENERGY, dihedrals_energy_avr);

  vtkm::cont::ArrayHandle<Vec3f> reduce_force_dihedrals;
  //reduce dihedrals force
  auto atom_id = _para.GetFieldAsArrayHandle<Id>(field::atom_id);
  auto atom_id_number = atom_id.GetNumberOfValues();
  auto original_number = dihedrals_list.GetNumberOfValues();

  std::vector<Id> new_dihedralslist(original_number);
  memcpy(&new_dihedralslist[0],
         dihedrals_list.ReadPortal().GetArray(),
         original_number * sizeof(vtkm::Id));

  std::vector<vtkm::Id> append_list(atom_id_number);
  memcpy(&append_list[0], atom_id.ReadPortal().GetArray(), atom_id_number * sizeof(vtkm::Id));

  //append key
  new_dihedralslist.insert(new_dihedralslist.end(), append_list.begin(), append_list.end());

  //append force bond
  std::vector<vtkm::Vec3f> new_force_dihedrals(original_number);
  memcpy(&new_force_dihedrals[0],
         force_dihedrals.ReadPortal().GetArray(),
         original_number * sizeof(vtkm::Vec3f));
  std::vector<Vec3f> value(atom_id_number, vtkm::Vec3f{ 0.0, 0.0, 0.0 });
  new_force_dihedrals.insert(new_force_dihedrals.end(), value.begin(), value.end());

  vtkm::worklet::Keys<vtkm::Id> keys_dihedrals(vtkm::cont::make_ArrayHandle(new_dihedralslist));
  Invoker{}(MolecularWorklet::ReduceForceWorklet{},
            keys_dihedrals,
            vtkm::cont::make_ArrayHandle(new_force_dihedrals),
            reduce_force_dihedrals);

  return reduce_force_dihedrals;
}

vtkm::cont::ArrayHandle<Vec3f> ExecutionNVT::SpecialCoulForce()
{
  auto vLength = _para.GetParameter<Real>(PARA_VLENGTH);
  auto source_array = _para.GetFieldAsArrayHandle<Id>(field::special_source_array);
  auto offsets_array = _para.GetFieldAsArrayHandle<Id>(field::special_offsets_array);
  auto groupVecArray = vtkm::cont::make_ArrayHandleGroupVecVariable(source_array, offsets_array);

  if (_para.GetParameter<bool>(PARA_DIHEDRALS_FORCE))
  {
    auto special_offsets = _para.GetFieldAsArrayHandle<Id>(field::special_offsets);
    auto special_weights = _para.GetFieldAsArrayHandle<Real>(field::special_weights);
    auto specoal_ids = _para.GetFieldAsArrayHandle<Id>(field::special_ids);
    auto ids_group = vtkm::cont::make_ArrayHandleGroupVecVariable(specoal_ids, special_offsets);
    auto weight_group =
      vtkm::cont::make_ArrayHandleGroupVecVariable(special_weights, special_offsets);

    RunWorklet::ComputeSpecialCoulGeneral(vLength,
                                             _atoms_id,
                                             groupVecArray,
                                             _force_function,
                                             _topology,
                                             _locator,
                                             ids_group,
                                             weight_group,
                                             _spec_coul_force);
  }
  else
  {  
    RunWorklet::ComputeSpecialCoul(vLength, _atoms_id, groupVecArray, _force_function, _topology, _locator, _spec_coul_force);
  }
  return _spec_coul_force;
}

vtkm::cont::ArrayHandle<Vec3f> ExecutionNVT::EleNewForce()
{
  _farforce_type = _para.GetParameter<std::string>(PARA_COULOMB_TYPE);
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

void ExecutionNVT::TempConTypeForce()
{
  vtkm::cont::ArrayHandle<Real> mass;
  mass.AllocateAndFill(_all_force.GetNumberOfValues(), 1.0);
  _temp_con_type = _para.GetParameter<std::string>(PARA_TEMP_CTRL_TYPE);
  if (_temp_con_type == "LANGEVIN")
  {
    // Underdamped Langevin
    Real kBT = 1.0;
    Real gamma = 100.0;
    RunWorklet::UnderdampedLangevin(_gaussian, kBT, gamma, _dt, mass, _velocity, _all_force);
  }
}

void ExecutionNVT::ComputeTempe()
{
  auto n = _position.GetNumberOfValues();
  ArrayHandle<Real> sq_velocity;
  sq_velocity.Allocate(n);

  RunWorklet::ComputerKineticEnergy(_velocity, _mass, sq_velocity);
  Invoker{}(MolecularWorklet::UnitRescaleWorklet{ _unit_factor._mvv2e }, sq_velocity);
  _tempT_sum =
    vtkm::cont::Algorithm::Reduce(sq_velocity, vtkm::TypeTraits<Real>::ZeroInitialization());

  //////////////////////////////////////////////////////////////////////
  //dof_shake is important!!!
  //It is used only for shake option, which may be modified in the future.
  //When shake is called, dof_shake = n(number of atoms of water); otherwise, dof_shake = 0
  // 3 is the extra dof
  ///////////////////////////////////////////////////////////////////////
  /*IdComponent dof_shake = _use_shake ? n : 0;
  auto field = _use_shake ? 3 : 0;*/
  auto shake = _para.GetParameter<std::string>(PARA_FIX_SHAKE);
  Real temperature_kB = _unit_factor._kB;
  if (shake == "null")
  {
    _tempT = 0.5 * _tempT_sum / ((3 * n - 3) * temperature_kB / 2.0);
  }
  else if (shake == "false")
  {
    _tempT = 0.5 * _tempT_sum / ((3 * n ) * temperature_kB / 2.0);
  }
  else
  {
    _tempT = 0.5 * _tempT_sum / ((3 * n - n - 3) * temperature_kB / 2.0);
  }
  _para.SetParameter(PARA_TEMPT_SUM, _tempT_sum);
  _para.SetParameter(PARA_TEMPT, _tempT);
}

void ExecutionNVT::SetForceFunction()
{
  InitERF();
  auto cut_off = _para.GetParameter<Real>(PARA_CUTOFF);
  auto alpha = _para.GetParameter<Real>(PARA_ALPHA);
  auto volume = _para.GetParameter<Real>(PARA_VOLUME);
  auto vlength = _para.GetParameter<Real>(PARA_VLENGTH);
  auto Kmax = _para.GetParameter<IdComponent>(PARA_KMAX);

  _force_function.SetParameters(cut_off, alpha, volume, vlength, Kmax);
}

void ExecutionNVT::SetTopology()
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

void ExecutionNVT::InitParameters()
{
  ExecutionMD::InitParameters();
  _RBE_P = _para.GetParameter<IdComponent>(PARA_COULOMB_SAMPLE_NUM);
  _kbT = _para.GetParameter<std::vector<Real>>(PARA_TEMPERATURE)[0]; 
  _para.SetParameter(PARA_TEMPT_SUM, Real{ 0.0 });
  _para.SetParameter(PARA_TEMPT, Real{ 0.0 });
  _alpha = _para.GetParameter<Real>(PARA_ALPHA);
}

//void ExecutionNVT::InitField()
//{
//  ExecutionMD::InitField();
//  _para.AddField(field::position_flag, ArrayHandle<Id3>{});
//  _para.AddField(field::bond_atom_id, ArrayHandle<Id>{});
//  _para.AddField(field::bond_type, ArrayHandle<Id>{});
//  _para.AddField(field::bond_coeffs_k, ArrayHandle<Real>{});
//  _para.AddField(field::bond_coeffs_equilibrium, ArrayHandle<Real>{});
//  _para.AddField(field::angle_atom_id, ArrayHandle<Id>{});
//  _para.AddField(field::angle_type, ArrayHandle<Id>{});
//  _para.AddField(field::angle_coeffs_k, ArrayHandle<Real>{});
//  _para.AddField(field::angle_coeffs_equilibrium, ArrayHandle<Real>{});
//  _para.AddField(field::atom_id_center, ArrayHandle<Id>{});
//  _para.AddField(field::atom_id_target, ArrayHandle<Id>{});
//  _para.AddField(field::pts_type, ArrayHandle<Id>{});
//  _para.AddField(field::center_position, ArrayHandle<Vec3f>{});
//  _para.AddField(field::target_position, ArrayHandle<Vec3f>{});
//  _para.AddField(field::epsilon, ArrayHandle<Real>{});
//  _para.AddField(field::sigma, ArrayHandle<Real>{});
//  _para.AddField(field::signal_atoms_id, ArrayHandle<Id>{});
//  _para.AddField(field::special_source_array, ArrayHandle<Id>{});
//  _para.AddField(field::special_offsets_array, ArrayHandle<Id>{});
//
//  _para.AddField(field::dihedrals_atom_id, ArrayHandle<Id>{});
//  _para.AddField(field::dihedrals_type, ArrayHandle<Id>{});
//  _para.AddField(field::dihedrals_coeffs_k, ArrayHandle<Real>{});
//  _para.AddField(field::dihedrals_coeffs_sign, ArrayHandle<vtkm::IdComponent>{});
//  _para.AddField(field::dihedrals_coeffs_multiplicity, ArrayHandle<vtkm::IdComponent>{});
//
//  _para.AddField(field::position_flag, ArrayHandle<Id3>{});
//}

void ExecutionNVT::TimeIntegration() {}

void ExecutionNVT::ConstraintA()
{
  //vtkm::cont::Timer timer4ConstraintA;
  //timer4ConstraintA.Start();
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
    MolecularWorklet::NewConstraintAWaterBondAngleWorklet{
      _Vlength, _dt, _unit_factor._fmt2v, range },
    anglelist_group,
    _old_position,
    _old_velocity,
    _all_force,
    _mass,
    position_flag,
    _locator);

  vtkm::cont::ArrayCopy(_old_position, _position);
  vtkm::cont::ArrayCopy(_old_velocity, _velocity);
}

void ExecutionNVT::ConstraintB()
{
  //vtkm::cont::Timer timer4ConstraintB;
  //timer4ConstraintB.Start();
  auto angle_list = _para.GetFieldAsArrayHandle<Id>(field::angle_atom_id);
  auto&& anglelist_group = vtkm::cont::make_ArrayHandleGroupVec<3>(angle_list);

  auto data_range = _para.GetParameter<vtkm::Vec<vtkm::Range, 3>>(PARA_RANGE);
  Vec<Vec2f, 3> range{
    { static_cast<Real>(data_range[0].Min), static_cast<Real>(data_range[0].Max) },
    { static_cast<Real>(data_range[1].Min), static_cast<Real>(data_range[1].Max) },
    { static_cast<Real>(data_range[2].Min), static_cast<Real>(data_range[2].Max) }
  };

  Invoker{}(
    MolecularWorklet::NewConstraintBWaterBondAngleWorklet{
      _Vlength, _dt, _unit_factor._fmt2v, range },
    anglelist_group,
    _position,
    _old_velocity,
    _all_force,
    _mass,
    _locator);



  vtkm::cont::ArrayCopy(_old_velocity, _velocity);
}