#include "ExecutionNVT.h"
#include "Application.h"
#include "FieldName.h"
#include "MeshFreeCondition.h"
#include "math/Math.h"
#include "math/Utiles.h"
#include "RBEPSample.h"
#include "run/worklet/MolecularWorklet.h"
#include "run/worklet/RunWorklet.h"
#include <cmath> 
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

#define MIN(A, B) ((A) < (B) ? (A) : (B))
#define MAX(A, B) ((A) > (B) ? (A) : (B))

ExecutionNVT::ExecutionNVT(const Configuration& cfg)
  : ExecutionMD(cfg)
  , _executioner((_app.GetExecutioner()))
{
  ExecutionMD::SetParameters();
  InitParameters();
}

void ExecutionNVT::Init()
{
  if (_para.GetParameter<std::string>(PARA_FILE_TYPE) == "EAM")
  {
    _potential_file.open(_para.GetParameter<std::string>(PARA_POTENTIAL_FILE));
    ReadPotentialFile(_potential_file);
    InitStyle();
  }

  ExecutionMD::Init();

  InitialCondition();

  ComputeForce(); 
}

void ExecutionNVT::PreSolve()
{
  _locator.SetPosition(_position);
  if (_para.GetParameter<bool>(PARA_FAR_FORCE))
  {
    PreForce();
  }
  else
  {
    std::shared_ptr<Executioner>& executioner = _app.GetExecutioner();
    _dt = executioner->Dt();
  }
}

void ExecutionNVT::Solve()
{
  // stage1:
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

  if (_para.GetParameter<std::string>(PARA_PRESS_CTRL_TYPE) == "BERENDSEN")
  {
    fix_press_berendsen();
    //fix_press_berendsen_scale();
  }
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

  if (_para.GetParameter<bool>(PARA_FAR_FORCE))
  {
    PreForce();
  }
  _nosehooverxi = 0.0;

    //
  //bulkmodulus = 100.0;
  p_start[0] = p_start[1] = p_start[2] = _Pstart;
  p_stop[0] = p_stop[1] = p_stop[2] = _Pstop;
  p_period[0] = p_period[1] = p_period[2] = _Pperiod;
  set_global_box();
}

void ExecutionNVT::ComputeForce()
{
  ComputeAllForce();
  TempConTypeForce();
}

void ExecutionNVT::ComputeAllForce()
{
  auto force_field = _para.GetParameter<std::string>(PARA_FORCE_FIELD_TYPE);
  if ("LJ/CUT" == force_field)
  {
    NearForceLJ();
  }

  else if ("LJ/CUT/COUL/LONG" == force_field)
  {
    RunWorklet::SumFarNearForce(EleNewForce(), NearForce(), _all_force);
  }

  else if ("CVFF" == force_field)
  {
    RunWorklet::SumFarNearForce(EleNewForce(), NearForce(), _all_force);

    Invoker{}(MolecularWorklet::AddForceWorklet{}, SpecialCoulForce(), _all_force);

    //all force with bondforce
    Invoker{}(MolecularWorklet::AddForceWorklet{}, BondForce(), _all_force);

    //all force with angleforce
    Invoker{}(MolecularWorklet::AddForceWorklet{}, AngleForce(), _all_force);

    if (_para.GetParameter<bool>(PARA_DIHEDRALS_FORCE) &&
        _para.GetParameter<bool>(PARA_FILE_DIHEDRALS))
    {
      //all force with dihedral+force
      Invoker{}(MolecularWorklet::AddForceWorklet{}, DihedralsForce(), _all_force);
    }
  }

  else if ("EAM" == force_field)
  {
    NearForceEAM();
  }
}

void ExecutionNVT::UpdateVelocity()
{
  try
  {
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
  vtkm::cont::ArrayCopy(_position, _old_position);

  if (_para.GetParameter<std::string>(PARA_FIX_SHAKE) == "null" || _init_way == "inbuild")
  {
    auto&& position_flag = _para.GetFieldAsArrayHandle<Id3>(field::position_flag);
    RunWorklet::UpdatePositionFlag(_dt, _velocity, _locator, _position, position_flag);
  }
  else
  {
    RunWorklet::UpdatePosition(_dt, _velocity, _locator, _position);
  }
  //pbc
  ApplyPbc();
  _locator.SetPosition(_position);
  SetCenterTargetPositions();
}

void ExecutionNVT::UpdateVelocityByTempConType()
{
    // ！注意：tauT 和 dt_divide_taut 与 temperature[3] 相关 目前暂定统一为一个常量
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
    
    //Real dt_divide_taut = 0.02; for PEO
    //Real dt_divide_taut = 0.1; // 注意：不同系统相差很大 LJ
    auto dt_divide_taut = _dt / _Tdamp;
    Real coeff_Berendsen = vtkm::Sqrt(1.0 + dt_divide_taut * (_kbT / _tempT - 1.0));
    RunWorklet::UpdateVelocityRescale(coeff_Berendsen, _velocity);
  }
}

void ExecutionNVT::SetCenterTargetPositions()
{
  if (_init_way == "inbuild")
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

    auto center_position = _para.GetFieldAsArrayHandle<Vec3f>(field::center_position);
    vtkm::cont::ArrayCopy(center_position_temp, center_position);

    auto target_position = _para.GetFieldAsArrayHandle<Vec3f>(field::target_position);
    vtkm::cont::ArrayCopy(target_position_temp, target_position);
  }
  else
  {
    auto rdf_id = _para.GetFieldAsArrayHandle<Id>(field::atom_pair_id);
    auto atoms_pair_type_offsets = _para.GetFieldAsArrayHandle<Id>(field::atoms_pair_type_offsets);
    vtkm::cont::ArrayHandle<Vec3f> atom_type_position;
    atom_type_position.Allocate(rdf_id.GetNumberOfValues());
    Invoker{}(MolecularWorklet::GetPositionByTypeWorklet{}, rdf_id, _position, atom_type_position);
    auto atom_pair_position = _para.GetFieldAsArrayHandle<Vec3f>(field::atom_pair_position);
    vtkm::cont::ArrayCopy(atom_type_position, atom_pair_position);
  }
 }

void ExecutionNVT::PreForce()
{
  _Vlength = _para.GetParameter<Real>(PARA_VLENGTH);
  _box = _para.GetParameter<Vec3f>(PARA_BOX);
  Vec3f sigma = { static_cast<Real>(vtkm::Sqrt(_alpha / 2.0) * _box[0] / vtkm::Pi()),
                  static_cast<Real>(vtkm::Sqrt(_alpha / 2.0) * _box[1] / vtkm::Pi()),
                  static_cast<Real>(vtkm::Sqrt(_alpha / 2.0) * _box[2] / vtkm::Pi()) };
  _dt = _executioner->Dt();
  // prepare for RBE force
  auto velocity_type = _para.GetParameter<std::string>(gtest::velocity_type);
  auto random = (velocity_type != "TEST") ? true : false;
  RBEPSAMPLE rbe_presolve_psample = { _alpha, _Vlength, _box, _RBE_P };
  rbe_presolve_psample._RBE_random = random;
  _psample = rbe_presolve_psample.Fetch_P_Sample(Real(0.0), sigma);

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

vtkm::cont::ArrayHandle<Vec3f> ExecutionNVT::NearForceLJ()
{
  _nearforce_type = _para.GetParameter<std::string>(PARA_NEIGHBOR_TYPE);
  if (_nearforce_type == "RBL")
  {
    ComputeRBLLJForce(_all_force);
  }
  else if (_nearforce_type == "VERLETLIST")
  {
    ComputeVerletlistLJForce(_all_force);
  }
  else if (_nearforce_type == "ORIGINAL")
  {
    ComputeOriginalLJForce(_all_force);
  }
  return _all_force;
}

vtkm::cont::ArrayHandle<Vec3f> ExecutionNVT::NearForceEAM()
{
  _nearforce_type = _para.GetParameter<std::string>(PARA_NEIGHBOR_TYPE);
  if (_nearforce_type == "RBL")
  {
    ComputeRBLEAMForce(_all_force);
  }
  else if (_nearforce_type == "VERLETLIST")
  {
    ComputeVerletlistEAMForce(_all_force);
  }
  else if (_nearforce_type == "ORIGINAL")
  {
    ComputeOriginalEAMForce(_all_force);
  }
  return _all_force;
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
  Invoker{}(MolecularWorklet::ComputeBondHarmonicWorklet{ _box },
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
  Invoker{}(MolecularWorklet::ComputeAngleHarmonicWorklet{ _box },
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

  Invoker{}(MolecularWorklet::ComputeDihedralHarmonicWorklet{ _box },
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
  auto box = _para.GetParameter<Vec3f>(PARA_BOX);
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

    RunWorklet::ComputeSpecialCoulGeneral(box,
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
    RunWorklet::ComputeSpecialCoul(box, _atoms_id, groupVecArray, _force_function, _topology, _locator, _spec_coul_force);
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

vtkm::cont::ArrayHandle<Vec6f> ExecutionNVT::LJVirial()
{
  ComputeVerletlistLJVirial(_virial_atom_lj);
  return _virial_atom_lj;
}

vtkm::cont::ArrayHandle<Vec6f> ExecutionNVT::EwaldVirial()
{
  ComputeEwaldLongVirial(_Kmax, _virial_atom_ewald_long);
  return _virial_atom_ewald_long;
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
  auto shake = _para.GetParameter<std::string>(PARA_FIX_SHAKE);
  Real temperature_kB = _unit_factor._kB;
  if (_init_way == "inbuild") 
  {
    _tempT = 0.5 * _tempT_sum / (3 * n / 2.0);
  }
  else if (shake == "null")
  {
    _tempT = 0.5 * _tempT_sum / ((3 * n -3 ) * temperature_kB / 2.0);
  }
  else if (shake == "true")
  {
    _tempT = 0.5 * _tempT_sum / ((3 * n - n - 3) * temperature_kB / 2.0);
  }
  _para.SetParameter(PARA_TEMPT_SUM, _tempT_sum);
  _para.SetParameter(PARA_TEMPT, _tempT);
}

void ExecutionNVT::SetForceFunction()
{
  InitERF();
  if (_para.GetParameter<bool>(PARA_FAR_FORCE))
  {
    auto cut_off = _para.GetParameter<Real>(PARA_CUTOFF);
    auto volume = _para.GetParameter<Real>(PARA_VOLUME);
    auto vlength = _para.GetParameter<Real>(PARA_VLENGTH);
    auto box = _para.GetParameter<Vec3f>(PARA_BOX);
    _force_function.SetParameters(cut_off, _alpha, volume, vlength, box, _Kmax);
  }

  if (_para.GetParameter<std::string>(PARA_FILE_TYPE) == "EAM")
  {
    auto rhomax = _para.GetParameter<Real>(EAM_PARA_RHOMAX);
    auto nrho = _para.GetParameter<Id>(EAM_PARA_NRHO);
    auto drho = _para.GetParameter<Real>(EAM_PARA_DRHO);
    auto nr = _para.GetParameter<Id>(EAM_PARA_NR);
    auto dr = _para.GetParameter<Real>(EAM_PARA_DR);
    _force_function.SetEAMParameters(rhomax, nrho, drho, nr, dr);
  }

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

  if (_init_way != "inbuild")
  {
    auto source_array = _para.GetFieldAsArrayHandle<Id>(field::special_source_array);
    auto offsets_array = _para.GetFieldAsArrayHandle<Id>(field::special_offsets_array);
    _topology.SetSourceAndOffsets(source_array, offsets_array);
  }
}

void ExecutionNVT::InitParameters()
{
  ExecutionMD::InitParameters();
  if (_para.GetParameter<bool>(PARA_FAR_FORCE))
  {
    _RBE_P = _para.GetParameter<IdComponent>(PARA_COULOMB_SAMPLE_NUM);
    _alpha = _para.GetParameter<Real>(PARA_ALPHA);
    _Kmax = _para.GetParameter<IdComponent>(PARA_KMAX); 
  }
  _kbT = _para.GetParameter<std::vector<Real>>(PARA_TEMPERATURE)[0]; 
  _Tdamp = _para.GetParameter<std::vector<Real>>(PARA_TEMPERATURE)[2];

  _Pstart = _para.GetParameter<std::vector<Real>>(PARA_PRESSURE_VECTOR)[0];
  _Pstop = _para.GetParameter<std::vector<Real>>(PARA_PRESSURE_VECTOR)[1];
  _Pperiod = _para.GetParameter<std::vector<Real>>(PARA_PRESSURE_VECTOR)[2];
  _bulkmodulus = _para.GetParameter<std::vector<Real>>(PARA_PRESSURE_VECTOR)[3];

  _para.SetParameter(PARA_TEMPT_SUM, Real{ 0.0 });
  _para.SetParameter(PARA_TEMPT, Real{ 0.0 });
  _para.SetParameter(PARA_PRESSURE, Real{ 0.0 });
  _init_way = _para.GetParameter<std::string>(PARA_INIT_WAY);
}

void ExecutionNVT::TimeIntegration() {}

void ExecutionNVT::ConstraintA()
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
    MolecularWorklet::NewConstraintAWaterBondAngleWorklet{ _box, _dt, _unit_factor._fmt2v, range },
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
  auto angle_list = _para.GetFieldAsArrayHandle<Id>(field::angle_atom_id);
  auto&& anglelist_group = vtkm::cont::make_ArrayHandleGroupVec<3>(angle_list);

  auto data_range = _para.GetParameter<vtkm::Vec<vtkm::Range, 3>>(PARA_RANGE);
  Vec<Vec2f, 3> range{
    { static_cast<Real>(data_range[0].Min), static_cast<Real>(data_range[0].Max) },
    { static_cast<Real>(data_range[1].Min), static_cast<Real>(data_range[1].Max) },
    { static_cast<Real>(data_range[2].Min), static_cast<Real>(data_range[2].Max) }
  };

  Invoker{}(
    MolecularWorklet::NewConstraintBWaterBondAngleWorklet{ _box, _dt, _unit_factor._fmt2v, range },
    anglelist_group,
    _position,
    _old_velocity,
    _all_force,
    _mass,
    _locator);



  vtkm::cont::ArrayCopy(_old_velocity, _velocity);
}

void ExecutionNVT::ReadPotentialFile(std::ifstream& input_file)
{
  if (!input_file.is_open())
  {
    std::cerr << "Unable to open the file." << std::endl;
  }
  for (int i = 0; i < 2; ++i)
  {
    input_file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  }
  input_file >> file.nrho >> file.drho >> file.nr >> file.dr >> file.cut_off;

  file.frho.resize(file.nrho + 1);
  file.zr.resize(file.nr + 1);
  file.rhor.resize(file.nrho + 1);

  for (int i = 0; i < file.nrho; ++i)
  {
    input_file >> file.frho[i];
  }

  for (int i = 0; i < file.nr; ++i)
  {
    input_file >> file.zr[i];
  }

  for (int i = 0; i < file.nrho; ++i)
  {
    input_file >> file.rhor[i];
  }
  input_file.close();
}

void ExecutionNVT::AllocateEAM() {}

void ExecutionNVT::file2array()
{
  Id i, j, k, m, n;
  Real sixth = 1.0 / 6.0;

  Real rmax;
  dr = drho = rmax = rhomax = 0.0;

  dr = MAX(dr, file.dr);
  drho = MAX(drho, file.drho);
  rmax = MAX(rmax, (file.nr - 1) * file.dr);
  rhomax = MAX(rhomax, (file.nrho - 1) * file.drho);

  // set nr,nrho from cutoff and spacings
  // 0.5 is for round-off in divide

  nr = static_cast<int>(rmax / dr + 0.5);
  nrho = static_cast<int>(rhomax / drho + 0.5);

  // ------------------------------------------------------------------
  // setup frho arrays
  // ------------------------------------------------------------------
  frho.resize(nrho + 1);

  Real r, p, cof1, cof2, cof3, cof4;
  for (m = 1; m <= nrho; m++)
  {
    r = (m - 1) * drho;
    p = r / file.drho + 1.0;
    k = static_cast<int>(p);
    k = MIN(k, file.nrho - 2);
    k = MAX(k, 2);
    p -= k;
    p = MIN(p, 2.0);
    cof1 = -sixth * p * (p - 1.0) * (p - 2.0);
    cof2 = 0.5 * (p * p - 1.0) * (p - 2.0);
    cof3 = -0.5 * p * (p + 1.0) * (p - 2.0);
    cof4 = sixth * p * (p * p - 1.0);
    frho[m] = cof1 * file.frho[k - 1] + cof2 * file.frho[k] + cof3 * file.frho[k + 1] +
      cof4 * file.frho[k + 2];
  }

  // ------------------------------------------------------------------
  // setup rhor arrays
  // ------------------------------------------------------------------
  rhor.resize(nrho + 1);
  for (m = 1; m <= nr; m++)
  {
    r = (m - 1) * dr;
    p = r / file.dr + 1.0;
    k = static_cast<int>(p);
    k = MIN(k, file.nr - 2);
    k = MAX(k, 2);
    p -= k;
    p = MIN(p, 2.0);
    auto cof1 = -sixth * p * (p - 1.0) * (p - 2.0);
    auto cof2 = 0.5 * (p * p - 1.0) * (p - 2.0);
    auto cof3 = -0.5 * p * (p + 1.0) * (p - 2.0);
    auto cof4 = sixth * p * (p * p - 1.0);
    rhor[m] = cof1 * file.rhor[k - 1] + cof2 * file.rhor[k] + cof3 * file.rhor[k + 1] +
      cof4 * file.rhor[k + 2];
  }

  // ------------------------------------------------------------------
  // setup z2r arrays
  // ------------------------------------------------------------------
  z2r.resize(nr + 1);

  double zri;
  for (m = 1; m <= nr; m++)
  {
    r = (m - 1) * dr;

    p = r / file.dr + 1.0;
    k = static_cast<int>(p);
    k = MIN(k, file.nr - 2);
    k = MAX(k, 2);
    p -= k;
    p = MIN(p, 2.0);
    cof1 = -sixth * p * (p - 1.0) * (p - 2.0);
    cof2 = 0.5 * (p * p - 1.0) * (p - 2.0);
    cof3 = -0.5 * p * (p + 1.0) * (p - 2.0);
    cof4 = sixth * p * (p * p - 1.0);
    zri = cof1 * file.zr[k - 1] + cof2 * file.zr[k] + cof3 * file.zr[k + 1] + cof4 * file.zr[k + 2];

    z2r[m] = 27.2 * 0.529 * zri * zri;
  }
}

void ExecutionNVT::interpolate(Id n, Real delta, std::vector<Real>& f, std::vector<Vec7f>& spline)
{
  for (int m = 1; m <= n; m++)
  {
    spline[m][6] = f[m];
  }

  spline[1][5] = spline[2][6] -
    spline[1][6]; //f'(x) = (f(x + h) - f(x)) / h    [5] is the coefficient of the first derivative(energy expression)
  spline[2][5] = 0.5 * (spline[3][6] - spline[1][6]);
  spline[n - 1][5] = 0.5 * (spline[n][6] - spline[n - 2][6]);
  spline[n][5] = spline[n][6] - spline[n - 1][6];

  for (int m = 3; m <= n - 2; m++)
  {
    spline[m][5] =
      ((spline[m - 2][6] - spline[m + 2][6]) + 8.0 * (spline[m + 1][6] - spline[m - 1][6])) /
      12.0; //further sample points for a more accurate estimate
  }

  for (int m = 1; m <= n - 1; m++)
  {
    spline[m][4] = 3.0 * (spline[m + 1][6] - spline[m][6]) - 2.0 * spline[m][5] -
      spline[m + 1][5]; //[4] is the coefficient of the second derivative
    spline[m][3] = spline[m][5] + spline[m + 1][5] -
      2.0 * (spline[m + 1][6] - spline[m][6]); // [3] is the coefficient of the third derivative
  }

  spline[n][4] = 0.0;
  spline[n][3] = 0.0; //The second and third derivative coefficients at the last sample point are zero,
                      // To make the interpolation curve smoother at both ends, the higher derivative coefficient at the boundary can be set to zero.
                      // This is because spline interpolation typically uses higher-order polynomial interpolation
                      // at the inner sample points and lower-order polynomials at the boundaries to ensure smoothness.

  for (int m = 1; m <= n; m++)
  {
    spline[m][2] = spline[m][5] / delta;       //The coefficient of the second derivative(force expression)  
    spline[m][1] = 2.0 * spline[m][4] / delta; //The coefficient of the first derivative
    spline[m][0] = 3.0 * spline[m][3] / delta; //The coefficient of the zero derivative (i.e. the value of the function).
  }
}

void ExecutionNVT::array2spline()
{
  frho_spline.resize(nrho + 1);
  rhor_spline.resize(nrho + 1);
  z2r_spline.resize(nr + 1);

  interpolate(nrho, drho, frho, frho_spline);
  interpolate(nr, dr, rhor, rhor_spline);
  interpolate(nr, dr, z2r, z2r_spline);
}

void ExecutionNVT::SetEAM()
{
  _para.SetParameter(EAM_PARA_CUTOFF, file.cut_off);
  _para.SetParameter(EAM_PARA_RHOMAX, rhomax);
  _para.SetParameter(EAM_PARA_NRHO, nrho);
  _para.SetParameter(EAM_PARA_DRHO, drho);
  _para.SetParameter(EAM_PARA_NR, nr);
  _para.SetParameter(EAM_PARA_DR, dr);

  _para.AddField(field::rhor_spline, ArrayHandle<Vec7f>{});
  auto rhor_spline_get = _para.GetFieldAsArrayHandle<Vec7f>(field::rhor_spline);
  vtkm::cont::ArrayCopy(vtkm::cont::make_ArrayHandle<Vec7f>(rhor_spline), rhor_spline_get);

  _para.AddField(field::frho_spline, ArrayHandle<Vec7f>{});
  auto frho_spline_get = _para.GetFieldAsArrayHandle<Vec7f>(field::frho_spline);
  vtkm::cont::ArrayCopy(vtkm::cont::make_ArrayHandle<Vec7f>(frho_spline), frho_spline_get);

  _para.AddField(field::z2r_spline, ArrayHandle<Vec7f>{});
  auto z2r_spline_get = _para.GetFieldAsArrayHandle<Vec7f>(field::z2r_spline);
  vtkm::cont::ArrayCopy(vtkm::cont::make_ArrayHandle<Vec7f>(z2r_spline), z2r_spline_get);
}

void ExecutionNVT::InitStyle()
{
  AllocateEAM();
  file2array();
  array2spline();
  SetEAM();
}

void ExecutionNVT::fix_press_berendsen()
{
  // compute new T,P

  //ComputeTempe();
  Compute_Pressure_Scalar();
  Couple();


  //Compute dilation
  auto currentstep = _app.GetExecutioner()->CurrentStep();
  auto beginstep = 0;
  auto endstep = _app.GetExecutioner()->NumStep();

  Real delta = currentstep - beginstep;
  if (delta != 0.0)
  {
    delta = delta / static_cast<Real>(endstep - beginstep);
  }

  for (int i = 0; i < 3; i++)
  {
    auto dt_over_period = _dt / p_period[i];
    auto bulkmodulus_inv = 1.0 / _bulkmodulus;
    p_target[i] = p_start[i] + delta * (p_stop[i] - p_start[i]);

    dilation[i] =
      vtkm::Pow(1.0 - dt_over_period * (p_target[0] - p_current[i]) * bulkmodulus_inv, 1.0 / 3.0);
  }

  //remap simulation box and atoms
  remap();
  ApplyPbc();
  _locator.SetPosition(_position);
}

void ExecutionNVT::fix_press_berendsen_scale()
{
  // compute new T,P

  //ComputeTempe();
  Compute_Pressure_Scalar();
  Couple();

  auto currentstep = _app.GetExecutioner()->CurrentStep();
  auto beginstep = 0;
  auto endstep = _app.GetExecutioner()->NumStep();

  Real delta = currentstep - beginstep;
  if (delta != 0.0)
  {
    delta = delta / static_cast<Real>(endstep - beginstep);
  }


  for (int i = 0; i < 3; i++)
  {
    p_target[i] = p_start[i] + delta * (p_stop[i] - p_start[i]);
  }
  _pressure_coupling = _dt / (p_period[0] * 3 * _bulkmodulus);
  _scale_factor = 1.0 -
    _pressure_coupling * (p_target[0] - (p_current[0] + p_current[1] + p_current[2]) * 0.3333);


  // change range
  auto range = _para.GetParameter<vtkm::Vec<vtkm::Range, 3>>(PARA_RANGE);
  for (int i = 0; i < 3; ++i)
  {
    Real center = (range[i].Max + range[i].Min) / 2.0;
    Real half_length = (range[i].Max - range[i].Min) / 2.0;
    Real new_half_length = half_length * _scale_factor;

    range[i].Min = center - new_half_length;
    range[i].Max = center + new_half_length;
  }
  std::cout << "scale_factor =" << _scale_factor << ",range.Min=" << range[0].Min
            << ",range.Max=" << range[0].Max << std::endl;
  _para.SetParameter(PARA_RANGE, range);

  Vec3f left_bottom{ Real(range[0].Min), Real(range[1].Min), Real(range[2].Min) };
  Vec3f right_top{ Real(range[0].Max), Real(range[1].Max), Real(range[2].Max) };
  _locator.SetRange(left_bottom, right_top);

  // change box
  Vec3f box{ Real(range[0].Max - range[0].Min),
             Real(range[1].Max - range[1].Min),
             Real(range[2].Max - range[2].Min)
  };
  _para.SetParameter(PARA_BOX, box);

  // change position
  RunWorklet::fix_press_berendsen(_scale_factor, _position, _locator);
  ApplyPbc();
  _locator.SetPosition(_position);
}

void ExecutionNVT::Compute_Pressure_Scalar()
{
  //compute temperature
  auto temperature = _para.GetParameter<Real>(PARA_TEMPT);

  // compute  virial
  auto range = _para.GetParameter<vtkm::Vec<vtkm::Range, 3>>(PARA_RANGE);
  auto volume =
    (range[0].Max - range[0].Min) * (range[1].Max - range[1].Min) * (range[2].Max - range[2].Min);
  auto inv_volume = 1.0 / volume;

  ComputeVirial();

  //compute dof
  auto n = _position.GetNumberOfValues();
  auto extra_dof = 3; //dimension =3
  auto dof = 3 * n - extra_dof;

  //compute pressure_scalar
  _pressure_scalar = (dof * _unit_factor._boltz * temperature + virial[0] + virial[1] + virial[2]) /
    3.0 * inv_volume * _unit_factor._nktv2p;

  _para.SetParameter(PARA_PRESSURE, _pressure_scalar);
  std::cout << " pressure=" << _pressure_scalar << std::endl;
}

void ExecutionNVT::ComputeVirial()
{
  ComputeVerletlistLJVirial(_virial_atom_lj); //_virial_atom
  auto virial_lj = vtkm::cont::Algorithm::Reduce(_virial_atom_lj, vtkm::TypeTraits<Vec6f>::ZeroInitialization());

  ComputeEwaldCoulVirial(_virial_atom_ewald_coul);
   auto virial_ewald_coul = vtkm::cont::Algorithm::Reduce(_virial_atom_ewald_coul, vtkm::TypeTraits<Vec6f>::ZeroInitialization()) *
                      _unit_factor._qqr2e;

  ComputeEwaldLongVirial(_Kmax, _virial_atom_ewald_long);
  auto virial_ewald_long = vtkm::cont::Algorithm::Reduce(_virial_atom_ewald_long, vtkm::TypeTraits<Vec6f>::ZeroInitialization())
                    * _unit_factor._qqr2e;

  virial = virial_lj + virial_ewald_coul + virial_ewald_long;

  //RunWorklet::SumVirial(LJVirial(), EwaldVirial(), _virial_atom);
  //reduce virial_atom
  //virial = { 0, 0, 0, 0, 0, 0 };
  //for (int i = 0; i < _virial_atom.GetNumberOfValues(); ++i)
  //{
  //  Vec6f vatom = _virial_atom.ReadPortal().Get(i);
  //  for (int j = 0; j < 6; ++j)
  //  {
  //    virial[j] += vatom[j];
  //  }
  //}

 // virial =
   // vtkm::cont::Algorithm::Reduce(_virial_atom, vtkm::TypeTraits<Vec6f>::ZeroInitialization());
}

void ExecutionNVT::Couple()
{
  p_current[0] = p_current[1] = p_current[2] = _pressure_scalar;
}

void ExecutionNVT::set_global_box()
{
  vtkm::Vec<vtkm::Range, 3> range = _para.GetParameter<vtkm::Vec<vtkm::Range, 3>>(PARA_RANGE);
  prd[0] = xprd = range[0].Max - range[0].Min;
  prd[1] = yprd = range[1].Max - range[1].Min;
  prd[2] = zprd = range[2].Max - range[2].Min;

  h[0] = xprd;
  h[1] = yprd;
  h[2] = zprd;
  h[3] = 0;
  h[4] = 0;
  h[5] = 0;

  Vec3f box{ h[0], h[1], h[2] };
  _para.SetParameter(PARA_BOX, box);


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

void ExecutionNVT::ApplyPbc()
{
  //pbc
  auto&& position_flag = _para.GetFieldAsArrayHandle<Id3>(field::position_flag);
  auto box = _para.GetParameter<Vec3f>(PARA_BOX); //
  auto range = _para.GetParameter<vtkm::Vec<vtkm::Range, 3>>(PARA_RANGE);
  Vec<Vec2f, 3> data_range{ { static_cast<Real>(range[0].Min), static_cast<Real>(range[0].Max) },
                            { static_cast<Real>(range[1].Min), static_cast<Real>(range[1].Max) },
                            { static_cast<Real>(range[2].Min), static_cast<Real>(range[2].Max) }};
  RunWorklet::ApplyPbcFlag(box, data_range, _position, _locator, position_flag);

}

void ExecutionNVT::x2lamda(Id n)
{
  //
  vtkm::Vec<vtkm::Range, 3> range = _para.GetParameter<vtkm::Vec<vtkm::Range, 3>>(PARA_RANGE);
  RunWorklet::X2Lamda(h_inv, range, _position);
  _locator.SetPosition(_position);


  //Vec3f delta;
  //vtkm::Vec<vtkm::Range, 3> range = _para.GetParameter<vtkm::Vec<vtkm::Range, 3>>(PARA_RANGE);

  //for (int i = 0; i < n; i++)
  //{
  //  delta[0] = _position.ReadPortal().Get(i)[0] - range[0].Min;
  //  delta[1] = _position.ReadPortal().Get(i)[1] - range[1].Min;
  //  delta[2] = _position.ReadPortal().Get(i)[2] - range[2].Min;

  //  Vec3f lamda_position = { h_inv[0] * delta[0] + h_inv[5] * delta[1] + h_inv[4] * delta[2],
  //                           h_inv[1] * delta[1] + h_inv[3] * delta[2],
  //                           h_inv[2] * delta[2] };

  //  _position.WritePortal().Set(i, lamda_position);
  //}
  //_locator.SetPosition(_position);
}

void ExecutionNVT::lamda2x(Id n)
{
  vtkm::Vec<vtkm::Range, 3> range = _para.GetParameter<vtkm::Vec<vtkm::Range, 3>>(PARA_RANGE);
  RunWorklet::Lamda2X(h, range, _position);
  _locator.SetPosition(_position);

  //vtkm::Vec<vtkm::Range, 3> range = _para.GetParameter<vtkm::Vec<vtkm::Range, 3>>(PARA_RANGE);
  //Vec3f position_base;

  //for (int i = 0; i < n; i++)
  //{
  //  position_base[0] = _position.ReadPortal().Get(i)[0];
  //  position_base[1] = _position.ReadPortal().Get(i)[1];
  //  position_base[2] = _position.ReadPortal().Get(i)[2];

  //  Vec3f x_position = { h[0] * position_base[0] + h[5] * position_base[1] +
  //                         h[4] * position_base[2] + static_cast<vtkm::FloatDefault>(range[0].Min),
  //                       h[1] * position_base[1] + h[3] * position_base[2] +
  //                         static_cast<vtkm::FloatDefault>(range[1].Min),
  //                       h[2] * position_base[2] + static_cast<vtkm::FloatDefault>(range[2].Min) };

  //  _position.WritePortal().Set(i, x_position);
  //}
  //_locator.SetPosition(_position);
}

void ExecutionNVT::remap()
{
  Real oldlo, oldhi, ctr;
  auto n = _position.GetNumberOfValues();

  //  convert lamda coords
  x2lamda(n);

  // change range and box 
  vtkm::Vec<vtkm::Range, 3> range = _para.GetParameter<vtkm::Vec<vtkm::Range, 3>>(PARA_RANGE);
  for (int i = 0; i < 3; i++)
  {
    oldlo = range[i].Min;
    oldhi = range[i].Max;
    ctr = 0.5 * (oldlo + oldhi);
    range[i].Min = (oldlo - ctr) * dilation[i] + ctr;
    range[i].Max = (oldhi - ctr) * dilation[i] + ctr;
  }
  std::cout << "range.Min=" << range[0].Min << ",range.Max=" << range[0].Max << std::endl;
  //
  _para.SetParameter(PARA_RANGE, range);


  Vec3f left_bottom{ Real(range[0].Min), Real(range[1].Min), Real(range[2].Min) };
  Vec3f right_top{ Real(range[0].Max), Real(range[1].Max), Real(range[2].Max) };
  _locator.SetRange(left_bottom, right_top);

  set_global_box();

  // convert real coords
  lamda2x(n);
}