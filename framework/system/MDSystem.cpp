#include "FieldName.h"
#include "MDSystem.h"
#include "system/worklet/MolecularWorklet.h"
#include "system/worklet/SystemWorklet.h"
#include <vtkm/cont/Algorithm.h>
#include <vtkm/cont/ArrayHandleGroupVec.h>
#include "Executioner.h"
#include "ERFTable.h"
#include <vtkm/worklet/Keys.h>
#include <vtkm/cont/ArrayHandleConstant.h>
//#include <vtkm/exec/ParallelExecutor.h>
//#include <vtkm/cont/BlockForEach.h>
#include <vtkm/cont/DeviceAdapter.h>
#include <vtkm/cont/ArrayHandleSOA.h>
#include <vtkm/cont/EnvironmentTracker.h>

struct SetIndex : vtkm::worklet::WorkletMapField
{
  SetIndex(const Id& num)
    : _pos_num(num)
  {
  }
  using ControlSignature = void(FieldIn ids, FieldOut p);
  using ExecutionSignature = void(_1, _2);

  //template<typename CoordType>
  VTKM_EXEC void operator()(const Id& ids, Id& p) const { p = ids / _pos_num; }

  Id _pos_num;
};

MDSystem::MDSystem(const Configuration& cfg)
  : MeshFreeSystem(cfg)
  , _use_erf(Get<bool>("use_erf"))
  , _RBE_P(Get<IdComponent>("rbeP"))
{

}

void MDSystem::Init()
{
  MeshFreeSystem::Init();
  SetForceFunction();
  SetTopology();

  vtkm::cont::Timer timer;
  timer.Start();
  auto N = _position.GetNumberOfValues();
  _rhok_Re.AllocateAndFill(_RBE_P * N, 0.0);
  _rhok_Im.AllocateAndFill(_RBE_P * N, 0.0);

  vtkm::cont::ArrayHandleIndex indexArray(N * _RBE_P);
  vtkm::cont::Invoker{}(SetIndex{ N }, indexArray, _psamplekey);
  std::cout << "init rhok: " << timer.GetElapsedTime() << std::endl;

  _EleNearPairtimer_counting = 0.0;
}

void MDSystem::InitialCondition()
{
  // 物理场初始化
  _charge = GetFieldAsArrayHandle<Real>(field::charge);
  _topology.SetCharge(_charge);
  _velocity = GetFieldAsArrayHandle<Vec3f>(field::velocity);
  _mass = GetFieldAsArrayHandle<Real>(field::mass);

  _molecule_id = GetFieldAsArrayHandle<Id>(field::molecule_id);
  _atoms_id = GetFieldAsArrayHandle<Id>(field::atom_id);
}

void MDSystem::SetForceFunction() {}

void MDSystem::SetTopology() {}

void MDSystem::InitERF() 
{
  if (_use_erf == true)
  {     
     _static_table.SetTableIndex1(vtkm::cont::make_ArrayHandle(erf_index_leq2_v));
     _static_table.SetTableIndex2(vtkm::cont::make_ArrayHandle(erf_index_geq2_v));
     _static_table.SetTableRij1(vtkm::cont::make_ArrayHandle(erf_dis_leq2_v));
     _static_table.SetTableRij2(vtkm::cont::make_ArrayHandle(erf_dis_geq2_v));
     _static_table.SetTabledRij1(vtkm::cont::make_ArrayHandle(erf_dis_dif_leq2_v));
     _static_table.SetTabledRij2(vtkm::cont::make_ArrayHandle(erf_dis_dif_geq2_v));
     _static_table.SetTableFunctionRij1(vtkm::cont::make_ArrayHandle(erf_gnear_leq2_v));
     _static_table.SetTableFunctionRij2(vtkm::cont::make_ArrayHandle(erf_gnear_geq2_v));
     _static_table.SetTabledFunctionRij1(vtkm::cont::make_ArrayHandle(erf_gnear_der_leq2_v));
     _static_table.SetTabledFunctionRij2(vtkm::cont::make_ArrayHandle(erf_gnear_der_geq2_v));
  }
}

void MDSystem::InitField()
{
  AddField(field::charge, ArrayHandle<Real>{});
  AddField(field::velocity, ArrayHandle<Vec3f>{});
  AddField(field::mass, ArrayHandle<Real>{});
  AddField(field::molecule_id, ArrayHandle<Id>{});
  AddField(field::atom_id, ArrayHandle<Id>{});
  AddField(field::position, ArrayHandle<Vec3f>{});
  AddField(field::special_offsets, ArrayHandle<Id>{});
  AddField(field::special_weights, ArrayHandle<Real>{});
  AddField(field::special_ids, ArrayHandle<Id>{});
}

void MDSystem::UpdateVerletList() 
{
  auto cut_off = GetParameter<Real>(PARA_CUTOFF);
  auto N = _position.GetNumberOfValues();
  auto verletlist_num = N * N;

  ArrayHandle<Id> id_verletlist;
  id_verletlist.Allocate(verletlist_num);

  std::vector<Id> offset_vec(N + 1);
  Id inc = 0;
  std::generate(
    offset_vec.begin(), offset_vec.end(), [&](void) -> Id { return (inc++) * N; });
  ArrayHandle<Id> temp_offset = vtkm::cont::make_ArrayHandle(offset_vec);

  auto id_verletlist_group =
    vtkm::cont::make_ArrayHandleGroupVecVariable(id_verletlist, temp_offset);

  ArrayHandle<Id> num_verletlist;
  ArrayHandle<Vec3f> offset_verletlist;
  offset_verletlist.Allocate(N * N);

  auto offset_verletlist_group =
    vtkm::cont::make_ArrayHandleGroupVecVariable(offset_verletlist, temp_offset);

  SystemWorklet::ComputeNeighbours(
    cut_off, _atoms_id, _locator, id_verletlist_group, num_verletlist, offset_verletlist_group);

  _locator.SetVerletListInfo(num_verletlist, id_verletlist_group, offset_verletlist_group);
}

void MDSystem::ComputeCorrForce(vtkm::cont::ArrayHandle<Vec3f>& corr_force)
{
  auto rc = GetParameter<Real>(PARA_CUTOFF);
  auto rs = GetParameter<Real>(PARA_RS);
  auto rho_system = GetParameter<Real>(PARA_RHO);
  auto N = _position.GetNumberOfValues();

  vtkm::cont::ArrayHandle<vtkm::Id> num_verletlist;
  vtkm::cont::ArrayHandle<vtkm::Id> id_verletlist;
  vtkm::cont::ArrayHandle<vtkm::Vec3f> offset_verletlist;

  Id rs_num = rho_system * vtkm::Ceil(4.0 / 3.0 * vtkm::Pif() * rs * rs * rs) + 1;
  Id rc_num = rho_system * vtkm::Ceil(4.0 / 3.0 * vtkm::Pif() * rc * rc * rc) + 1;
  Id rcs_num = rc_num - rs_num;
  Real random_num = GetParameter<Real>(PARA_RANDOM_NUM);
  Real random_rate = random_num / (rcs_num);
  Id pice_num = std::ceil(1.0 / random_rate);
  Id cutoff_num = rs_num + random_num;
  Id verletlist_num = N * cutoff_num;

  num_verletlist.Allocate(N * 2);
  id_verletlist.Allocate(verletlist_num);
  offset_verletlist.Allocate(verletlist_num);

  std::vector<Id> stride_vec(N + 1);
  Id inc = 0;
  std::generate(stride_vec.begin(), stride_vec.end(), [&](void) -> Id { return (inc++) * cutoff_num; });
  vtkm::cont::ArrayHandle<vtkm::Id> stride_array = vtkm::cont::make_ArrayHandle(stride_vec);

  auto id_verletlist_group = vtkm::cont::make_ArrayHandleGroupVecVariable(id_verletlist, stride_array);
  auto offset_verletlist_group = vtkm::cont::make_ArrayHandleGroupVecVariable(offset_verletlist, stride_array);
  auto num_verletlist_group = vtkm::cont::make_ArrayHandleGroupVec<2>(num_verletlist);

  SystemWorklet::ComputeRBLNeighboursOnce(rs_num,
                                          pice_num,
                                          _atoms_id,
                                          _locator,
                                          id_verletlist_group,
                                          num_verletlist_group,
                                          offset_verletlist_group);

  SystemWorklet::NearForceRBLERF(rs_num,
                                 pice_num,
                                 _unit_factor._qqr2e,
                                 _atoms_id,
                                 _locator,
                                 _topology,
                                 _force_function,
                                 _static_table,
                                 id_verletlist_group,
                                 num_verletlist_group,
                                 offset_verletlist_group,
                                 corr_force);
}

std::vector<Vec2f> MDSystem::ComputeChargeStructureFactorRBE(Real& _Vlength, ArrayHandle<Vec3f>& _psample)
{
  std::vector<Vec2f> rhok;
  auto p_number = _psample.GetNumberOfValues();
  ArrayHandle<Real> density_real;
  ArrayHandle<Real> density_image;

  //auto sum = 0.0f;
  //vtkm::cont::Timer timer;
  //#pragma omp parallel for
  for (Id i = 0; i < p_number; i++)
  {
    auto kl = _psample.ReadPortal().Get(i);
    kl = 2 * vtkm::Pi() * kl / _Vlength;
    SystemWorklet::ComputeChargeStructureFactorComponent(
      kl, _position, _charge, density_real, density_image);
    //timer.Start();
    Real value_Re =
      vtkm::cont::Algorithm::Reduce(density_real, vtkm::TypeTraits<Real>::ZeroInitialization());
    Real value_Im =
      vtkm::cont::Algorithm::Reduce(density_image, vtkm::TypeTraits<Real>::ZeroInitialization());
    //sum += timer.GetElapsedTime();
    rhok.push_back({ value_Re, value_Im });
  }
  //std::cout << "timer: " << sum << std::endl;
  return rhok;
}

std::vector<Vec2f> MDSystem::ComputeChargeStructureFactorEwald(Real& _Vlength, IdComponent& Kmax)
{
  std::vector<Vec2f> rhok;
  ArrayHandle<Real> density_real;
  ArrayHandle<Real> density_image;
  for (Id i = -Kmax; i <= Kmax; i++)
  {
    for (Id j = -Kmax; j <= Kmax; j++)
    {
      for (Id k = -Kmax; k <= Kmax; k++)
      {
        if (!(i == 0 && j == 0 && k == 0))
        {
          Vec3f K = { Real(i), Real(j), Real(k) };
          K = 2 * vtkm::Pi() * K / _Vlength;
          SystemWorklet::ComputeChargeStructureFactorComponent(
            K, _position, _charge, density_real, density_image);
          Real value_Re = vtkm::cont::Algorithm::Reduce(
            density_real, vtkm::TypeTraits<Real>::ZeroInitialization());
          Real value_Im = vtkm::cont::Algorithm::Reduce(
            density_image, vtkm::TypeTraits<Real>::ZeroInitialization());

          rhok.push_back({ value_Re, value_Im });
        }
      }
    }
  }
  return rhok;
}

void MDSystem::ComputeRBEEleForce(ArrayHandle<Vec3f>& psample,
                                  IdComponent& RBE_P,
                                  ArrayHandle<Vec3f>& RBE_ele_force)
{
 

  auto Vlength = GetParameter<Real>(PARA_VLENGTH);
 

  ArrayHandle<Vec2f> new_whole_rhok;
  ComputeNewChargeStructureFactorRBE(Vlength, psample, new_whole_rhok);
  
  SystemWorklet::ComputeNewRBEForce(RBE_P,
                                    _atoms_id,
                                    psample,
                                    /*whole_rhok,*/ new_whole_rhok, _force_function,
                                    _topology,
                                    _locator,
                                    RBE_ele_force);
 

}

struct ReduceWorklet : vtkm::worklet::WorkletReduceByKey
{
  using ControlSignature = void(KeysIn p, ValuesIn value, ReducedValuesOut reduce_value);

  using ExecutionSignature = _3(_2);
  using InputDomain = _1;
  template<typename ForceVecType>
  VTKM_EXEC typename ForceVecType::ComponentType operator()(const ForceVecType& value) const
  {
    typename ForceVecType::ComponentType sum = 0;
    for (vtkm::IdComponent index = 0; index < value.GetNumberOfComponents(); index++)
    {
      sum = sum + value[index];
    }
    return sum;
  }
};


struct SetValue : vtkm::worklet::WorkletMapField
{
  using ControlSignature = void(FieldInOut array_);
  using ExecutionSignature = void(_1);

  template<typename CoordType>
  VTKM_EXEC void operator()(CoordType& array_) const
  {
    array_ = 0;
  }
};


struct SetValue_New : vtkm::worklet::WorkletMapField
{
  SetValue_New(const Id& pos, const Id& num)
    : _pos_num(pos)
    , _pnumber(num)
  {

  }

  using ControlSignature = void(FieldIn ids, WholeArrayInOut array_);
  using ExecutionSignature = void(_1,_2);

  template<typename CoordType>
  VTKM_EXEC void operator()(const Id& id, CoordType& array_) const
  {
    for (int i = 0; i < _pnumber; ++i)
    {
      auto index = id + i * _pos_num;
      array_.Set(index, 0);
    }
  }

  Id _pos_num;
  Id _pnumber;
};

void  MDSystem::ComputeNewChargeStructureFactorRBE(Real& _Vlength,
                                                                ArrayHandle<Vec3f>& _psample,
                                                                ArrayHandle<Vec2f>& new_rhok)
{
 
  auto N = _position.GetNumberOfValues();
  
  const Id p_number = _psample.GetNumberOfValues();
  


  SystemWorklet::ComputePnumberChargeStructureFactor(_Vlength,
                                                     p_number,
                                                     N,
                                                     _atoms_id,
                                                     _position,
                                                     _charge,
                                                     _psample,
                                                     //psamplekey,
                                                     /*rhok_Re,*/_rhok_Re,
                                                     /*rhok_Im*/_rhok_Im);
 

  vtkm::cont::ArrayHandle<Id> psamplekey_out;
  vtkm::cont::ArrayHandle<Real> rhok_Re_reduce;
  vtkm::cont::ArrayHandle<Real> rhok_Im_reduce;

   
  vtkm::cont::Algorithm::ReduceByKey<Id, Real>(
    /*psamplekey,*/ _psamplekey, /*rhok_Re,*/ _rhok_Re, psamplekey_out, rhok_Re_reduce, vtkm::Add());
  vtkm::cont::Algorithm::ReduceByKey<Id, Real>(
    /*psamplekey,*/ _psamplekey, /*rhok_Im,*/ _rhok_Im, psamplekey_out, rhok_Im_reduce, vtkm::Add());
 


  SystemWorklet::ChangePnumberChargeStructureFactor(rhok_Re_reduce, rhok_Im_reduce, new_rhok);
 

  
}

void MDSystem::ComputeEwaldEleForce(IdComponent& Kmax, ArrayHandle<Vec3f>& Ewald_ele_force)
{
  auto Vlength = GetParameter<Real>(PARA_VLENGTH);
  auto rhok = ComputeChargeStructureFactorEwald(Vlength, Kmax);
  ArrayHandle<Vec2f> whole_rhok = vtkm::cont::make_ArrayHandle(rhok);
  SystemWorklet::ComputeNewFarElectrostatics(
    Kmax, _atoms_id, whole_rhok, _force_function, _topology, _locator, Ewald_ele_force);
}

void MDSystem::ComputeRBLNearForce(ArrayHandle<Vec3f>& nearforce)
{
  auto N = _position.GetNumberOfValues();
  vtkm::cont::ArrayHandle<Vec3f> corr_force;
  corr_force.Allocate(N);

  auto rc = GetParameter<Real>(PARA_CUTOFF);
  auto rs = GetParameter<Real>(PARA_RS);
  auto rho_system = GetParameter<Real>(PARA_RHO);

  vtkm::cont::ArrayHandle<vtkm::Id> num_verletlist;
  vtkm::cont::ArrayHandle<vtkm::Id> id_verletlist;
  vtkm::cont::ArrayHandle<vtkm::Vec3f> offset_verletlist;

  Id rs_num =  6* rho_system * vtkm::Ceil(4.0 / 3.0 * vtkm::Pif() * rs * rs * rs) + 1;
  Id rc_num =  6* rho_system * vtkm::Ceil(4.0 / 3.0 * vtkm::Pif() * rc * rc * rc) + 1;
  Id rcs_num = rc_num - rs_num;
  Real random_num = GetParameter<Real>(PARA_RANDOM_NUM);
  Real random_rate = random_num / (rcs_num);
  Id pice_num = std::ceil(1.0 / random_rate);
  //std::cout << "rs_num = " << rs_num << " rc_num = " << rc_num << " rcs_num = " << rcs_num
  //          << " random_num = " << random_num << " random_rate = " << random_rate << " pice_num = " << pice_num  << std::endl;
  //Id cutoff_num = rs_num + random_num;
  Id cutoff_num = rc_num + random_num;
  Id verletlist_num = N * cutoff_num;

  num_verletlist.Allocate(N * 2);
  id_verletlist.Allocate(verletlist_num);
  offset_verletlist.Allocate(verletlist_num);

  std::vector<Id> stride_vec(N + 1);
  Id inc = 0;
  std::generate(
    stride_vec.begin(), stride_vec.end(), [&](void) -> Id { return (inc++) * cutoff_num; });
  vtkm::cont::ArrayHandle<vtkm::Id> stride_array = vtkm::cont::make_ArrayHandle(stride_vec);

  auto id_verletlist_group =
    vtkm::cont::make_ArrayHandleGroupVecVariable(id_verletlist, stride_array);
  auto offset_verletlist_group =
    vtkm::cont::make_ArrayHandleGroupVecVariable(offset_verletlist, stride_array);
  auto num_verletlist_group = vtkm::cont::make_ArrayHandleGroupVec<2>(num_verletlist);

  SystemWorklet::ComputeRBLNeighboursOnce(rs_num,
                                          pice_num,
                                          _atoms_id,
                                          _locator,
                                          id_verletlist_group,
                                          num_verletlist_group,
                                          offset_verletlist_group);

  auto special_offsets = GetFieldAsArrayHandle<Id>(field::special_offsets);
  auto special_weights = GetFieldAsArrayHandle<Real>(field::special_weights);
  auto specoal_ids = GetFieldAsArrayHandle<Id>(field::special_ids);
  auto ids_group = vtkm::cont::make_ArrayHandleGroupVecVariable(specoal_ids, special_offsets);
  auto weight_group =
    vtkm::cont::make_ArrayHandleGroupVecVariable(special_weights, special_offsets);
  
  
  _EleNearPairtimer.Start();
  if (_use_erf == true)
  {

    /*SystemWorklet::NearForceRBLERF(rs_num,
                                   pice_num,
                                   _unit_factor._qqr2e,
                                   _atoms_id,
                                   _locator,
                                   _topology,
                                   _force_function,
                                   _static_table,
                                   id_verletlist_group,
                                   num_verletlist_group,
                                   offset_verletlist_group,
                                   corr_force);*/
    SystemWorklet::NearForceRBLERFSpecialBonds(rs_num,
                                   pice_num,
                                   _unit_factor._qqr2e,
                                   _atoms_id,
                                   _locator,
                                   _topology,
                                   _force_function,
                                   _static_table,
                                   id_verletlist_group,
                                   num_verletlist_group,
                                   offset_verletlist_group,
                                   ids_group,
                                   weight_group,
                                   corr_force);
  }
  else
  {
    SystemWorklet::NearForceRBL(rs_num,
                                pice_num,
                                _unit_factor._qqr2e,
                                _atoms_id,
                                _locator,
                                _topology,
                                _force_function,
                                id_verletlist_group,
                                num_verletlist_group,
                                offset_verletlist_group,
                                corr_force);
  }
  
  
  Vec3f corr_value = vtkm::cont::Algorithm::Reduce(corr_force, vtkm::TypeTraits<Vec3f>::ZeroInitialization()) / N;
  SystemWorklet::SumRBLCorrForce(corr_value, corr_force, nearforce);
  _EleNearPairtimer_counting =
    _EleNearPairtimer_counting + _EleNearPairtimer.GetElapsedTime();
  std::cout << "RBL pair time: " << _EleNearPairtimer_counting << std::endl;

  //auto N = _position.GetNumberOfValues();
  //vtkm::cont::ArrayHandle<Vec3f> corr_force;
  //corr_force.Allocate(N);
  //ComputeCorrForce(corr_force);
  //Vec3f corr_value =
  //  vtkm::cont::Algorithm::Reduce(corr_force, vtkm::TypeTraits<Vec3f>::ZeroInitialization()) / N;
  //SystemWorklet::SumRBLCorrForce(corr_value, corr_force, nearforce);
}

void MDSystem::ComputeRBLLJForce(ArrayHandle<Vec3f>& LJforce)
{
  auto N = _position.GetNumberOfValues();
  vtkm::cont::ArrayHandle<Vec3f> corr_ljforce;
  corr_ljforce.Allocate(N);

  auto rc = GetParameter<Real>(PARA_CUTOFF);
  auto rs = GetParameter<Real>(PARA_RS);
  auto rho_system = GetParameter<Real>(PARA_RHO);

  vtkm::cont::ArrayHandle<vtkm::Id> num_verletlist;
  vtkm::cont::ArrayHandle<vtkm::Id> id_verletlist;
  vtkm::cont::ArrayHandle<vtkm::Vec3f> offset_verletlist;

  Id rs_num = 6*rho_system * vtkm::Ceil(4.0 / 3.0 * vtkm::Pif() * rs * rs * rs) + 1;
  Id rc_num = 6*rho_system * vtkm::Ceil(4.0 / 3.0 * vtkm::Pif() * rc * rc * rc) + 1;
  Id rcs_num = rc_num - rs_num;
  Real random_num = GetParameter<Real>(PARA_RANDOM_NUM);
  Real random_rate = random_num / (rcs_num);
  Id pice_num = std::ceil(1.0 / random_rate);
  //Id cutoff_num = rs_num + random_num;
  Id cutoff_num = rc_num + random_num;
  Id verletlist_num = N * cutoff_num;

  num_verletlist.Allocate(N * 2);
  id_verletlist.Allocate(verletlist_num);
  offset_verletlist.Allocate(verletlist_num);

  std::vector<Id> stride_vec(N + 1);
  Id inc = 0;
  std::generate(
    stride_vec.begin(), stride_vec.end(), [&](void) -> Id { return (inc++) * cutoff_num; });
  vtkm::cont::ArrayHandle<vtkm::Id> stride_array = vtkm::cont::make_ArrayHandle(stride_vec);

  auto id_verletlist_group =
    vtkm::cont::make_ArrayHandleGroupVecVariable(id_verletlist, stride_array);
  auto offset_verletlist_group =
    vtkm::cont::make_ArrayHandleGroupVecVariable(offset_verletlist, stride_array);
  auto num_verletlist_group = vtkm::cont::make_ArrayHandleGroupVec<2>(num_verletlist);

  SystemWorklet::ComputeRBLNeighboursOnce(rs_num,
                                          pice_num,
                                          _atoms_id,
                                          _locator,
                                          id_verletlist_group,
                                          num_verletlist_group,
                                          offset_verletlist_group);

  SystemWorklet::LJForceRBL(rs_num,
                            pice_num,
                            _atoms_id,
                            _locator,
                            _topology,
                            _force_function,
                            id_verletlist_group,
                            num_verletlist_group,
                            offset_verletlist_group,
                            corr_ljforce);

  Vec3f corr_value =
    vtkm::cont::Algorithm::Reduce(corr_ljforce, vtkm::TypeTraits<Vec3f>::ZeroInitialization()) / N;
  SystemWorklet::SumRBLCorrForce(corr_value, corr_ljforce, LJforce);
}

void MDSystem::ComputeVerletlistNearForce(ArrayHandle<Vec3f>& nearforce)
{
  auto cut_off = GetParameter<Real>(PARA_CUTOFF);
  auto N = _position.GetNumberOfValues();
  auto rho_system = GetParameter<Real>(PARA_RHO);
  auto max_j_num = rho_system * vtkm::Ceil(4.0 / 3.0 * vtkm::Pif() * cut_off * cut_off * cut_off) + 1;
  auto verletlist_num = N * max_j_num;

  ArrayHandle<Id> id_verletlist;
  ArrayHandle<Vec3f> offset_verletlist;
  offset_verletlist.Allocate(verletlist_num);
  id_verletlist.Allocate(verletlist_num);

  std::vector<Id> temp_vec(N + 1);
  Id inc = 0;
  std::generate(temp_vec.begin(), temp_vec.end(), [&](void) -> Id { return (inc++) * max_j_num; });
  vtkm::cont::ArrayHandle<vtkm::Id> temp_offset = vtkm::cont::make_ArrayHandle(temp_vec);

  auto id_verletlist_group = vtkm::cont::make_ArrayHandleGroupVecVariable(id_verletlist, temp_offset);
  auto offset_verletlist_group = vtkm::cont::make_ArrayHandleGroupVecVariable(offset_verletlist, temp_offset);
  vtkm::cont::ArrayHandle<vtkm::Id> num_verletlist;

  SystemWorklet::ComputeNeighbours( cut_off, _atoms_id, _locator, id_verletlist_group, num_verletlist, offset_verletlist_group);

  if (_use_erf == true)
  {
    SystemWorklet::NearForceVerletERF(cut_off,
                                      _atoms_id,
                                      _locator,
                                      _topology,
                                      _force_function,
                                      _static_table,
                                      id_verletlist_group,
                                      num_verletlist,
                                      offset_verletlist_group,
                                      nearforce);
  }
  else
  {
    SystemWorklet::NearForceVerlet(cut_off,
                                   _atoms_id,
                                   _locator,
                                   _topology,
                                   _force_function,
                                   id_verletlist_group,
                                   num_verletlist,
                                   offset_verletlist_group,
                                   nearforce);
  }
}

void MDSystem::ComputeVerletlistLJForce(ArrayHandle<Vec3f>& ljforce)
{
  auto cut_off = GetParameter<Real>(PARA_CUTOFF);
  auto N = _position.GetNumberOfValues();
  auto rho_system = GetParameter<Real>(PARA_RHO);
  auto max_j_num = rho_system * vtkm::Ceil(4.0 / 3.0 * vtkm::Pif() * cut_off * cut_off * cut_off) + 1;
  auto verletlist_num = N * max_j_num;

  ArrayHandle<Id> id_verletlist;
  ArrayHandle<Vec3f> offset_verletlist;
  offset_verletlist.Allocate(verletlist_num);
  id_verletlist.Allocate(verletlist_num);

  std::vector<Id> temp_vec(N + 1);
  Id inc = 0;
  std::generate(temp_vec.begin(), temp_vec.end(), [&](void) -> Id { return (inc++) * max_j_num; });
  vtkm::cont::ArrayHandle<vtkm::Id> temp_offset = vtkm::cont::make_ArrayHandle(temp_vec);

  vtkm::cont::ArrayHandle<vtkm::Id> num_verletlist;
  auto id_verletlist_group = vtkm::cont::make_ArrayHandleGroupVecVariable(id_verletlist, temp_offset);
  auto offset_verletlist_group = vtkm::cont::make_ArrayHandleGroupVecVariable(offset_verletlist, temp_offset);

  SystemWorklet::ComputeNeighbours( cut_off, _atoms_id, _locator, id_verletlist_group, num_verletlist, offset_verletlist_group);

  SystemWorklet::LJForceVerlet(cut_off,
                               _atoms_id,
                               _locator,
                               _topology,
                               _force_function,
                               id_verletlist_group,
                               num_verletlist,
                               offset_verletlist_group,
                               ljforce);
}
void MDSystem::ComputeOriginalLJForce(ArrayHandle<Vec3f>& ljforce)
{
  auto cut_off = GetParameter<Real>(PARA_CUTOFF);

  SystemWorklet::LJForceWithPeriodicBC(
    cut_off, _atoms_id, _locator, _topology, _force_function, ljforce);
}

void MDSystem::ComputeRBLEAMForce(ArrayHandle<Vec3f>& force)
{
  ArrayHandle<Real> fp;
  auto Vlength = GetParameter<Real>(PARA_VLENGTH);
  auto rhor_spline = GetFieldAsArrayHandle<Vec7f>(field::rhor_spline);
  auto frho_spline = GetFieldAsArrayHandle<Vec7f>(field::frho_spline);
  auto z2r_spline = GetFieldAsArrayHandle<Vec7f>(field::z2r_spline);
  auto N = _position.GetNumberOfValues();

  auto rc = GetParameter<Real>(PARA_CUTOFF);
  auto rs = GetParameter<Real>(PARA_RS);
  auto rho_system = GetParameter<Real>(PARA_RHO);
  vtkm::cont::ArrayHandle<Vec3f> corr_force;
  corr_force.Allocate(N);

  vtkm::cont::ArrayHandle<vtkm::Id> num_verletlist;
  vtkm::cont::ArrayHandle<vtkm::Id> id_verletlist;
  vtkm::cont::ArrayHandle<vtkm::Vec3f> offset_verletlist;

  Id rs_num = rho_system * vtkm::Ceil(4.0 / 3.0 * vtkm::Pif() * rs * rs * rs) + 1;
  Id rc_num = rho_system * vtkm::Ceil(4.0 / 3.0 * vtkm::Pif() * rc * rc * rc) + 1;
  Id rcs_num = rc_num - rs_num;
  Real random_num = GetParameter<Real>(PARA_RANDOM_NUM);
  Real random_rate = random_num / (rcs_num);
  Id pice_num = std::ceil(1.0 / random_rate);
  //Id cutoff_num = rs_num + random_num;
  Id cutoff_num = rc_num + random_num;
  Id verletlist_num = N * cutoff_num;

  num_verletlist.Allocate(N * 2);
  id_verletlist.Allocate(verletlist_num);
  offset_verletlist.Allocate(verletlist_num);

  std::vector<Id> stride_vec(N + 1);
  Id inc = 0;
  std::generate(
    stride_vec.begin(), stride_vec.end(), [&](void) -> Id { return (inc++) * cutoff_num; });
  vtkm::cont::ArrayHandle<vtkm::Id> stride_array = vtkm::cont::make_ArrayHandle(stride_vec);

  auto id_verletlist_group =
    vtkm::cont::make_ArrayHandleGroupVecVariable(id_verletlist, stride_array);
  auto offset_verletlist_group =
    vtkm::cont::make_ArrayHandleGroupVecVariable(offset_verletlist, stride_array);
  auto num_verletlist_group = vtkm::cont::make_ArrayHandleGroupVec<2>(num_verletlist);

  SystemWorklet::ComputeRBLNeighboursOnce(rs_num,
                                          pice_num,
                                          _atoms_id,
                                          _locator,
                                          id_verletlist_group,
                                          num_verletlist_group,
                                          offset_verletlist_group);

  SystemWorklet::EAMfp(rc,
                       Vlength,
                       rs,
                       rs_num,
                       pice_num,
                       _atoms_id,
                       rhor_spline,
                       frho_spline,
                       _locator,
                       _force_function,
                       id_verletlist_group,
                       num_verletlist_group,
                       offset_verletlist_group,
                       fp);

 SystemWorklet::EAMRBLForce(rc,
                            Vlength,
                            rs,
                            rs_num,                            
                            pice_num,
                            _atoms_id,
                            rhor_spline,
                            z2r_spline,
                            fp,
                            _locator,
                            _force_function,
                            id_verletlist_group,
                            num_verletlist_group,
                            offset_verletlist_group,
                            corr_force);

  Vec3f corr_value =
    vtkm::cont::Algorithm::Reduce(corr_force, vtkm::TypeTraits<Vec3f>::ZeroInitialization()) / N;
 SystemWorklet::SumRBLCorrForce(corr_value, corr_force, force);
}

void MDSystem::ComputeVerletlistEAMForce(ArrayHandle<Vec3f>& force)
{
 ArrayHandle<Real> fp;
 auto Vlength = GetParameter<Real>(PARA_VLENGTH);
 auto rhor_spline = GetFieldAsArrayHandle<Vec7f>(field::rhor_spline);
 auto frho_spline = GetFieldAsArrayHandle<Vec7f>(field::frho_spline);
 auto z2r_spline = GetFieldAsArrayHandle<Vec7f>(field::z2r_spline);

  auto cut_off = GetParameter<Real>(PARA_CUTOFF);
  auto N = _position.GetNumberOfValues();
  auto rho_system = GetParameter<Real>(PARA_RHO);
  auto max_j_num = rho_system * vtkm::Ceil(4.0 / 3.0 * vtkm::Pif() * cut_off * cut_off * cut_off) + 1;
  auto verletlist_num = N * max_j_num;

  ArrayHandle<Id> id_verletlist;
  ArrayHandle<Vec3f> offset_verletlist;
  offset_verletlist.Allocate(verletlist_num);
  id_verletlist.Allocate(verletlist_num);

  std::vector<Id> temp_vec(N + 1);
  Id inc = 0;
  std::generate(temp_vec.begin(), temp_vec.end(), [&](void) -> Id { return (inc++) * max_j_num; });
  vtkm::cont::ArrayHandle<vtkm::Id> temp_offset = vtkm::cont::make_ArrayHandle(temp_vec);

  vtkm::cont::ArrayHandle<vtkm::Id> num_verletlist;
  auto id_verletlist_group =
    vtkm::cont::make_ArrayHandleGroupVecVariable(id_verletlist, temp_offset);
  auto offset_verletlist_group =
    vtkm::cont::make_ArrayHandleGroupVecVariable(offset_verletlist, temp_offset);

  SystemWorklet::ComputeNeighbours(
    cut_off, _atoms_id, _locator, id_verletlist_group, num_verletlist, offset_verletlist_group);

  SystemWorklet::EAMfpVerlet(cut_off,
                             Vlength,
                             _atoms_id,
                             rhor_spline,
                             frho_spline,
                             _locator,
                             _force_function,
                             id_verletlist_group,
                             num_verletlist,
                             offset_verletlist_group,
                             fp);

  SystemWorklet::EAMForceVerlet(cut_off,
                                Vlength,
                               _atoms_id,
                                rhor_spline,
                                z2r_spline,
                                fp,
                               _locator,
                               _force_function,
                               id_verletlist_group,
                               num_verletlist,
                               offset_verletlist_group,
                               force);
}

void MDSystem::ComputeOriginalEAMForce(ArrayHandle<Vec3f>& force)
{
  auto cut_off = GetParameter<Real>(EAM_PARA_CUTOFF);
  auto Vlength = GetParameter<Real>(PARA_VLENGTH);

  auto rhor_spline = GetFieldAsArrayHandle<Vec7f>(field::rhor_spline);
  auto frho_spline = GetFieldAsArrayHandle<Vec7f>(field::frho_spline);
  auto z2r_spline = GetFieldAsArrayHandle<Vec7f>(field::z2r_spline);

  ArrayHandle<Real> EAM_rho;
  ArrayHandle<Real> fp;

  //1:compute _EAM_rho   = density at each atom
  SystemWorklet::EAM_rho(
    cut_off, Vlength, _atoms_id, rhor_spline, _locator, _topology, _force_function, EAM_rho);

  // 2:compute fp    = derivative of embedding energy at each atom
  SystemWorklet::EAM_fp(_atoms_id, EAM_rho, frho_spline, _locator, _topology, _force_function, fp);

  // 3:compute force  = EAM_force
  SystemWorklet::EAM_force(cut_off,
                           Vlength,
                           _atoms_id,
                           fp,
                           rhor_spline,
                           z2r_spline,
                           _locator,
                           _topology,
                           _force_function,
                           force);
}

void MDSystem::ComputeSpecialBondsLJForce(ArrayHandle<Vec3f>& ljforce) 
{
  auto cut_off = GetParameter<Real>(PARA_CUTOFF);

  auto special_offsets = GetFieldAsArrayHandle<Id>(field::special_offsets);
  auto special_weights = GetFieldAsArrayHandle<Real>(field::special_weights);
  auto specoal_ids = GetFieldAsArrayHandle<Id>(field::special_ids);
  auto ids_group = vtkm::cont::make_ArrayHandleGroupVecVariable(specoal_ids, special_offsets);
  auto weight_group = vtkm::cont::make_ArrayHandleGroupVecVariable(special_weights, special_offsets);

  SystemWorklet::SpecicalBondsLJForce(
    cut_off, _atoms_id, _locator, _topology, _force_function, ids_group, weight_group, ljforce);
}
