#include "FieldName.h"
#include "ExecutionMD.h"
#include "run/worklet/MolecularWorklet.h"
#include "run/worklet/RunWorklet.h"
#include <vtkm/cont/Algorithm.h>
#include <vtkm/cont/ArrayHandleGroupVec.h>
#include "Executioner.h"
#include "ERFTable.h"
#include <vtkm/worklet/Keys.h>
#include <vtkm/cont/ArrayHandleConstant.h>
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
  using ExecutionSignature = void(_1, _2);

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

ExecutionMD::ExecutionMD(const Configuration& cfg)
  : Execution(cfg)
{

}

void ExecutionMD::Init()
{
  _position = _para.GetFieldAsArrayHandle<Vec3f>(field::position);
  InitPointLocator(); 
  SetForceFunction();
  SetTopology();

  if (_para.GetParameter<bool>(PARA_FAR_FORCE))
  {
    auto N = _position.GetNumberOfValues();
    _rhok_Re.AllocateAndFill(_RBE_P * N, 0.0);
    _rhok_Im.AllocateAndFill(_RBE_P * N, 0.0);

    vtkm::cont::ArrayHandleIndex indexArray(N * _RBE_P);
    vtkm::cont::Invoker{}(SetIndex{ N }, indexArray, _psamplekey);
  }
  _EleNearPairtimer_counting = 0.0;
}

void ExecutionMD::Execute() 
{
  _timer.Start();
  PreSolve();
  Solve();
  PostSolve();
}

void ExecutionMD::SetParameters()
{
  _para.SetParameter(PARA_ENSEMBLE, Get<std::string>("ensemble"));
  _para.SetParameter(PARA_TEMP_CTRL_TYPE, Get<std::string>("temp_ctrl_type"));
  _para.SetParameter(PARA_PRESS_CTRL_TYPE, Get<std::string>("press_ctrl_type"));
  _para.SetParameter(PARA_TIMESTEP, Get<Real>("timestep"));
  _para.SetParameter(PARA_NUM_STEPS, Get<Real>("num_steps"));
  _para.SetParameter(PARA_TEMPERATURE, GetVectorValue<Real>("temperature"));
  _para.SetParameter(PARA_PRESSURE, GetVectorValue<Real>("pressure"));
}

void ExecutionMD::InitialCondition()
{
  // 物理场初始化
  if (_init_way == "inbuild")
  {
    if (_para.GetParameter<bool>(PARA_FAR_FORCE))
    {
      _charge = _para.GetFieldAsArrayHandle<Real>(field::charge);
      auto n = _position.GetNumberOfValues();
      _charge.Allocate(n);
      _charge.Fill(-1.0, 0);
      _charge.Fill(1.0, n / 2);
      _topology.SetCharge(_charge);
    }
    else
    {
      _charge = _para.GetFieldAsArrayHandle<Real>(field::charge);
      auto n = _position.GetNumberOfValues();
      _charge.AllocateAndFill(n, 0);
      _topology.SetCharge(_charge);
    }
  }
  else if (_para.GetParameter<std::string>(PARA_FILE_TYPE) == "EAM")
  {
    _charge = _para.GetFieldAsArrayHandle<Real>(field::charge);
    auto n = _position.GetNumberOfValues();
    _charge.AllocateAndFill(n, 0);
    _topology.SetCharge(_charge);
  }
  else
  {
    _charge = _para.GetFieldAsArrayHandle<Real>(field::charge);
    _topology.SetCharge(_charge);
  }

  _velocity = _para.GetFieldAsArrayHandle<Vec3f>(field::velocity);
  _mass = _para.GetFieldAsArrayHandle<Real>(field::mass);

  _molecule_id = _para.GetFieldAsArrayHandle<Id>(field::molecule_id);
  _atoms_id = _para.GetFieldAsArrayHandle<Id>(field::atom_id);
}

void ExecutionMD::SetForceFunction() {}

void ExecutionMD::SetTopology() {}

void ExecutionMD::InitERF() 
{
  if (_para.GetParameter<bool>(PARA_FAR_FORCE)/*_use_erf == true*/)
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

//void ExecutionMD::InitField()
//{
//  _para.AddField(field::charge, ArrayHandle<Real>{});
//  _para.AddField(field::velocity, ArrayHandle<Vec3f>{});
//  _para.AddField(field::mass, ArrayHandle<Real>{});
//  _para.AddField(field::molecule_id, ArrayHandle<Id>{});
//  _para.AddField(field::atom_id, ArrayHandle<Id>{});
//  _para.AddField(field::position, ArrayHandle<Vec3f>{});
//  _para.AddField(field::special_offsets, ArrayHandle<Id>{});
//  _para.AddField(field::special_weights, ArrayHandle<Real>{});
//  _para.AddField(field::special_ids, ArrayHandle<Id>{});
//}

void ExecutionMD::UpdateVerletList() 
{
  auto cut_off = _para.GetParameter<Real>(PARA_CUTOFF);
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

  RunWorklet::ComputeNeighbours(
    cut_off, _atoms_id, _locator, id_verletlist_group, num_verletlist, offset_verletlist_group);

  _locator.SetVerletListInfo(num_verletlist, id_verletlist_group, offset_verletlist_group);
}

void ExecutionMD::ComputeCorrForce(vtkm::cont::ArrayHandle<Vec3f>& corr_force)
{
  auto rc = _para.GetParameter<Real>(PARA_CUTOFF);
  auto rs = _para.GetParameter<Real>(PARA_R_CORE);
  auto rho_system = _para.GetParameter<Real>(PARA_RHO);
  auto N = _position.GetNumberOfValues();

  vtkm::cont::ArrayHandle<vtkm::Id> num_verletlist;
  vtkm::cont::ArrayHandle<vtkm::Id> id_verletlist;
  vtkm::cont::ArrayHandle<vtkm::Vec3f> offset_verletlist;

  Id rs_num = rho_system * vtkm::Ceil(4.0 / 3.0 * vtkm::Pif() * rs * rs * rs) + 1;
  Id rc_num = rho_system * vtkm::Ceil(4.0 / 3.0 * vtkm::Pif() * rc * rc * rc) + 1;
  Id rcs_num = rc_num - rs_num;
  Real random_num = _para.GetParameter<Real>(PARA_NEIGHBOR_SAMPLE_NUM);
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

  RunWorklet::ComputeRBLNeighboursOnce(rs_num,
                                          pice_num,
                                          _atoms_id,
                                          _locator,
                                          id_verletlist_group,
                                          num_verletlist_group,
                                          offset_verletlist_group);

  RunWorklet::NearForceRBLERF(rs_num,
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

std::vector<Vec2f> ExecutionMD::ComputeChargeStructureFactorRBE(Real& _Vlength, ArrayHandle<Vec3f>& _psample)
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
    RunWorklet::ComputeChargeStructureFactorComponent(
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

void ExecutionMD::InitParameters() 
{
  _unit = _para.GetParameter<std::string>(PARA_UNIT);
  if (_unit == "REAL")
  {
    _unit_factor._kB = 1.9872067 * vtkm::Pow(10.0, -3);
    _unit_factor._fmt2v = 4.186 * vtkm::Pow(10.0, -4);
    _unit_factor._mvv2e = 1.0 / (4.186 * vtkm::Pow(10.0, -4));
    _unit_factor._qqr2e = 332.06371;
    _unit_factor._boltz = 0.0019872067;
    _unit_factor._nktv2p = 68568.415;
  }
  else if (_unit == "LJ")
  {
    _unit_factor._kB = 1.0;
    _unit_factor._fmt2v = 1.0;
    _unit_factor._mvv2e = 1.0;
    _unit_factor._qqr2e = 1.0;
    _unit_factor._boltz = 1.0;
    _unit_factor._nktv2p = 1.0;
  }
  else if (_unit == "METAL")
  {
    _unit_factor._kB = 8.617343e-5;
    _unit_factor._fmt2v = 1.0 / 1.0364269e-4;
    _unit_factor._mvv2e = 1.0364269e-4;
    _unit_factor._qqr2e = 14.399645;
    _unit_factor._boltz = 8.617343e-5;
    _unit_factor._nktv2p = 1.6021765e6;
  }
  _para.SetParameter(PARA_UNIT_FACTOR, _unit_factor);
}

std::vector<Vec2f> ExecutionMD::ComputeChargeStructureFactorEwald(Real& _Vlength,
                                                                      IdComponent& Kmax)
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
          RunWorklet::ComputeChargeStructureFactorComponent(
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

void ExecutionMD::ComputeRBEEleForce(ArrayHandle<Vec3f>& psample,
                                  IdComponent& RBE_P,
                                  ArrayHandle<Vec3f>& RBE_ele_force)
{
 

  auto Vlength = _para.GetParameter<Real>(PARA_VLENGTH);
 

  ArrayHandle<Vec2f> new_whole_rhok;
  ComputeNewChargeStructureFactorRBE(Vlength, psample, new_whole_rhok);
  
  RunWorklet::ComputeNewRBEForce(RBE_P,
                                    _atoms_id,
                                    psample,
                                    /*whole_rhok,*/ new_whole_rhok, _force_function,
                                    _topology,
                                    _locator,
                                    RBE_ele_force);
 

}

void  ExecutionMD::ComputeNewChargeStructureFactorRBE(Real& _Vlength,
                                                                ArrayHandle<Vec3f>& _psample,
                                                                ArrayHandle<Vec2f>& new_rhok)
{
 
  auto N = _position.GetNumberOfValues();
  
  const Id p_number = _psample.GetNumberOfValues();
  


  RunWorklet::ComputePnumberChargeStructureFactor(_Vlength,
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
 


  RunWorklet::ChangePnumberChargeStructureFactor(rhok_Re_reduce, rhok_Im_reduce, new_rhok);
 

  
}

void ExecutionMD::ComputeEwaldEleForce(IdComponent& Kmax, ArrayHandle<Vec3f>& Ewald_ele_force)
{
  auto Vlength = _para.GetParameter<Real>(PARA_VLENGTH);
  auto rhok = ComputeChargeStructureFactorEwald(Vlength, Kmax);
  ArrayHandle<Vec2f> whole_rhok = vtkm::cont::make_ArrayHandle(rhok);
  RunWorklet::ComputeNewFarElectrostatics(
    Kmax, _atoms_id, whole_rhok, _force_function, _topology, _locator, Ewald_ele_force);
}

void ExecutionMD::ComputeRBLNearForce(ArrayHandle<Vec3f>& nearforce)
{
  auto N = _position.GetNumberOfValues();
  vtkm::cont::ArrayHandle<Vec3f> corr_force;
  corr_force.Allocate(N);

  auto rc = _para.GetParameter<Real>(PARA_CUTOFF);
  auto rs = _para.GetParameter<Real>(PARA_R_CORE);
  auto rho_system = _para.GetParameter<Real>(PARA_RHO);

  vtkm::cont::ArrayHandle<vtkm::Id> num_verletlist;
  vtkm::cont::ArrayHandle<vtkm::Id> id_verletlist;
  vtkm::cont::ArrayHandle<vtkm::Vec3f> offset_verletlist;

  //0.05 is a temp parameter, to fit with the low density of system requirement
  //so as the RBL for LJ

  Real coeff_rcs = 1.0 + (0.05 / rho_system - 0.05);
  Id Id_coeff_rcs = vtkm::Round(coeff_rcs);
  Id rs_num = Id_coeff_rcs * rho_system * vtkm::Ceil(4.0 / 3.0 * vtkm::Pif() * rs * rs * rs) + 1;
  Id rc_num = Id_coeff_rcs * rho_system * vtkm::Ceil(4.0 / 3.0 * vtkm::Pif() * rc * rc * rc) + 1;
  Id rcs_num = rc_num - rs_num;
  Real random_num = _para.GetParameter<Real>(PARA_NEIGHBOR_SAMPLE_NUM);
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

  RunWorklet::ComputeRBLNeighboursOnce(rs_num,
                                       pice_num,
                                       _atoms_id,
                                       _locator,
                                       id_verletlist_group,
                                       num_verletlist_group,
                                       offset_verletlist_group);

  auto special_offsets = _para.GetFieldAsArrayHandle<Id>(field::special_offsets);
  auto special_weights = _para.GetFieldAsArrayHandle<Real>(field::special_weights);
  auto specoal_ids = _para.GetFieldAsArrayHandle<Id>(field::special_ids);
  auto ids_group = vtkm::cont::make_ArrayHandleGroupVecVariable(specoal_ids, special_offsets);
  auto weight_group =
    vtkm::cont::make_ArrayHandleGroupVecVariable(special_weights, special_offsets);
  
  
  _EleNearPairtimer.Start();

  // 当前默认使用ERF 注意：这里要看NearForceRBLERFSpecialBonds和NearForceRBLERF
  //if (_use_erf == true)
  //{
   // RunWorklet::NearForceRBLERF(rs_num,
   //                                pice_num,
   //                                _unit_factor._qqr2e,
   //                                _atoms_id,
   //                                _locator,
   //                                _topology,
   //                                _force_function,
   //                                _static_table,
   //                                id_verletlist_group,
   //                                num_verletlist_group,
   //                                offset_verletlist_group,
   //                                corr_force);
  // 这里的判断待定！！
  if (_para.GetParameter<bool>(PARA_DIHEDRALS_FORCE))
  {
    RunWorklet::NearForceRBLERFSpecialBonds(rs_num,
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
    RunWorklet::NearForceRBLERF(rs_num,
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

  //}
  //else
  //{
  //  RunWorklet::NearForceRBL(rs_num,
  //                              pice_num,
  //                              _unit_factor._qqr2e,
  //                              _atoms_id,
  //                              _locator,
  //                              _topology,
  //                              _force_function,
  //                              id_verletlist_group,
  //                              num_verletlist_group,
  //                              offset_verletlist_group,
  //                              corr_force);
  //}

  Vec3f corr_value = vtkm::cont::Algorithm::Reduce(corr_force, vtkm::TypeTraits<Vec3f>::ZeroInitialization()) / N;
  RunWorklet::SumRBLCorrForce(corr_value, corr_force, nearforce);

  //auto N = _position.GetNumberOfValues();
  //vtkm::cont::ArrayHandle<Vec3f> corr_force;
  //corr_force.Allocate(N);
  //ComputeCorrForce(corr_force);
  //Vec3f corr_value =
  //  vtkm::cont::Algorithm::Reduce(corr_force, vtkm::TypeTraits<Vec3f>::ZeroInitialization()) / N;
  //RunWorklet::SumRBLCorrForce(corr_value, corr_force, nearforce);
}

void ExecutionMD::ComputeRBLLJForce(ArrayHandle<Vec3f>& LJforce)
{
  auto N = _position.GetNumberOfValues();
  vtkm::cont::ArrayHandle<Vec3f> corr_ljforce;
  corr_ljforce.Allocate(N);

  auto rc = _para.GetParameter<Real>(PARA_CUTOFF);
  auto rs = _para.GetParameter<Real>(PARA_R_CORE);
  auto range = _para.GetParameter<vtkm::Vec<vtkm::Range, 3>>(PARA_RANGE);
  auto vol = (range[0].Max - range[0].Min) * (range[1].Max - range[1].Min) * (range[2].Max - range[2].Min);
  auto rho_system = N / vol;
  //auto rho_system = _para.GetParameter<Real>(PARA_RHO);

  vtkm::cont::ArrayHandle<vtkm::Id> num_verletlist;
  vtkm::cont::ArrayHandle<vtkm::Id> id_verletlist;
  vtkm::cont::ArrayHandle<vtkm::Vec3f> offset_verletlist;

  //0.05 is a temp parameter, to fit with the low density of system requirement
  //so as the RBL for LJ
  Real coeff_rcs = 1.0 + (0.05 / rho_system - 0.05);
  Id Id_coeff_rcs = vtkm::Round(coeff_rcs);
  Id rs_num = Id_coeff_rcs * rho_system * vtkm::Ceil(4.0 / 3.0 * vtkm::Pif() * rs * rs * rs) + 1;
  Id rc_num = Id_coeff_rcs * rho_system * vtkm::Ceil(4.0 / 3.0 * vtkm::Pif() * rc * rc * rc) + 1;
  Id rcs_num = rc_num - rs_num;
  Real random_num = _para.GetParameter<Real>(PARA_NEIGHBOR_SAMPLE_NUM);
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

  RunWorklet::ComputeRBLNeighboursOnce(rs_num,
                                          pice_num,
                                          _atoms_id,
                                          _locator,
                                          id_verletlist_group,
                                          num_verletlist_group,
                                          offset_verletlist_group);

  RunWorklet::LJForceRBL(rs_num,
                            pice_num,
                            _atoms_id,
                            _locator,
                            _topology,
                            _force_function,
                            id_verletlist_group,
                            num_verletlist_group,
                            offset_verletlist_group,
                            corr_ljforce);

    //pbc
  //auto range = _para.GetParameter<vtkm::Vec<vtkm::Range, 3>>(PARA_RANGE);
  //Vec3f Vlength;
  //for (int i = 0; i < 3; ++i)
  //{
  //  Vlength[i] = range[i].Max - range[i].Min;
  //}
  //RunWorklet::LJForceRBLPBC(rs_num,
  //                         pice_num,
  //                          Vlength,
  //                       _atoms_id,
  //                       _locator,
  //                       _topology,
  //                       _force_function,
  //                       id_verletlist_group,
  //                       num_verletlist_group,
  //                       offset_verletlist_group,
  //                       corr_ljforce);

  Vec3f corr_value =
    vtkm::cont::Algorithm::Reduce(corr_ljforce, vtkm::TypeTraits<Vec3f>::ZeroInitialization()) / N;
  RunWorklet::SumRBLCorrForce(corr_value, corr_ljforce, LJforce);
}

void ExecutionMD::ComputeRBLLJForce71(ArrayHandle<Vec3f>& LJforce,
                                      ArrayHandle<Vec6f>& LJvirial,
                                      Real& LJenergy)
{
  auto N = _position.GetNumberOfValues();
  vtkm::cont::ArrayHandle<Vec3f> corr_ljforce;
  vtkm::cont::ArrayHandle<Vec6f> corr_ljvirial;
  vtkm::cont::ArrayHandle<Real> corr_ljenergy;
  vtkm::cont::ArrayHandle<Real> ljenergy_atom;
  corr_ljforce.Allocate(N);

  auto rc = _para.GetParameter<Real>(PARA_CUTOFF);
  auto rs = _para.GetParameter<Real>(PARA_R_CORE);
  auto range = _para.GetParameter<vtkm::Vec<vtkm::Range, 3>>(PARA_RANGE);
  auto vol =
    (range[0].Max - range[0].Min) * (range[1].Max - range[1].Min) * (range[2].Max - range[2].Min);
  auto rho_system = N / vol;
  //auto rho_system = _para.GetParameter<Real>(PARA_RHO);

  vtkm::cont::ArrayHandle<vtkm::Id> num_verletlist;
  vtkm::cont::ArrayHandle<vtkm::Id> id_verletlist;
  vtkm::cont::ArrayHandle<vtkm::Vec3f> offset_verletlist;

  //0.05 is a temp parameter, to fit with the low density of system requirement
  //so as the RBL for LJ
  Real coeff_rcs = 1.0 + (0.05 / rho_system - 0.05);
  Id Id_coeff_rcs = vtkm::Round(coeff_rcs);
  Id rs_num = Id_coeff_rcs * rho_system * vtkm::Ceil(4.0 / 3.0 * vtkm::Pif() * rs * rs * rs) + 1;
  Id rc_num = Id_coeff_rcs * rho_system * vtkm::Ceil(4.0 / 3.0 * vtkm::Pif() * rc * rc * rc) + 1;
  Id rcs_num = rc_num - rs_num;
  Real random_num = _para.GetParameter<Real>(PARA_NEIGHBOR_SAMPLE_NUM);
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

  RunWorklet::ComputeRBLNeighboursOnce(rs_num,
                                       pice_num,
                                       _atoms_id,
                                       _locator,
                                       id_verletlist_group,
                                       num_verletlist_group,
                                       offset_verletlist_group);
  //pbc
  Vec3f box{ 0, 0, 0 };
  for (int i = 0; i < 3; ++i)
  {
    box[i] = range[i].Max - range[i].Min;
  }
  RunWorklet::LJForceRBL71(rs_num,
                           pice_num,
                           box,
                           _atoms_id,
                           _locator,
                           _topology,
                           _force_function,
                           id_verletlist_group,
                           num_verletlist_group,
                           offset_verletlist_group,
                           corr_ljforce,
                           corr_ljvirial,
                           corr_ljenergy);

  //pbc
  //auto range = _para.GetParameter<vtkm::Vec<vtkm::Range, 3>>(PARA_RANGE);
  //Vec3f Vlength;
  //for (int i = 0; i < 3; ++i)
  //{
  //  Vlength[i] = range[i].Max - range[i].Min;
  //}
  //RunWorklet::LJForceRBLPBC(rs_num,
  //                         pice_num,
  //                          Vlength,
  //                       _atoms_id,
  //                       _locator,
  //                       _topology,
  //                       _force_function,
  //                       id_verletlist_group,
  //                       num_verletlist_group,
  //                       offset_verletlist_group,
  //                       corr_ljforce);

  Vec3f corr_value =
    vtkm::cont::Algorithm::Reduce(corr_ljforce, vtkm::TypeTraits<Vec3f>::ZeroInitialization()) / N;
  RunWorklet::SumRBLCorrForce(corr_value, corr_ljforce, LJforce);

  Vec6f corr_virial_value =
    vtkm::cont::Algorithm::Reduce(corr_ljvirial, vtkm::TypeTraits<Vec6f>::ZeroInitialization()) / N;
  RunWorklet::SumRBLCorrLJvirial(corr_virial_value, corr_ljvirial, LJvirial);

    Real corr_energy_value =
    vtkm::cont::Algorithm::Reduce(corr_ljenergy, vtkm::TypeTraits<Real>::ZeroInitialization()) / N;
  RunWorklet::SumRBLCorrLJenergy(corr_energy_value, corr_ljenergy, ljenergy_atom);
   LJenergy = vtkm::cont::Algorithm::Reduce(ljenergy_atom, vtkm::TypeTraits<Real>::ZeroInitialization()) / N;

}

void ExecutionMD::ComputeVerletlistNearForce(ArrayHandle<Vec3f>& nearforce)
{
  auto cut_off = _para.GetParameter<Real>(PARA_CUTOFF);
  auto N = _position.GetNumberOfValues();
  auto rho_system = _para.GetParameter<Real>(PARA_RHO);
  auto max_j_num = rho_system * vtkm::Ceil(4.0 / 3.0 * vtkm::Pif() * cut_off * cut_off * cut_off) + 1;
  auto verletlist_num = N * N;

  ArrayHandle<Id> id_verletlist;
  ArrayHandle<Vec3f> offset_verletlist;
  offset_verletlist.Allocate(verletlist_num);
  id_verletlist.Allocate(verletlist_num);

  std::vector<Id> temp_vec(N + 1);
  Id inc = 0;
  std::generate(temp_vec.begin(), temp_vec.end(), [&](void) -> Id { return (inc++) * N; });
  vtkm::cont::ArrayHandle<vtkm::Id> temp_offset = vtkm::cont::make_ArrayHandle(temp_vec);

  auto id_verletlist_group = vtkm::cont::make_ArrayHandleGroupVecVariable(id_verletlist, temp_offset);
  auto offset_verletlist_group = vtkm::cont::make_ArrayHandleGroupVecVariable(offset_verletlist, temp_offset);
  vtkm::cont::ArrayHandle<vtkm::Id> num_verletlist;

  RunWorklet::ComputeNeighbours( cut_off, _atoms_id, _locator, id_verletlist_group, num_verletlist, offset_verletlist_group);

  // 默认使用erf
  //if (_use_erf == true)
  //{
    //RunWorklet::NearForceVerletERF(cut_off,
    //                                  _atoms_id,
    //                                  _locator,
    //                                  _topology,
    //                                  _force_function,
    //                                  _static_table,
    //                                  id_verletlist_group,
    //                                  num_verletlist,
    //                                  offset_verletlist_group,
    //                                  nearforce);
  //}
  //else
  //{
    RunWorklet::NearForceVerlet(cut_off,
                                   _atoms_id,
                                   _locator,
                                   _topology,
                                   _force_function,
                                   id_verletlist_group,
                                   num_verletlist,
                                   offset_verletlist_group,
                                   nearforce);
  //}
}

void ExecutionMD::ComputeVerletlistLJForce(ArrayHandle<Vec3f>& ljforce)
{
  auto cut_off = _para.GetParameter<Real>(PARA_CUTOFF);
  auto N = _position.GetNumberOfValues();
  auto rho_system = _para.GetParameter<Real>(PARA_RHO);
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

  RunWorklet::ComputeNeighbours( cut_off, _atoms_id, _locator, id_verletlist_group, num_verletlist, offset_verletlist_group);

  RunWorklet::LJForceVerlet(cut_off,
                            _atoms_id,
                            _locator,
                            _topology,
                            _force_function,
                            id_verletlist_group,
                            num_verletlist,
                            offset_verletlist_group,
                            ljforce);
}
void ExecutionMD::ComputeOriginalLJForce(ArrayHandle<Vec3f>& ljforce)
{
  auto cut_off = _para.GetParameter<Real>(PARA_CUTOFF);

  RunWorklet::LJForceWithPeriodicBC(
    cut_off, _atoms_id, _locator, _topology, _force_function, ljforce);
}

void ExecutionMD::ComputeRBLEAMForce(ArrayHandle<Vec3f>& force)
{
  ArrayHandle<Real> fp;
  auto Vlength = _para.GetParameter<Real>(PARA_VLENGTH);
  auto rhor_spline = _para.GetFieldAsArrayHandle<Vec7f>(field::rhor_spline);
  auto frho_spline = _para.GetFieldAsArrayHandle<Vec7f>(field::frho_spline);
  auto z2r_spline = _para.GetFieldAsArrayHandle<Vec7f>(field::z2r_spline);
  auto N = _position.GetNumberOfValues();

  auto rc = _para.GetParameter<Real>(PARA_CUTOFF);
  auto rs = _para.GetParameter<Real>(PARA_R_CORE);
  auto rho_system = _para.GetParameter<Real>(PARA_RHO);
  vtkm::cont::ArrayHandle<Vec3f> corr_force;
  corr_force.Allocate(N);

  vtkm::cont::ArrayHandle<vtkm::Id> num_verletlist;
  vtkm::cont::ArrayHandle<vtkm::Id> id_verletlist;
  vtkm::cont::ArrayHandle<vtkm::Vec3f> offset_verletlist;

  Id rs_num = rho_system * vtkm::Ceil(4.0 / 3.0 * vtkm::Pif() * rs * rs * rs) + 1;
  Id rc_num = rho_system * vtkm::Ceil(4.0 / 3.0 * vtkm::Pif() * rc * rc * rc) + 1;
  Id rcs_num = rc_num - rs_num;
  Real random_num = _para.GetParameter<Real>(PARA_NEIGHBOR_SAMPLE_NUM);
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

  RunWorklet::ComputeRBLNeighboursOnce(rs_num,
                                          pice_num,
                                          _atoms_id,
                                          _locator,
                                          id_verletlist_group,
                                          num_verletlist_group,
                                          offset_verletlist_group);

  RunWorklet::EAMfp(rc,
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

 RunWorklet::EAMRBLForce(rc,
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
 RunWorklet::SumRBLCorrForce(corr_value, corr_force, force);
}

void ExecutionMD::ComputeVerletlistEAMForce(ArrayHandle<Vec3f>& force)
{
 ArrayHandle<Real> fp;
 auto Vlength = _para.GetParameter<Real>(PARA_VLENGTH);
 auto rhor_spline = _para.GetFieldAsArrayHandle<Vec7f>(field::rhor_spline);
 auto frho_spline = _para.GetFieldAsArrayHandle<Vec7f>(field::frho_spline);
 auto z2r_spline = _para.GetFieldAsArrayHandle<Vec7f>(field::z2r_spline);

  auto cut_off = _para.GetParameter<Real>(PARA_CUTOFF);
  auto N = _position.GetNumberOfValues();
  auto rho_system = _para.GetParameter<Real>(PARA_RHO);
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

  RunWorklet::ComputeNeighbours(
    cut_off, _atoms_id, _locator, id_verletlist_group, num_verletlist, offset_verletlist_group);

  RunWorklet::EAMfpVerlet(cut_off,
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

  RunWorklet::EAMForceVerlet(cut_off,
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

void ExecutionMD::ComputeOriginalEAMForce(ArrayHandle<Vec3f>& force)
{
  auto cut_off = _para.GetParameter<Real>(EAM_PARA_CUTOFF);
  auto Vlength = _para.GetParameter<Real>(PARA_VLENGTH);

  auto rhor_spline = _para.GetFieldAsArrayHandle<Vec7f>(field::rhor_spline);
  auto frho_spline = _para.GetFieldAsArrayHandle<Vec7f>(field::frho_spline);
  auto z2r_spline = _para.GetFieldAsArrayHandle<Vec7f>(field::z2r_spline);

  ArrayHandle<Real> EAM_rho;
  ArrayHandle<Real> fp;

  //1:compute _EAM_rho   = density at each atom
  RunWorklet::EAM_rho(
    cut_off, Vlength, _atoms_id, rhor_spline, _locator, _topology, _force_function, EAM_rho);

  // 2:compute fp    = derivative of embedding energy at each atom
  RunWorklet::EAM_fp(_atoms_id, EAM_rho, frho_spline, _locator, _topology, _force_function, fp);

  // 3:compute force  = EAM_force
  RunWorklet::EAM_force(cut_off,
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

void ExecutionMD::ComputeSpecialBondsLJForce(ArrayHandle<Vec3f>& ljforce) 
{
  auto cut_off = _para.GetParameter<Real>(PARA_CUTOFF);

  auto special_offsets = _para.GetFieldAsArrayHandle<Id>(field::special_offsets);
  auto special_weights = _para.GetFieldAsArrayHandle<Real>(field::special_weights);
  auto specoal_ids = _para.GetFieldAsArrayHandle<Id>(field::special_ids);
  auto ids_group = vtkm::cont::make_ArrayHandleGroupVecVariable(specoal_ids, special_offsets);
  auto weight_group = vtkm::cont::make_ArrayHandleGroupVecVariable(special_weights, special_offsets);

  RunWorklet::SpecicalBondsLJForce(
    cut_off, _atoms_id, _locator, _topology, _force_function, ids_group, weight_group, ljforce);
}

void ExecutionMD::InitPointLocator()
{
  vtkm::Vec<vtkm::Range, 3> range = _para.GetParameter<vtkm::Vec<vtkm::Range, 3>>(PARA_RANGE);
  vtkm::Vec3f left_bottom{ {
                             static_cast<vtkm::FloatDefault>(range[0].Min),
                           },
                           {
                             static_cast<vtkm::FloatDefault>(range[1].Min),
                           },
                           {
                             static_cast<vtkm::FloatDefault>(range[2].Min),
                           } };
  vtkm::Vec3f right_top{ {
                           static_cast<vtkm::FloatDefault>(range[0].Max),
                         },
                         {
                           static_cast<vtkm::FloatDefault>(range[1].Max),
                         },
                         {
                           static_cast<vtkm::FloatDefault>(range[2].Max),
                         } };
  _locator.SetRange(left_bottom, right_top);

  _locator.SetCutOff(_para.GetParameter<Real>(PARA_CUTOFF));

  _locator.SetRs(_para.GetParameter<Real>(PARA_R_CORE));

  _locator.SetPosition(_position);

}