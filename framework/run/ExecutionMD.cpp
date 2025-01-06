//==================================================================================
//  RBMD 2.1.0 is developed for random batch molecular dynamics calculation.
//
//  Copyright(C) 2024 SHANGHAI JIAOTONG UNIVERSITY CHONGQING RESEARCH INSTITUTE
//
//  This program is free software : you can redistribute it and /or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with this program.If not, see < https://www.gnu.org/licenses/>.
//
//  The post-processing data produced by VASPKIT may be used in any publications 
//  provided that its use is explicitly acknowledged. A suitable reference for VASPKIT is:
//  [1] Gao W, Zhao T, Guo Y, et al.RBMD: A molecular dynamics package enabling to simulate 
//  10 million all - atom particles in a single graphics processing unit[J].arXiv preprint arXiv : 2407.09315, 2024.
// 
//  Contact Email : [support_wz@sciai.com.cn]
//==================================================================================

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
  _para.SetParameter(PARA_TEMPERATUREE_VECTOR, GetVectorValue<Real>("temperature"));
  _para.SetParameter(PARA_PRESSURE_VECTOR, GetVectorValue<Real>("pressure"));
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

void ExecutionMD::UpdateVerletList() 
{
  auto cut_off = _para.GetParameter<Real>(PARA_CUTOFF);
  auto box = _para.GetParameter<Vec3f>(PARA_BOX);
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
    cut_off, box,_atoms_id, _locator, id_verletlist_group, num_verletlist, offset_verletlist_group);

  _locator.SetVerletListInfo(num_verletlist, id_verletlist_group, offset_verletlist_group);
}

void ExecutionMD::ComputeCorrForce(vtkm::cont::ArrayHandle<Vec3f>& corr_force, vtkm::cont::ArrayHandle<Vec6f>& corr_virial)
{
  auto rc = _para.GetParameter<Real>(PARA_CUTOFF);
  auto rs = _para.GetParameter<Real>(PARA_R_CORE);
  auto rho_system = _para.GetParameter<Real>(PARA_RHO);
  auto N = _position.GetNumberOfValues();
  auto box = _para.GetParameter<Vec3f>(PARA_BOX);

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
                              box,
                              _atoms_id,
                              _locator,
                              _topology,
                              _force_function,
                              _static_table,
                              id_verletlist_group,
                              num_verletlist_group,
                              offset_verletlist_group,
                              corr_force, corr_virial);
}

std::vector<Vec2f> ExecutionMD::ComputeChargeStructureFactorRBE(Real& _Vlength, ArrayHandle<Vec3f>& _psample)
{
  std::vector<Vec2f> rhok;
  auto p_number = _psample.GetNumberOfValues();
  ArrayHandle<Real> density_real;
  ArrayHandle<Real> density_image;

  for (Id i = 0; i < p_number; i++)
  {
    auto kl = _psample.ReadPortal().Get(i);
    kl = 2 * vtkm::Pi() * kl / _Vlength;
    RunWorklet::ComputeChargeStructureFactorComponent(
      kl, _position, _charge, density_real, density_image);
    Real value_Re =
      vtkm::cont::Algorithm::Reduce(density_real, vtkm::TypeTraits<Real>::ZeroInitialization());
    Real value_Im =
      vtkm::cont::Algorithm::Reduce(density_image, vtkm::TypeTraits<Real>::ZeroInitialization());
    rhok.push_back({ value_Re, value_Im });
  }
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

std::vector<Vec2f> ExecutionMD::ComputeChargeStructureFactorEwald(Vec3f& box,
                                                                  Id3& Kmax,
                                                                  Real&  alpha,
                                                                  Real& ewald_energy_total)
{
  std::vector<Vec2f> rhok;
  ArrayHandle<Real> density_real;
  ArrayHandle<Real> density_image;
  for (Id i = -Kmax[0]; i <= Kmax[0]; i++)
  {
    for (Id j = -Kmax[1]; j <= Kmax[1]; j++)
    {
      for (Id k = -Kmax[2]; k <= Kmax[2]; k++)
      {
        if (!(i == 0 && j == 0 && k == 0))
        {
          Vec3f K = { Real(i), Real(j), Real(k) };
          //
          K = { Real(2 * vtkm::Pi() * K[0] / box[0]),
                Real(2 * vtkm::Pi() * K[1] / box[1]),
                Real(2 * vtkm::Pi() * K[2] / box[2]) };
          Real Range_K = vtkm::Magnitude(K);

          RunWorklet::ComputeChargeStructureFactorComponent(
            K, _position, _charge, density_real, density_image);
          Real value_Re = vtkm::cont::Algorithm::Reduce(
            density_real, vtkm::TypeTraits<Real>::ZeroInitialization());
          Real value_Im = vtkm::cont::Algorithm::Reduce(
            density_image, vtkm::TypeTraits<Real>::ZeroInitialization());

          Real Range_density2 = vtkm::Pow(value_Re, 2) + vtkm::Pow(value_Im, 2);

          ewald_energy_total +=
              vtkm::Exp(-Range_K * Range_K / (4 * alpha)) * Range_density2 / (Range_K * Range_K);

          rhok.push_back({ value_Re, value_Im });
        }
      }
    }
  }
  return rhok;
}

void  ExecutionMD::ComputeEwaldEnergy(Vec3f& box,
    Id3& Kmax,
    Real& alpha,
    Real& ewald_energy_total)
{
    std::vector<Vec2f> rhok;
    ArrayHandle<Real> density_real;
    ArrayHandle<Real> density_image;
    for (Id i = -Kmax[0]; i <= Kmax[0]; i++)
    {
        for (Id j = -Kmax[1]; j <= Kmax[1]; j++)
        {
            for (Id k = -Kmax[2]; k <= Kmax[2]; k++)
            {
                if (!(i == 0 && j == 0 && k == 0))
                {
                    Vec3f K = { Real(i), Real(j), Real(k) };
                    //
                    K = { Real(2 * vtkm::Pi() * K[0] / box[0]),
                          Real(2 * vtkm::Pi() * K[1] / box[1]),
                          Real(2 * vtkm::Pi() * K[2] / box[2]) };
                    Real Range_K = vtkm::Magnitude(K);

                    RunWorklet::ComputeChargeStructureFactorComponent(
                        K, _position, _charge, density_real, density_image);
                    Real value_Re = vtkm::cont::Algorithm::Reduce(
                        density_real, vtkm::TypeTraits<Real>::ZeroInitialization());
                    Real value_Im = vtkm::cont::Algorithm::Reduce(
                        density_image, vtkm::TypeTraits<Real>::ZeroInitialization());

                    Real Range_density2 = vtkm::Pow(value_Re, 2) + vtkm::Pow(value_Im, 2);

                    ewald_energy_total +=
                        vtkm::Exp(-Range_K * Range_K / (4 * alpha)) * Range_density2 / (Range_K * Range_K);
                }
            }
        }
    }
}


void ExecutionMD::ComputeRBEEleForce(ArrayHandle<Vec3f>& psample,
                                  IdComponent& RBE_P,
                                  ArrayHandle<Vec3f>& RBE_ele_force,
                                  ArrayHandle<Vec6f>& ewald_long_virial_atom)
{
  auto box = _para.GetParameter<Vec3f>(PARA_BOX);

  ArrayHandle<Vec2f> new_whole_rhok;
  ComputeNewChargeStructureFactorRBE(box, psample, new_whole_rhok);
  
  RunWorklet::ComputeNewRBEForce(RBE_P, box,
                                    _atoms_id,
                                    psample,
                                    new_whole_rhok,
                                    _force_function,
                                    _topology,
                                    _locator,
                                    RBE_ele_force,
                                    ewald_long_virial_atom);
 

}

void ExecutionMD::ComputeNewChargeStructureFactorRBE(Vec3f& _box,
                                                                ArrayHandle<Vec3f>& _psample,
                                                                ArrayHandle<Vec2f>& new_rhok)
{
 
  auto N = _position.GetNumberOfValues();
  
  const Id p_number = _psample.GetNumberOfValues();
  


  RunWorklet::ComputePnumberChargeStructureFactor(_box,
                                                     p_number,
                                                     N,
                                                     _atoms_id,
                                                     _position,
                                                     _charge,
                                                     _psample,
                                                     _rhok_Re,
                                                     _rhok_Im);
 

  vtkm::cont::ArrayHandle<Id> psamplekey_out;
  vtkm::cont::ArrayHandle<Real> rhok_Re_reduce;
  vtkm::cont::ArrayHandle<Real> rhok_Im_reduce;

   
  vtkm::cont::Algorithm::ReduceByKey<Id, Real>(_psamplekey, _rhok_Re, psamplekey_out, rhok_Re_reduce, vtkm::Add());
  vtkm::cont::Algorithm::ReduceByKey<Id, Real>(_psamplekey, _rhok_Im, psamplekey_out, rhok_Im_reduce, vtkm::Add());
 


  RunWorklet::ChangePnumberChargeStructureFactor(rhok_Re_reduce, rhok_Im_reduce, new_rhok);
 

  //energy
  Real self_potential_energy_ave;
  ComputeSelfEnergy(self_potential_energy_ave);

  auto Kmax_vec = _para.GetParameter<std::vector<Id>>(PARA_KMAX);
  Id3 kmax = { Kmax_vec[0],Kmax_vec[1], Kmax_vec[2] };
  auto alpha = _para.GetParameter<Real>(PARA_ALPHA);

  auto box = _para.GetParameter<Vec3f>(PARA_BOX);
  auto volume = box[0] * box[1] * box[2];
  Real ewald_energy_total = 0.0;
  ComputeEwaldEnergy(box, kmax, alpha, ewald_energy_total);
  Real ewald_energy_ave = _unit_factor._qqr2e * ewald_energy_total * ((2 * vtkm::Pi() / volume)) / N;
  ewald_energy_ave = ewald_energy_ave + self_potential_energy_ave;

  _para.SetParameter(PARA_EWALD_LONG_ENERGY, ewald_energy_ave);
}

void ExecutionMD::ComputeEwaldEleForce(Id3& Kmax, ArrayHandle<Vec3f>& Ewald_ele_force, ArrayHandle<Vec6f>& ewald_long_virial_atom)
{ 
  auto N = _position.GetNumberOfValues();
  auto alpha = _para.GetParameter<Real>(PARA_ALPHA);
  auto Vlength = _para.GetParameter<Real>(PARA_VLENGTH);
  auto box = _para.GetParameter<Vec3f>(PARA_BOX);
  auto volume = box[0] * box[1] * box[2];
  Real ewald_energy_total =0.0;
  auto rhok = ComputeChargeStructureFactorEwald(box, Kmax, alpha, ewald_energy_total);
  Real ewald_energy_ave = _unit_factor._qqr2e * ewald_energy_total * ((2 * vtkm::Pi() / volume)) /N;

  ArrayHandle<Vec2f> whole_rhok = vtkm::cont::make_ArrayHandle(rhok);
  RunWorklet::ComputeNewFarElectrostatics(
    Kmax, box,_atoms_id, whole_rhok, _force_function, _topology, _locator, Ewald_ele_force, ewald_long_virial_atom);

  Real  self_potential_energy_ave;
  ComputeSelfEnergy(self_potential_energy_ave);

  ewald_energy_ave = ewald_energy_ave + self_potential_energy_ave;
  _para.SetParameter(PARA_EWALD_LONG_ENERGY, ewald_energy_ave);
} 

void ExecutionMD::ComputeSelfEnergy(Real& self_potential_energy_ave)
{
    auto N = _position.GetNumberOfValues();
    auto alpha = _para.GetParameter<Real>(PARA_ALPHA);
    auto charge = _para.GetFieldAsArrayHandle<Real>(field::charge);

    ArrayHandle<Real> _self_energy;
    RunWorklet::ComputeSqCharge(charge, _self_energy);
    auto self_potential_energy_total = -vtkm::Sqrt(alpha / vtkm::Pi()) *
        vtkm::cont::Algorithm::Reduce(_self_energy, vtkm::TypeTraits<Real>::ZeroInitialization());
    self_potential_energy_ave = self_potential_energy_total / N;
    self_potential_energy_ave = self_potential_energy_ave * _unit_factor._qqr2e;

    _para.SetParameter(PARA_SELF_ENERGY, self_potential_energy_ave);
}

void ExecutionMD::ComputeRBLNearForce(ArrayHandle<Vec3f>& nearforce,
     ArrayHandle<Vec6f>& lj_coul_rbl_virial_atom)
{
  auto N = _position.GetNumberOfValues();
  vtkm::cont::ArrayHandle<Vec3f> corr_force;
  corr_force.Allocate(N);

  vtkm::cont::ArrayHandle<Vec6f> corr_virial;
  corr_virial.Allocate(N);

  auto rc = _para.GetParameter<Real>(PARA_CUTOFF);
  auto rs = _para.GetParameter<Real>(PARA_R_CORE);
  auto box = _para.GetParameter<Vec3f>(PARA_BOX);
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

  if (_para.GetParameter<bool>(PARA_DIHEDRALS_FORCE))
  {
    RunWorklet::NearForceRBLERFSpecialBonds(rs_num,
                                            pice_num,
                                            _unit_factor._qqr2e,
                                            box,
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
                                            corr_force, corr_virial);
  }
  else
  {
    RunWorklet::NearForceRBLERF(rs_num,
                                pice_num,
                                _unit_factor._qqr2e,
                                box,
                                _atoms_id,
                                _locator,
                                _topology,
                                _force_function,
                                _static_table,
                                id_verletlist_group,
                                num_verletlist_group,
                                offset_verletlist_group,
                                corr_force, corr_virial);
  }

  Vec3f corr_value = vtkm::cont::Algorithm::Reduce(corr_force, vtkm::TypeTraits<Vec3f>::ZeroInitialization()) / N;
  RunWorklet::SumRBLCorrForce(corr_value, corr_force, nearforce);

  Vec6f corr_virial_value = vtkm::cont::Algorithm::Reduce(corr_virial, vtkm::TypeTraits<Vec6f>::ZeroInitialization()) / N;
  RunWorklet::SumRBLCorrVirial(corr_virial_value, corr_virial, lj_coul_rbl_virial_atom);

  //energy
  ComputeLJCoulEnergy();

}

void ExecutionMD::ComputeLJCoulEnergy() 
{
    auto cut_off = _para.GetParameter<Real>(PARA_CUTOFF);
    auto box = _para.GetParameter<Vec3f>(PARA_BOX);
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

    RunWorklet::ComputeNeighbours(cut_off, box, _atoms_id, _locator, id_verletlist_group,
        num_verletlist, offset_verletlist_group);

    //energy
    ArrayHandle<Real> energy_lj;
    ArrayHandle<Real>energy_coul;
    if (_para.GetParameter<bool>(PARA_DIHEDRALS_FORCE)) 
    {
        auto special_offsets = _para.GetFieldAsArrayHandle<Id>(field::special_offsets);
        auto special_weights = _para.GetFieldAsArrayHandle<Real>(field::special_weights);
        auto specoal_ids = _para.GetFieldAsArrayHandle<Id>(field::special_ids);
        auto ids_group = vtkm::cont::make_ArrayHandleGroupVecVariable(specoal_ids, special_offsets);
        auto weight_group =
            vtkm::cont::make_ArrayHandleGroupVecVariable(special_weights, special_offsets);

        RunWorklet::SpecialLJCoulEnergyVerletWeightVirial(cut_off,
            box,
            _unit_factor._qqr2e,
            _atoms_id,
            _locator,
            _topology,
            _force_function,
            id_verletlist_group,
            num_verletlist,
            offset_verletlist_group,
            ids_group,
            weight_group,
            energy_lj, energy_coul);
    }
    else {
        RunWorklet::LJCoulEnergyVerletVirial(cut_off,
            box,
            _unit_factor._qqr2e,
            _atoms_id,
            _locator,
            _topology,
            _force_function,
            id_verletlist_group,
            num_verletlist,
            offset_verletlist_group,
            energy_lj, energy_coul);
    }

    auto lj_potential_energy_total = vtkm::cont::Algorithm::Reduce(
        energy_lj, vtkm::TypeTraits<Real>::ZeroInitialization());
    auto lj_potential_energy_avr = lj_potential_energy_total / N;

    auto coul_potential_energy_total = vtkm::cont::Algorithm::Reduce(
        energy_coul, vtkm::TypeTraits<Real>::ZeroInitialization());
    auto coul_potential_energy_avr = coul_potential_energy_total / N;

    _para.SetParameter(PARA_LJ_ENERGY, lj_potential_energy_avr);
    _para.SetParameter(PARA_COUL_ENERGY, coul_potential_energy_avr);
}

void ExecutionMD::ComputeRBLLJForce(ArrayHandle<Vec3f>& LJforce, ArrayHandle<Vec6f>& lj_rbl_virial_atom)
{
  auto N = _position.GetNumberOfValues();
  vtkm::cont::ArrayHandle<Vec3f> corr_ljforce;
  corr_ljforce.Allocate(N);
  vtkm::cont::ArrayHandle<Vec6f> corr_ljvirial;
  corr_ljvirial.Allocate(N);

  auto box = _para.GetParameter<Vec3f>(PARA_BOX);

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
                         box,
                         _atoms_id,
                         _locator,
                         _topology,
                         _force_function,
                         id_verletlist_group,
                         num_verletlist_group,
                         offset_verletlist_group,
                         corr_ljforce, corr_ljvirial);

  Vec3f corr_value =
    vtkm::cont::Algorithm::Reduce(corr_ljforce, vtkm::TypeTraits<Vec3f>::ZeroInitialization()) / N;
  RunWorklet::SumRBLCorrForce(corr_value, corr_ljforce, LJforce);

  Vec6f corr_virial_value =
      vtkm::cont::Algorithm::Reduce(corr_ljvirial, vtkm::TypeTraits<Vec6f>::ZeroInitialization()) / N;
  RunWorklet::SumRBLCorrVirial(corr_virial_value, corr_ljvirial, lj_rbl_virial_atom);

  //energy
  ComputeLJEnergy();
}

void ExecutionMD::ComputeLJEnergy()
{
    auto rc = _para.GetParameter<Real>(PARA_CUTOFF);
    auto box = _para.GetParameter<Vec3f>(PARA_BOX);
    auto N = _position.GetNumberOfValues();
    auto rho_system = _para.GetParameter<Real>(PARA_RHO);
    auto max_j_num = rho_system * vtkm::Ceil(4.0 / 3.0 * vtkm::Pif() * rc * rc * rc) + 1;
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

    RunWorklet::ComputeNeighbours(rc, box, _atoms_id, _locator, id_verletlist_group,
        num_verletlist, offset_verletlist_group);


    //
    ArrayHandle<Real> energy_lj;
    RunWorklet::LJEnergyVerlet(rc,
        box,
        _atoms_id,
        _locator,
        _topology,
        _force_function,
        id_verletlist_group,
        num_verletlist,
        offset_verletlist_group,
        energy_lj);

    auto lj_potential_energy_total = vtkm::cont::Algorithm::Reduce(
        energy_lj, vtkm::TypeTraits<Real>::ZeroInitialization());
    auto lj_potential_energy_avr = lj_potential_energy_total / N;
    _para.SetParameter(PARA_LJ_ENERGY, lj_potential_energy_avr);
}

void ExecutionMD::ComputeVerletlistNearForce(ArrayHandle<Vec3f>& nearforce,ArrayHandle<Vec6f>& lj_coul_virial_atom)
{
  auto cut_off = _para.GetParameter<Real>(PARA_CUTOFF);
  auto box = _para.GetParameter<Vec3f>(PARA_BOX);
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

  RunWorklet::ComputeNeighbours( cut_off, box,_atoms_id, _locator, id_verletlist_group,
      num_verletlist, offset_verletlist_group);

  ArrayHandle<Real> energy_lj, energy_coul;
  if (_para.GetParameter<bool>(PARA_DIHEDRALS_FORCE))
  {
      //
      auto special_offsets = _para.GetFieldAsArrayHandle<Id>(field::special_offsets);
      auto special_weights = _para.GetFieldAsArrayHandle<Real>(field::special_weights);
      auto specoal_ids = _para.GetFieldAsArrayHandle<Id>(field::special_ids);
      auto ids_group = vtkm::cont::make_ArrayHandleGroupVecVariable(specoal_ids, special_offsets);
      auto weight_group = vtkm::cont::make_ArrayHandleGroupVecVariable(special_weights, special_offsets);

      RunWorklet::NearForceVerletWeightVirial(cut_off,
          box,
          _unit_factor._qqr2e,
          _atoms_id,
          _locator,
          _topology,
          _force_function,
          id_verletlist_group,
          num_verletlist,
          offset_verletlist_group,
          ids_group,
          weight_group,
          nearforce,
          lj_coul_virial_atom,
          energy_lj, energy_coul);
  }
  else
  {
      RunWorklet::NearForceVerletVirial(cut_off,
          box,
          _unit_factor._qqr2e,
          _atoms_id,
          _locator,
          _topology,
          _force_function,
          id_verletlist_group,
          num_verletlist,
          offset_verletlist_group,
          nearforce,
          lj_coul_virial_atom,
          energy_lj, energy_coul);
  }

  auto lj_potential_energy_total = vtkm::cont::Algorithm::Reduce(
      energy_lj, vtkm::TypeTraits<Real>::ZeroInitialization());
  auto lj_potential_energy_avr = lj_potential_energy_total / N;

  auto coul_potential_energy_total = vtkm::cont::Algorithm::Reduce(
      energy_coul, vtkm::TypeTraits<Real>::ZeroInitialization());
  auto coul_potential_energy_avr = coul_potential_energy_total / N;
  
  _para.SetParameter(PARA_LJ_ENERGY, lj_potential_energy_avr);
  _para.SetParameter(PARA_COUL_ENERGY, coul_potential_energy_avr);
}

void ExecutionMD::ComputeVerletlistLJForce(ArrayHandle<Vec3f>& ljforce, ArrayHandle<Vec6f>& lj_virial_atom)
{
  auto cut_off = _para.GetParameter<Real>(PARA_CUTOFF);
  auto box = _para.GetParameter<Vec3f>(PARA_BOX);
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

  RunWorklet::ComputeNeighbours( cut_off, box, _atoms_id, _locator, id_verletlist_group, 
      num_verletlist, offset_verletlist_group);

  ArrayHandle<Real> energy_lj;
  RunWorklet::LJForceVerlet(cut_off,
                            box,
                            _atoms_id,
                            _locator,
                            _topology,
                            _force_function,
                            id_verletlist_group,
                            num_verletlist,
                            offset_verletlist_group,
                            ljforce,
                            lj_virial_atom,
                            energy_lj);
  auto lj_potential_energy_total = vtkm::cont::Algorithm::Reduce(
      energy_lj, vtkm::TypeTraits<Real>::ZeroInitialization());
  auto lj_potential_energy_avr = lj_potential_energy_total / N;

  _para.SetParameter(PARA_LJ_ENERGY, lj_potential_energy_avr);
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
  auto box = _para.GetParameter<Vec3f>(PARA_BOX);
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
                    box,
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
                          box,
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
 auto box = _para.GetParameter<Vec3f>(PARA_BOX);
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
    cut_off, box,_atoms_id, _locator, id_verletlist_group, num_verletlist, offset_verletlist_group);

  RunWorklet::EAMfpVerlet(cut_off,
                          box,
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
                             box,
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
  auto box = _para.GetParameter<Vec3f>(PARA_BOX);

  auto rhor_spline = _para.GetFieldAsArrayHandle<Vec7f>(field::rhor_spline);
  auto frho_spline = _para.GetFieldAsArrayHandle<Vec7f>(field::frho_spline);
  auto z2r_spline = _para.GetFieldAsArrayHandle<Vec7f>(field::z2r_spline);

  ArrayHandle<Real> EAM_rho;
  ArrayHandle<Real> fp;

  //1:compute _EAM_rho   = density at each atom
  RunWorklet::EAM_rho(
    cut_off, box, _atoms_id, rhor_spline, _locator, _topology, _force_function, EAM_rho);

  // 2:compute fp    = derivative of embedding energy at each atom
  RunWorklet::EAM_fp(_atoms_id, EAM_rho, frho_spline, _locator, _topology, _force_function, fp);

  // 3:compute force  = EAM_force
  RunWorklet::EAM_force(cut_off,
                        box,
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

void ExecutionMD::ApplyPbc()
{
    //pbc
    auto&& position_flag = _para.GetFieldAsArrayHandle<Id3>(field::position_flag);
    auto box = _para.GetParameter<Vec3f>(PARA_BOX); //
    auto range = _para.GetParameter<vtkm::Vec<vtkm::Range, 3>>(PARA_RANGE);
    Vec<Vec2f, 3> data_range{ { static_cast<Real>(range[0].Min), static_cast<Real>(range[0].Max) },
                              { static_cast<Real>(range[1].Min), static_cast<Real>(range[1].Max) },
                              { static_cast<Real>(range[2].Min), static_cast<Real>(range[2].Max) } };
    RunWorklet::ApplyPbcFlag(box, data_range, _position, _locator, position_flag);
}

void ExecutionMD::Computedof()
{
    auto n = _position.GetNumberOfValues();
    auto extra_dof = 3; //dimension =3
    tdof = 3 * n - extra_dof;
    auto shake = _para.GetParameter<std::string>(PARA_FIX_SHAKE);
    if (shake == "true")
    {
        tdof = tdof - n;
    }
}

void ExecutionMD::ComputeVirial()
{
    Vec6f lj_virial = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
    Vec6f lj_coul_virial = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

    Vec6f spec_lj_virial, ewald_long_virial, spec_coul_virial, spec_ljcoul_virial, bond_virial, angle_virial, dihedral_virial;
    Vec6f shake_virial;

    spec_lj_virial = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
    spec_ljcoul_virial = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
    ewald_long_virial = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
    spec_coul_virial = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
    bond_virial = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
    angle_virial = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
    dihedral_virial = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
    shake_virial = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

    auto force_field = _para.GetParameter<std::string>(PARA_FORCE_FIELD_TYPE);
    auto nearforce_type = _para.GetParameter<std::string>(PARA_NEIGHBOR_TYPE);
    if ("LJ/CUT" == force_field)
    {
        if ("RBL" == nearforce_type)
        {
            lj_virial = vtkm::cont::Algorithm::Reduce(_lj_rbl_virial_atom,
                vtkm::TypeTraits<Vec6f>::ZeroInitialization());
        }
        else if ("VERLETLIST" == nearforce_type)
        {
            lj_virial = vtkm::cont::Algorithm::Reduce(_lj_virial_atom,
                vtkm::TypeTraits<Vec6f>::ZeroInitialization());
        }
        virial = lj_virial;
    }
    //
    if ("LJ/CUT/COUL/LONG" == force_field)
    {
        if ("RBL" == nearforce_type)
        {
            lj_coul_virial = vtkm::cont::Algorithm::Reduce(_lj_coul_rbl_virial_atom,
                vtkm::TypeTraits<Vec6f>::ZeroInitialization());
        }
        else if ("VERLETLIST" == nearforce_type)
        {
            lj_coul_virial = vtkm::cont::Algorithm::Reduce(_lj_coul_virial_atom,
                vtkm::TypeTraits<Vec6f>::ZeroInitialization());
        }
        ewald_long_virial = vtkm::cont::Algorithm::Reduce(
            _ewald_long_virial_atom, vtkm::TypeTraits<Vec6f>::ZeroInitialization()) *
            _unit_factor._qqr2e;

        virial = lj_coul_virial + ewald_long_virial;
    }

    if ("CVFF" == force_field)
    {
        if ("RBL" == nearforce_type)
        {
            spec_lj_virial = vtkm::cont::Algorithm::Reduce(_lj_coul_rbl_virial_atom,
                vtkm::TypeTraits<Vec6f>::ZeroInitialization());

            spec_coul_virial = vtkm::cont::Algorithm::Reduce(_spec_coul_virial_atom,
                vtkm::TypeTraits<Vec6f>::ZeroInitialization());
        }
        else if ("VERLETLIST" == nearforce_type)
        {
            spec_ljcoul_virial = vtkm::cont::Algorithm::Reduce(_lj_coul_virial_atom,
                vtkm::TypeTraits<Vec6f>::ZeroInitialization());
        }
        //ComputeEwaldLongVirial(_Kmax, _ewald_long_virial_atom);
        ewald_long_virial = vtkm::cont::Algorithm::Reduce(
            _ewald_long_virial_atom, vtkm::TypeTraits<Vec6f>::ZeroInitialization()) *
            _unit_factor._qqr2e;

        bond_virial = vtkm::cont::Algorithm::Reduce(_bond_virial_atom,
            vtkm::TypeTraits<Vec6f>::ZeroInitialization());

        angle_virial = vtkm::cont::Algorithm::Reduce(_angle_virial_atom,
            vtkm::TypeTraits<Vec6f>::ZeroInitialization());

        if (_para.GetParameter<bool>(PARA_DIHEDRALS_FORCE) &&
            _para.GetParameter<bool>(PARA_FILE_DIHEDRALS))
        {
            dihedral_virial = vtkm::cont::Algorithm::Reduce(
                _dihedral_virial_atom, vtkm::TypeTraits<Vec6f>::ZeroInitialization());
        }
        if (_para.GetParameter<std::string>(PARA_FIX_SHAKE) == "true")
        {
            shake_virial = vtkm::cont::Algorithm::Reduce(
                _shake_first_virial_atom, vtkm::TypeTraits<Vec6f>::ZeroInitialization());
        }

        //total virial
        if ("RBL" == nearforce_type)
        {
            virial = spec_lj_virial + spec_coul_virial + ewald_long_virial + bond_virial + angle_virial +
                dihedral_virial + shake_virial;
        }
        else if ("VERLETLIST" == nearforce_type)
        {
            virial = spec_ljcoul_virial + ewald_long_virial + bond_virial + angle_virial +
                dihedral_virial + shake_virial;
        }
    }
}
