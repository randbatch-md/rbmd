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

#pragma once
#include "Execution.h"
#include "locator/ContPointLocator.h"
#include "forceFunction/ContForceFunction.h"
#include "topology/ContTopology.h"
#include "UnitFactor.h"
#include <vtkm/cont/ArrayCopy.h>
#include "staticTable/ContStaticTable.h"
#include <vtkm/cont/ArrayHandle.h>
class ExecutionMD : public Execution
{
public:
  ExecutionMD(const Configuration& cfg);
  virtual ~ExecutionMD() = default;

  void Init() override;
  void Execute() override;
  void SetParameters() override;
  virtual void PreSolve(){};
  virtual void Solve(){};
  virtual void PostSolve(){};
  virtual void InitParameters();

protected:
  std::vector<Vec2f> ComputeChargeStructureFactorEwald(Vec3f& _box, IdComponent& Kmax, Real& alpha,Real& ewald_energy_total);
  std::vector<Vec2f> ComputeChargeStructureFactorRBE(Real& _Vlength, ArrayHandle<Vec3f>& _psample);
  void ComputeRBEEleForce(ArrayHandle<Vec3f>& psample,
                          IdComponent& RBE_P,
                          ArrayHandle<Vec3f>& RBE_ele_force,
                          ArrayHandle<Vec6f>& ewald_long_virial_atom);
  void ComputeEwaldEleForce(IdComponent& Kmax, ArrayHandle<Vec3f>& Ewald_ele_force, ArrayHandle<Vec6f>& ewald_long_virial_atom);
  void ComputeLJCoulEnergy(ArrayHandle<Vec6f>& nearVirial_atom);

  void ComputeLJEnergy(ArrayHandle<Vec6f>& nearVirial_atom);

  void ComputeSelfEnergy(Real& self_potential_energy_avr);
  void ComputeEwaldEnergy(Vec3f& _box, IdComponent& Kmax, Real& alpha, Real& ewald_energy_total);

  void ComputeRBLNearForce(ArrayHandle<Vec3f>& nearforce, ArrayHandle<Vec6f>& nearVirial_atom);
  void ComputeRBLLJForce(ArrayHandle<Vec3f>& LJforce, ArrayHandle<Vec6f>& nearVirial_atom);
  void ComputeVerletlistNearForce(ArrayHandle<Vec3f>& nearforce, ArrayHandle<Vec6f>& nearVirial_atom);

  void ComputeVerletlistLJForce(ArrayHandle<Vec3f>& ljforce, ArrayHandle<Vec6f>& nearVirial_atom);
  void ComputeOriginalLJForce(ArrayHandle<Vec3f>& ljforce);
  void ComputeSpecialBondsLJForce(ArrayHandle<Vec3f>& ljforce);
  void ComputeRBLEAMForce(ArrayHandle<Vec3f>& force);
  void ComputeVerletlistEAMForce(ArrayHandle<Vec3f>& force);
  void ComputeOriginalEAMForce(ArrayHandle<Vec3f>& force);


  void UpdateVerletList();

  void ComputeCorrForce(vtkm::cont::ArrayHandle<Vec3f>& corr_force);
  virtual void InitialCondition();
  virtual void SetForceFunction();
  virtual void SetTopology();
  virtual void SetCharge(){};
  void InitERF();

  void ComputeNewChargeStructureFactorRBE(Vec3f& _box,
                                          ArrayHandle<Vec3f>& _psample,
                                          ArrayHandle<Vec2f>& new_rhok);
  void InitPointLocator();


  //
  void ComputeVerletlistLJVirial(ArrayHandle<Vec6f>& lj_virial);
  void ComputeCoulVirial(ArrayHandle<Vec6f>& Coul_virial);
  void ComputeEwaldLongVirial(IdComponent& Kmax, ArrayHandle<Vec6f>& Ewald_long_virial);

  void ApplyPbc();
protected:
  ArrayHandle<Id> _molecule_id;
  ArrayHandle<Id> _atoms_id;
  ArrayHandle<Real> _charge;
  ArrayHandle<Real> _mass;
  ArrayHandle<Real> _rhok_Re;
  ArrayHandle<Real> _rhok_Im;
  ArrayHandle<Vec3f> _velocity;
  ArrayHandle<Id> _psamplekey;
  IdComponent _RBE_P;
  ContForceFunction _force_function;
  ContTopology _topology;
  bool _use_erf;

  vtkm::cont::Timer _EleFartimer;
  Real _Elefartimer_counting;
  vtkm::cont::Timer _EleNearPairtimer;
  Real _EleNearPairtimer_counting;

  std::string _unit;
  UnitFactor _unit_factor;
  ArrayHandle<Vec3f> _position;
  ContPointLocator _locator;

  ContStaticTable _static_table;
  std::string _init_way;

  ArrayHandle<Vec6f> _nearVirial_atom;
  ArrayHandle<Vec6f> _lj_virial_atom;
  ArrayHandle<Vec6f> _coul_virial_atom;
  ArrayHandle<Vec6f> _ewald_long_virial_atom;
  ArrayHandle<Vec6f> _ewald_long_virial_atom_1;
  ArrayHandle<Vec6f> _bond_virial_atom;
  ArrayHandle<Real> _bond_virial_atom_1;
  ArrayHandle<Vec6f> _angle_virial_atom;
  ArrayHandle<Real> _angle_virial_atom_1;
  ArrayHandle<Vec6f> _dihedral_virial_atom;
  ArrayHandle<Vec6f> _spec_coul_virial_atom;

  ArrayHandle<Vec6f> _shake_first_virial_atom;

  ArrayHandle<Real> _energy_lj, _energy_coul;
};
