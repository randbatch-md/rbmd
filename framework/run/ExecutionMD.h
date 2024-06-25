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
  std::vector<Vec2f> ComputeChargeStructureFactorEwald(Real& _Vlength, IdComponent& Kmax);
  std::vector<Vec2f> ComputeChargeStructureFactorRBE(Real& _Vlength, ArrayHandle<Vec3f>& _psample);
  void ComputeRBEEleForce(ArrayHandle<Vec3f>& psample,
                          IdComponent& RBE_P,
                          ArrayHandle<Vec3f>& RBE_ele_force);
  void ComputeEwaldEleForce(IdComponent& Kmax, ArrayHandle<Vec3f>& Ewald_ele_force);

  void ComputeRBLNearForce(ArrayHandle<Vec3f>& nearforce);
  void ComputeRBLLJForce(ArrayHandle<Vec3f>& LJforce);
  void ComputeVerletlistNearForce(ArrayHandle<Vec3f>& nearforce);
  void ComputeVerletlistLJForce(ArrayHandle<Vec3f>& ljforce);
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

  void  ComputeNewChargeStructureFactorRBE(Real& _Vlength,
                                                        ArrayHandle<Vec3f>& _psample,
                                                        ArrayHandle<Vec2f>& new_rhok);
  void InitPointLocator();
  void set_global_box();

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

};
