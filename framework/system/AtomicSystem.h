#pragma once
#include "Executioner.h"
#include "MeshFreeSystem.h"
#include <vtkm/cont/ArrayHandleGroupVecVariable.h>

class AtomicSystem : public MeshFreeSystem
{
public:
  AtomicSystem(const Configuration& cfg);
  virtual ~AtomicSystem() = default;

  void Init() override;

protected:
  void InitField() override;

  virtual void InitialCondition(); 
  virtual void SetForceFunction(){};
  virtual void SetTopology(){};
  virtual void SetCharge(){};

  std::vector<Vec2f> ComputeChargeStructureFactorEwald(Real& _Vlength, IdComponent& Kmax);
  std::vector<Vec2f> ComputeChargeStructureFactorRBE(Real& _Vlength,ArrayHandle<Vec3f>& _psample);
  void ComputeRBEEleForce(ArrayHandle<Vec3f>& psample,
                          IdComponent& RBE_P,
                          ArrayHandle<Vec3f>& RBE_ele_force);
  void ComputeEwaldEleForce(IdComponent& Kmax,
                            ArrayHandle<Vec3f>& Ewald_ele_force);

protected:
  ContForceFunction _force_function;
  ContTopology _topology;
  ArrayHandle<Id> _atoms_id;
  ArrayHandle<Real> _charge;
  ArrayHandle<Real> _mass;
  ArrayHandle<Vec3f> _velocity;
};