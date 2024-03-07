#pragma once
#include "MeshFreeSystem.h"
class MolecularSystem : public MeshFreeSystem
{
public:
  MolecularSystem(const Configuration& cfg);
  virtual ~MolecularSystem() = default;

  void Init() override;

protected:
  std::vector<Vec2f> ComputeChargeStructureFactorEwald(Real& _Vlength, IdComponent& Kmax);
  std::vector<Vec2f> ComputeChargeStructureFactorRBE(Real& _Vlength,ArrayHandle<Vec3f>& _psample);
  void ComputeRBEEleForce(ArrayHandle<Vec3f>& psample,
                          IdComponent & RBE_P,
                          ArrayHandle<Vec3f>& RBE_ele_force);
  void ComputeEwaldEleForce(IdComponent& Kmax,
                            ArrayHandle<Vec3f>& Ewald_ele_force);

  void InitField() override;

  virtual void InitialCondition(); 
  virtual void SetForceFunction();
  virtual void SetTopology();

protected:
  ArrayHandle<Id> _molecule_id;
  ArrayHandle<Id> _atoms_id;
  ArrayHandle<Real> _charge;
  ArrayHandle<Real> _mass;
  ArrayHandle<Vec3f> _velocity;
  ContForceFunction _force_function;
  ContTopology _topology;
};