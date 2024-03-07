#pragma once
#include "MDSystem.h"
#include "Executioner.h"
#include <vtkm/cont/ArrayHandleGroupVecVariable.h>

class NaClSystem : public MDSystem
{
public:
  NaClSystem(const Configuration& cfg);
  virtual ~NaClSystem(){};

  void Init() override;

private:
  void PreSolve() override;
  void Solve() override;
  void PostSolve() override;
  void virtual TimeIntegration();
  vtkm::cont::ArrayHandle<Vec3f> LJForce();
  vtkm::cont::ArrayHandle<Vec3f> EleNearForce();
  vtkm::cont::ArrayHandle<Vec3f> NearForce(); // NearForce: EleNearForce and LJforce
  vtkm::cont::ArrayHandle<Vec3f> BondForce();
  vtkm::cont::ArrayHandle<Vec3f> AngleForce();
  vtkm::cont::ArrayHandle<Vec3f> SpecialCoulForce();
  vtkm::cont::ArrayHandle<Vec3f> EleNewForce();
  void TempConTypeForce();
  void ComputeTempe();
  void InitialCondition() override;
  void SetForceFunction() override;
  void SetTopology() override;
  void InitField() override;
  void ComputeForce();
  void ComputeAllForce();
  void UpdateVelocity();
  void UpdatePosition();
  void UpdateVelocityByTempConType();
  void SetCenterTargetPositions();
  void PreForce();
  void Rattle();
  void ConstraintA();
  void ConstraintB();

private:
  ArrayHandle<Vec3f> _nearforce;
  ArrayHandle<Vec3f> _LJforce;
  ArrayHandle<Vec3f> _ele_near_force;
  ArrayHandle<Vec3f> _ele_new_force;
  ArrayHandle<Vec3f> _all_force;
  ArrayHandle<Vec3f> _spec_coul_force;
  Real _dt;
  Real _kbT;
  std::string _farforce_type;
  std::string _nearforce_type;
  std::string _temp_con_type;
  bool _use_shake;
  
  IdComponent _RBE_P;
  ArrayHandle<Vec3f> _psample;
  Vec3f _gaussian;

  std::shared_ptr<Executioner>& _executioner;

  ArrayHandle<Vec3f> _old_velocity;
  ArrayHandle<Vec3f> _old_position;

  IdComponent _Kmax;
  Real _cut_off;
  Real _nosehooverxi;
  Real _Vlength;
  Real _Volume;
  Real _alpha;
  Real _tempT_sum;
  Real _tempT;
};