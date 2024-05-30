#pragma once
#include "MDSystem.h"
#include "Executioner.h"
#include <vtkm/cont/ArrayHandleGroupVecVariable.h>

class SalineSolutionSystem : public MDSystem
{
public:
  SalineSolutionSystem(const Configuration& cfg);
  virtual ~SalineSolutionSystem(){};

  void Init() override;

private:
  void PreSolve() override;
  void Solve() override;
  void PostSolve() override;

  void InitField() override; 
  void UpdateVelocityByTempConType();
  void ComputeAllForce();
  void TimeIntegration() ;
  void InitialCondition() override;
  void SetCenterTargetPositions();
  void ComputeForce();
  void UpdateVelocity();
  void UpdatePosition();
  void SetForceFunction() override;
  void SetTopology() override;
  void PreForce();
  void SetCharge() override;
  vtkm::cont::ArrayHandle<Vec3f> LJForce();
  vtkm::cont::ArrayHandle<Vec3f> EleNearForce();
  vtkm::cont::ArrayHandle<Vec3f> EleNewForce();
  vtkm::cont::ArrayHandle<Vec3f> NearForce();   // NearForce: EleNearForce and LJforce
  void TempConTypeForce();
  void ComputeTempe();
 
private:
  std::shared_ptr<Executioner>& _executioner;
  ArrayHandle<Vec3f> _nearforce;
  ArrayHandle<Vec3f> _LJforce;
  ArrayHandle<Vec3f> _ele_near_force;
  ArrayHandle<Vec3f> _ele_new_force;
  ArrayHandle<Vec3f> _all_force;
  ArrayHandle<Vec3f> _psample;
  Vec3f _gaussian;
  std::string _farforce_type;
  std::string _nearforce_type;
  std::string _temp_con_type;
  IdComponent _RBE_P;
  IdComponent _Kmax;
  Real _dt;
  Real _kbT;
  Real _nosehooverxi;
  Real _Vlength;
  Real _alpha;
  Real _rho;
  Real _tempT_sum;
  Real _tempT;
};