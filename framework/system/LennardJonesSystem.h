#pragma once
#include "MDSystem.h"
#include "Executioner.h"
#include <vtkm/cont/ArrayHandleGroupVecVariable.h>
#include "ContPointLocator.h"

class LennardJonesSystem : public MDSystem
{

public:
  LennardJonesSystem(const Configuration& cfg);
  virtual ~LennardJonesSystem(){};

  void Init() override;

private:
  void PreSolve() override;
  void Solve() override;
  void PostSolve() override;

  void InitField() override; 
  void TimeIntegration();
  void InitialCondition();
  void SetCenterTargetPositions();
  void ComputeForce();
  void UpdateVelocity();
  void UpdatePosition();
  void SetCharge() override;
  void ComputeTempe();
  void SetForceFunction() override;
  void SetTopology() override;  
  void ESS_RBL();

private:
  std::string _nearforce_type;
  vtkm::cont::Invoker _invoker;
  ArrayHandle<Vec3f> _force;
  Real _rho;
  Real _tempT_sum;
  Real _tempT;
  Real _dt;
  Real _nosehooverxi;
  vtkm::Float64 _rbl_pe_ave;

  //vtkm::cont::ArrayHandle<vtkm::Id> _num_verletlist;
  //vtkm::cont::ArrayHandle<vtkm::Id> _id_verletlist;
  //vtkm::cont::ArrayHandle<vtkm::Vec3f> _offset_verletlist;
};