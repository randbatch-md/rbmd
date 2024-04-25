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

    //
  void ComputeVirial();
  void Compute_Pressure_Scalar();
  void Compute_Temp_Scalar();
  void Couple();
  void fix_press_berendsen();
  void x2lamda(Id n);
  void lamda2x(Id n);
  void set_global_box();
  void remap();

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

    //
  ArrayHandle<Vec6f> _virial_atom;
  Vec6f virial;          // accumulated virial: xx,yy,zz,xy,xz,yz
  Real _pressure_scalar; // computed global pressure scalar

  //
  Vec3f p_start, p_stop;
  Vec3f p_period, p_target;

  Vec3f p_current, dilation;
  Real bulkmodulus;


  //
  Vec6f h, h_inv; // shape matrix in Voigt ordering

  // orthogonal box

  Real xprd, yprd, zprd;                // global box dimensions
  Real xprd_half, yprd_half, zprd_half; // half dimensions
  Vec3f prd;                            // array form of dimensions
  Vec3f prd_half;                       // array form of half dimensions
};