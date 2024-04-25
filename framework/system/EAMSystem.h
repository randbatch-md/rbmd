#pragma once
#include "MDSystem.h"
#include "Executioner.h"
#include <vtkm/cont/ArrayHandleGroupVecVariable.h>

class EAMSystem : public MDSystem
{
public:
  EAMSystem(const Configuration& cfg);
  virtual ~EAMSystem(){};

  void Init() override;

private:
  void PreSolve() override;
  void Solve() override;
  void PostSolve() override;

  void InitField() override; 
  void UpdateVelocityByTempConType();
  void TimeIntegration() ;
  void InitialCondition() override;
  void SetCenterTargetPositions();
  void ComputeForce();
  void UpdateVelocity();
  void UpdatePosition();
  void SetForceFunction() override;
  void SetTopology() override;
  void SetCharge() override;
  void ComputeTempe();
  void PreForce();

  void ReadPotentialFile(std::ifstream& file);
  void AllocateEAM();
  void file2array();
  void interpolate(Id n, Real delta, std::vector<Real>& f, std::vector<Vec7f>& spline);
  void array2spline();
  void SetEAM();
  void InitStyle();

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
  std::shared_ptr<Executioner>& _executioner;
  ArrayHandle<Vec3f> _force;
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
  std::ifstream _potential_file;


  //
  ArrayHandle<Vec6f> _virial_atom; 
  Vec6f virial;                      // accumulated virial: xx,yy,zz,xy,xz,yz    
  Real _pressure_scalar;              // computed global pressure scalar

  //
  Vec3f p_start, p_stop;
  Vec3f p_period ,p_target;

  Vec3f p_current, dilation;
  Real bulkmodulus;


  //
  Vec6f h, h_inv; // shape matrix in Voigt ordering

   // orthogonal box

  Real  xprd, yprd, zprd;                  // global box dimensions
  Real  xprd_half, yprd_half, zprd_half; // half dimensions
  Vec3f prd;                           // array form of dimensions
  Vec3f prd_half;                      // array form of half dimensions

  
  ArrayHandle<Real> _EAM_rho;
  ArrayHandle<Real> _fp;
  
  struct Funcfl
  {
    Id nrho, nr;
    Real drho, dr, cut_off;
    std::vector<Real> frho;
    std::vector<Real> zr;
    std::vector<Real> rhor;
  };
  Funcfl file;

  // potentials as array data
  Id nrho, nr;
  std::vector<Real> frho;
  std::vector<Real> z2r;
  std::vector<Real> rhor;

  std::vector<Id> type2frho;
  std::vector<Id2> type2rhor;
  std::vector<Id2> type2z2r;
  std::vector<Vec2f> scale;

  // potentials in spline form used for force computation
  Real dr, rdr, drho, rdrho, rhomax, rhomin;
  
  std::vector<Vec7f> frho_spline;
  std::vector<Vec7f> rhor_spline;
  std::vector<Vec7f> z2r_spline;
  //  per-atom arrays
  std::vector<Real> rho;
  std::vector<Real> fp;
  std::vector<Id> numforce;

  std::vector<Id> map; // mapping from atom types to elements
};