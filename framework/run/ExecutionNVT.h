#pragma once
#include "Executioner.h"
#include "ExecutionMD.h"
#include <vtkm/cont/ArrayHandleGroupVecVariable.h>

class ExecutionNVT : public ExecutionMD
{
public:
  ExecutionNVT(const Configuration& cfg);
  virtual ~ExecutionNVT(){};

  void Init() override;

private:
  void PreSolve() override;
  void Solve() override;
  void PostSolve() override;
  void virtual TimeIntegration();
  vtkm::cont::ArrayHandle<Vec3f> LJForce();
  vtkm::cont::ArrayHandle<Vec3f> EleNearForce();
  vtkm::cont::ArrayHandle<Vec3f> NearForce(); // NearForce: EleNearForce and LJforce
  vtkm::cont::ArrayHandle<Vec3f> NearForceLJ(); // NearForce: EleNearForce and LJforce
  vtkm::cont::ArrayHandle<Vec3f> NearForceEAM(); // NearForce: EleNearForce and LJforce


  vtkm::cont::ArrayHandle<Vec3f> BondForce();
  vtkm::cont::ArrayHandle<Vec3f> AngleForce();
  vtkm::cont::ArrayHandle<Vec3f> DihedralsForce();
  vtkm::cont::ArrayHandle<Vec3f> SpecialCoulForce();
  vtkm::cont::ArrayHandle<Vec3f> EleNewForce();
  void TempConTypeForce();
  void ComputeTempe();
  void InitialCondition() override;
  void SetForceFunction() override;
  void SetTopology() override;
  void InitParameters() override;
  //void InitField() override;
  void ComputeForce();
  void ComputeAllForce();
  void UpdateVelocity();
  void UpdatePosition();
  void UpdateVelocityByTempConType();
  void SetCenterTargetPositions();
  void PreForce();

  void ConstraintA();
  void ConstraintB();

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
  ArrayHandle<Vec3f> _nearforce;
  ArrayHandle<Vec3f> _LJforce;
  ArrayHandle<Vec3f> _ele_near_force;
  ArrayHandle<Vec3f> _ele_new_force;
  ArrayHandle<Vec3f> _all_force;
  ArrayHandle<Vec3f> _spec_coul_force;
  Real _dt;
  Real _kbT ,_Tdamp;
  Real _Tstart, _Tstop, _Tperiod, _Ttarget;
  Real _Pstart, _Pstop, _Pperiod, _Ptarget;
  std::string _farforce_type;
  std::string _nearforce_type;
  std::string _temp_con_type;
  bool _use_shake;

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


  
  std::ifstream _potential_file;


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