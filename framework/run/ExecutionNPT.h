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
#include "ExecutionMD.h"
#include "Executioner.h"
#include <vtkm/cont/ArrayHandleGroupVecVariable.h>

class ExecutionNPT : public ExecutionMD
{
public:
    ExecutionNPT(const Configuration& cfg);
    virtual ~ExecutionNPT() {};

    void Init() override;

private:
    void PreSolve() override;
    void Solve() override;
    void PostSolve() override;
    void virtual TimeIntegration();
    vtkm::cont::ArrayHandle<Vec3f> LJForce();
    vtkm::cont::ArrayHandle<Vec3f> EleNearForce();
    vtkm::cont::ArrayHandle<Vec3f> NearForce();
    vtkm::cont::ArrayHandle<Vec3f> NearForceLJ();
    vtkm::cont::ArrayHandle<Vec3f> NearForceEAM();


    vtkm::cont::ArrayHandle<Vec3f> BondForce();
    vtkm::cont::ArrayHandle<Vec3f> AngleForce();
    vtkm::cont::ArrayHandle<Vec3f> DihedralsForce();
    vtkm::cont::ArrayHandle<Vec3f> SpecialCoulForce();
    vtkm::cont::ArrayHandle<Vec3f> EleNewForce();
    vtkm::cont::ArrayHandle<Vec6f> LJVirial();
    vtkm::cont::ArrayHandle<Vec6f> EwaldVirial();
    void TempConTypeForce();
    void ComputeTempe();
    void InitialCondition() override;
    void SetForceFunction() override;
    void SetTopology() override;
    void InitParameters() override;
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
    void ComputeVirial_r();
    void Compute_Pressure_Scalar();

    void Compute_Temp_Scalar();
    void Couple();
    void fix_press_berendsen();
    void fix_press_berendsen_scale();
    void x2lamda(Id n);
    void lamda2x(Id n);
    void set_global_box();
    void remap();

    //
    void Computedof();
    void SetUp();
    void NoseHooverChain();
    void ComputeTempTarget();
    void ComputePressTarget();
    void InitialIntegrate();
    void FinalIntegrate();

    void NHCTempIntegrate();
    void NHCPressIntegrate();
    void NVE_v();
    void NVE_x();
    void ComputeScalar();

    //
    void NHOmegaDot();
    void NH_V_Press();
    void ReMap();

    void ComputePCOM();
    void ComputeVCOM();
    void FixMomentum();

    //Shake
    void ShakeSetUp();
    void correct_coordinates();
    void correct_velocities();
    void PostForce();
    void Shake_Run();
private:
    ArrayHandle<Vec3f>   _com;
    ArrayHandle<Vec3f>   _vcom;
    ArrayHandle<Vec3f> _nearforce;
    ArrayHandle<Vec3f> _LJforce;
    ArrayHandle<Vec6f> _LJVirial;
    ArrayHandle<Vec3f> _ele_near_force;
    ArrayHandle<Vec3f> _ele_new_force;
    ArrayHandle<Vec3f> _all_force;
    ArrayHandle<Vec3f> _spec_coul_force;
    Real _dt;
    Real _kbT, _Tdamp;

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
    Vec3f _box;
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
    Id nrho, nr;
    std::vector<Real> frho;
    std::vector<Real> z2r;
    std::vector<Real> rhor;

    std::vector<Id> type2frho;
    std::vector<Id2> type2rhor;
    std::vector<Id2> type2z2r;
    std::vector<Vec2f> scale;
    Real dr, rdr, drho, rdrho, rhomax, rhomin;

    std::vector<Vec7f> frho_spline;
    std::vector<Vec7f> rhor_spline;
    std::vector<Vec7f> z2r_spline;
    std::vector<Real> rho;
    std::vector<Real> fp;
    std::vector<Id> numforce;

    std::vector<Id> map; // mapping from atom types to elements


    //
    Vec6f virial;          // accumulated virial: xx,yy,zz,xy,xz,yz
    Real _pressure_scalar; // computed global pressure scalar

    //pressure
    Real _Pstart, _Pstop, _Pperiod, _bulkmodulus; //read
    Real _Ptarget, _scale_factor, _pressure_coupling;

    Vec3f p_start, p_stop, p_period, p_target;
    Vec3f p_freq;
    Real p_freq_max;
    Real p_hydro; // hydrostatic target pressure

    Id6 p_flag;
    Id pdim; // number of barostatted dims

    Vec3f p_current, dilation;

    Vec6f omega, omega_dot;
    Vec6f omega_mass;
    Real mtk_term1, mtk_term2; // Martyna-Tobias-Klein corrections
    Id mtk_flag;              // 0 if using Hoover barostat

    //temp
    Real t_start, t_stop, t_period, t_target, ke_target;
    Real t_freq;
    Real tdof;

    std::vector<Real> eta, eta_dot; // chain thermostat for particles
    std::vector<Real> eta_dotdot;
    std::vector<Real> eta_mass;
    Id mtchain;              // length of chain
    Id mtchain_default_flag; // 1 = mtchain is default

    std::vector<Real> etap; // chain thermostat for barostat
    std::vector<Real> etap_dot;
    std::vector<Real> etap_dotdot;
    std::vector<Real> etap_mass;
    Id mpchain; // length of chain

    Id nc_tchain, nc_pchain;
    Real factor_eta;

    Real dtv, dtf, dthalf, dt4, dt8, dto;
    Id eta_mass_flag;
    Real drag, tdrag_factor; // drag factor on particle thermostat
    Real pdrag_factor;       // drag factor on barostat

    //
    Vec6f h, h_inv; // shape matrix in Voigt ordering

    // orthogonal box

    Real xprd, yprd, zprd;                // global box dimensions
    Real xprd_half, yprd_half, zprd_half; // half dimensions
    Vec3f prd;                            // array form of dimensions
    Vec3f prd_half;                       // array form of half dimensions

    ArrayHandle<Vec3f> _omega;
    ArrayHandle<Vec3f> _unwarp_position;
    Vec3f _com_all;
    Vec3f _vcom_all;
    Real _mass_all;

    ArrayHandle<Vec3f>  _force_tmp;
    ArrayHandle<Vec3f>  _velocity_tmp;

    ArrayHandle<Vec3f> _shake_force;
    ArrayHandle<Vec3f> _shake_position;
    ArrayHandle<Vec3f>  _velocity_p;
};