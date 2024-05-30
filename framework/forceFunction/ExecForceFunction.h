#pragma once
#include <vtkm/VectorAnalysis.h>
#include <vtkm/cont/ArrayHandle.h>
#include <vtkm/cont/CoordinateSystem.h>
#include <vtkm/Deprecated.h>
#include <vtkm/Pair.h>
#include <vtkm/Math.h>
#include "Types.h"
#include "math/Math.h"

class ExecForceFunction
{
public:
  using IdPortalTypeId = typename vtkm::cont::ArrayHandle<vtkm::Id>::ReadPortalType;
  using IdPortalTypeReal = typename vtkm::cont::ArrayHandle<Real>::ReadPortalType;
  using IdPortalTypeVec2 = typename vtkm::cont::ArrayHandle<vtkm::Vec2i>::ReadPortalType;

  ExecForceFunction(const Real cut_off,
                    const Real& alpha,
                    const Real& volume,
                    const Real& vlength,
                    const IdComponent& kmax,
                    const IdComponent& rbeP)
    : _cut_Off(cut_off)
    , _alpha(alpha)
    , _volume(volume)
    , _Vlength(vlength)
    , Kmax(kmax)
    , RBEP(rbeP)
  { 
      _sum_gauss = Compute_S();
  }
 
  VTKM_EXEC vtkm::Vec3f ComputeLJForce(const Vec3f& r_ij,
                                       const Real& eps_i,
                                       const Real& eps_j,
                                       const Real& sigma_i,
                                       const Real& sigma_j,
                                       const Real& cut_off) const
  {
    const Real small_value = 0.0001;
    Vec3f force{ 0, 0, 0 };
    const Real dis_2 = r_ij[0] * r_ij[0] + r_ij[1] * r_ij[1] + r_ij[2] * r_ij[2];
    const Real cut_off_2 = cut_off * cut_off;

    if (dis_2 < cut_off_2 && dis_2 > small_value)
    {
      Real sigma_ij = (sigma_i + sigma_j) / 2;

      Real sigmaij_6 = sigma_ij * sigma_ij * sigma_ij * sigma_ij * sigma_ij * sigma_ij;
      Real dis_6 = dis_2 * dis_2 * dis_2;
      Real sigmaij_dis_6 = sigmaij_6 / dis_6;
      force = -24 * vtkm::Sqrt(eps_i * eps_j) * ((2 * sigmaij_dis_6 - 1) * sigmaij_dis_6) / dis_2 * r_ij;
    }
    return force;
  }

  VTKM_EXEC vtkm::Vec3f ComputeLJForceRcs(const Vec3f& r_ij,
                                          const Real& eps_i,
                                          const Real& eps_j,
                                          const Real& sigma_i,
                                          const Real& sigma_j,
                                          const Real& cut_off,
                                          const Real& rs) const
  {
    Vec3f force{ 0, 0, 0 };
    const Real dis_2 = r_ij[0] * r_ij[0] + r_ij[1] * r_ij[1] + r_ij[2] * r_ij[2];

    const Real cut_off_2 = cut_off * cut_off;
    const Real rs_2 = rs * rs;

    if (dis_2 < cut_off_2 && dis_2 > rs_2)
    {
      Real sigma_ij = (sigma_i + sigma_j) / 2;

      Real sigmaij_6 = sigma_ij * sigma_ij * sigma_ij * sigma_ij * sigma_ij * sigma_ij;
      Real dis_6 = dis_2 * dis_2 * dis_2;
      Real sigmaij_dis_6 = sigmaij_6 / dis_6;

      force = -24 * vtkm::Sqrt(eps_i * eps_j) * ((2 * sigmaij_dis_6 - 1) * sigmaij_dis_6) / dis_2 * r_ij;
    }
    return force;
  }

  VTKM_EXEC Vec3f ComputeNearEnergyForce(const Vec3f& p_i,
                                         const Vec3f& p_j,
                                         const Real& charge_pi,
                                         const Real& charge_pj) const
  {
    const Real small_value = 0.01;
    Vec3f force{ 0, 0, 0 };
    auto r_ij = p_j - p_i;
    Real dis = vtkm::Magnitude(r_ij);

    if (dis < _cut_Off && dis > small_value)
    {
      force = -charge_pi * charge_pj * Gnear(_alpha, dis) * r_ij / dis;
    }
    return force;
  }

  VTKM_EXEC Vec3f ComputeNearEnergyForceRcs(const Vec3f& p_i,
                                         const Vec3f& p_j,
                                         const Real& charge_pi,
                                         const Real& charge_pj,
                                         const Real& rs) const
  {
    const Real small_value = 0.01;
    Vec3f force{ 0, 0, 0 };
    auto r_ij = p_j - p_i;
    Real dis = vtkm::Magnitude(r_ij);

    if (dis < _cut_Off && dis > rs)
    {
      force = -charge_pi * charge_pj * Gnear(_alpha, dis) * r_ij / dis;
    }
    return force;
  }

  VTKM_EXEC Vec3f ComputeNearEnergyForceERF(const Vec3f& p_i,
                                         const Vec3f& p_j,
                                         const Real& charge_pi,
                                         const Real& charge_pj,
                                         const Real& table_pij) const
  {
    const Real small_value = 0.0001;
    Vec3f force{ 0, 0, 0 };
    auto r_ij = p_j - p_i;
    Real dis = vtkm::Magnitude(r_ij);
    Real dis_2 = vtkm::MagnitudeSquared(r_ij);
    auto rc_2 = _cut_Off * _cut_Off;
    if (dis_2 < rc_2 && dis_2 > small_value)
    {
      force = -charge_pi * charge_pj * table_pij * r_ij / dis;
    }
    return force;
  }

  VTKM_EXEC Vec3f ComputeNearEnergyForceRcsERF(const Vec3f& p_i,
                                          const Vec3f& p_j,
                                          const Real& charge_pi,
                                          const Real& charge_pj,
                                          const Real& table_pij,
                                          const Real& rs) const
  {
    const Real small_value = 0.0001;
    Vec3f force{ 0, 0, 0 };
    auto r_ij = p_j - p_i;
    Real dis = vtkm::Magnitude(r_ij);
    Real dis_2 = vtkm::MagnitudeSquared(r_ij);
    auto rc_2 = _cut_Off * _cut_Off;
    auto rs_2 = rs * rs;
    if (dis_2 < rc_2 && dis_2 > rs_2)
    {
      force = -charge_pi * charge_pj * table_pij * r_ij / dis;
    }
    return force;
  }

  VTKM_EXEC Real Gnear(const Real& _nearalpha, const Real& dis) const
  {
    Real erfcx = vtkm::Sqrt(_nearalpha) * dis;
    Real expx = -_nearalpha * dis * dis;
    Real Gnearvalue = (1.0 - vtkm::ERF(erfcx)) / (dis * dis) +
      2 * vtkm::Sqrt(_nearalpha) * vtkm::Exp(expx) / (vtkm::Sqrt(vtkm::Pif()) * dis);
    return Gnearvalue;
  }

  VTKM_EXEC Vec3f ComputeFarEle(const Vec3f& M,
                                const Vec3f& r_i,
                                const Real& charge_p_i,
                                const Vec2f& rhok_ri)
  {
    Vec3f force = { 0, 0, 0 };
  
    Vec3f K = 2 * vtkm::Pi() * M / _Vlength; // TODO: Lx
    Real range_K_2 = K[0] * K[0] + K[1] * K[1] + K[2] * K[2];
    auto factor_a = -4 * vtkm::Pi() * charge_p_i * K;
    auto factor_b = vtkm::Exp(-range_K_2 / (4 * _alpha));
    auto factor_c = vtkm::Cos(vtkm::Dot(K, r_i)) * rhok_ri[1];
    auto factor_d = vtkm::Sin(vtkm::Dot(K, r_i)) * rhok_ri[0];

    force = factor_a / (_volume * range_K_2) * factor_b * (factor_c - factor_d);  
    return force;
  }

  VTKM_EXEC Vec3f ComputeRBEForceSum(Vec3f& kl,
                                     const Vec3f& current_pts,
                                     const Real& current_charge,
                                     const Vec2f& rhok_ri)
  {
    Vec3f force = { 0, 0, 0 };
    kl = 2.0 * vtkm::Pi() * kl / _Vlength;
    Real range_kl_2 = kl[0] * kl[0] + kl[1] * kl[1] + kl[2] * kl[2];

    auto factor_a = -4 * vtkm::Pi() * current_charge * kl;
    auto factor_b = vtkm::Cos(vtkm::Dot(kl, current_pts)) * rhok_ri[1];
    auto factor_c = vtkm::Sin(vtkm::Dot(kl, current_pts)) * rhok_ri[0];

    force = (factor_a / (_volume * range_kl_2)) * (factor_b - factor_c);
    return force;
  }

  VTKM_EXEC Vec3f ComputeRBEForce(const Id& p, Vec3f& force)
  {
    return force * _sum_gauss / p;
  }

  VTKM_EXEC Real Compute_S() const
  {
    const Real& factor = Compute_H();
    Real factor_3 = factor * factor * factor;
    Real S = factor_3 - 1;
    return S;
  }

  VTKM_EXEC Real Compute_H() const
  {
    Real H = 0.0;
    const Real factor = -(_alpha * _Vlength * _Vlength);
    for (int m = -10; m <= 10; m++)
    {
      Real expx = m * m * factor;
      H += vtkm::Exp(expx);
    }
    H = H * vtkm::Sqrt(-(factor) / vtkm::Pi());
    return H;
  }

  VTKM_EXEC Real ComputePotentialEn(const Vec3f& r_ij,
                                    const Real& eps_ij,
                                    const Real& sigma_ij,
                                    const Real& cut_off) const
  {
    const Real small_value = 0.0001;
    Real cut_off_2 = cut_off * cut_off;
    Real ComptePE_ij = 0;
    Real dis_2 = r_ij[0] * r_ij[0] + r_ij[1] * r_ij[1] + r_ij[2] * r_ij[2];

    if (dis_2 < cut_off_2 && dis_2 > small_value)
    {
      Real sigmaij_6 = sigma_ij * sigma_ij * sigma_ij * sigma_ij * sigma_ij * sigma_ij;
      Real dis_6 = dis_2 * dis_2 * dis_2;
      ComptePE_ij = 4 * eps_ij * (sigmaij_6 / dis_6 - 1) * (sigmaij_6 / dis_6);
    }
    return 0.5 * ComptePE_ij;
  }

  VTKM_EXEC Real ComputeNearEleEnergy(const Vec3f& r_ij,
                                    const Real& charge_i,
                                    const Real& charge_j,
                                    const Real& cut_off,
                                    const Real& nearalpha) const
  {
    const Real small_value = 0.01;
    Real ComputePE_ij = 0;
    Real dis = vtkm::Magnitude(r_ij);
    if (dis < cut_off && dis > small_value)
    {
      ComputePE_ij = charge_i * charge_j * (1.0 - vtkm::ERF(vtkm::Sqrt(nearalpha) * dis)) / dis;
    }
    return 0.5 * ComputePE_ij;
  }

  private:
  Real _cut_Off;
  Real _alpha;
  Real _volume;
  Real _Vlength;
  Real _sum_gauss;
  IdComponent Kmax;
  IdComponent RBEP;

};
