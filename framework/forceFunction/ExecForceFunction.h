#pragma once
#include <vtkm/VectorAnalysis.h>
#include <vtkm/cont/ArrayHandle.h>
#include <vtkm/cont/CoordinateSystem.h>
#include <vtkm/Deprecated.h>
#include <vtkm/Pair.h>
#include <vtkm/Math.h>
#include "Types.h"
#include "math/Math.h"
#define MIN(A, B) ((A) < (B) ? (A) : (B))
#define MAX(A, B) ((A) > (B) ? (A) : (B))

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
                    const Vec3f& box,
                    const IdComponent& kmax,
                    const IdComponent& rbeP,
                    const Real& rhomax,
                    const Id& nrho,
                    const Real& drho,
                    const Id& nr,
                    const Real& dr)
    : _cut_Off(cut_off)
    , _alpha(alpha)
    , _volume(volume)
    , _Vlength(vlength)
    , _box(box)
    , Kmax(kmax)
    , RBEP(rbeP)
    , _rhomax(rhomax)
    , _nrho(nrho)
    , _drho(drho)
    , _nr(nr)
    , _dr(dr)
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

  VTKM_EXEC Vec3f ComputeNearEnergyForceERF_box(const Vec3f& r_ij,
                                                const Real& charge_pi,
                                                const Real& charge_pj,
                                                const Real& table_pij) const
  {
    const Real small_value = 0.0001;
    Vec3f force{ 0, 0, 0 };
    Real dis = vtkm::Magnitude(r_ij);
    Real dis_2 = vtkm::MagnitudeSquared(r_ij);
    auto rc_2 = _cut_Off * _cut_Off;
    if (dis_2 < rc_2 && dis_2 > small_value)
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
    Vec3f K{ 0, 0, 0 };
    K[0] = 2 * vtkm::Pi() * M[0] / _box[0];
    K[1] = 2 * vtkm::Pi() * M[1] / _box[1];
    K[2] = 2 * vtkm::Pi() * M[2] / _box[2];
    //Vec3f K = 2 * vtkm::Pi() * M / _Vlength; // TODO: Lx
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
    kl[0] = 2.0 * vtkm::Pi() * kl[0] / _box[0];
    kl[1] = 2.0 * vtkm::Pi() * kl[1] / _box[1];
    kl[2] = 2.0 * vtkm::Pi() * kl[2] / _box[2];
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
    const Vec3f& factor = Compute_H();
    Real factor_3 = factor[0] * factor[1] * factor[2];
    Real S = factor_3 - 1;
    return S;
  }

  VTKM_EXEC Vec3f Compute_H() const
  {
    Vec3f H = { 0, 0, 0 };
    for (Id i = 0; i < 3; ++i)
    {
      const Real factor = -(_alpha * _box[i] * _box[i]);
      for (int m = -10; m <= 10; m++)
      {
        Real expx = m * m * factor;
        H[i] += vtkm::Exp(expx);
      }
      H[i] *= vtkm::Sqrt(-(factor) / vtkm::Pi());
    }

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

    //compute EAMforce;
  template<typename rhor_splineType>
  VTKM_EXEC Real ComputeEAMrho(const Real& rc,
                               const Vec3f& r_ij,
                               const rhor_splineType& rhor_spline) const
  {
    Real EAM_rho = 0;
    auto rdr = 1 / _dr;

    auto rsq = r_ij[0] * r_ij[0] + r_ij[1] * r_ij[1] + r_ij[2] * r_ij[2];
    auto cutsq = rc * rc;

    if (rsq < cutsq && rsq > 0.01)
    {
      auto p = vtkm::Sqrt(rsq) * rdr + 1.0;
      auto m = static_cast<Id>(p);
      // Id m = p;
      m = MIN(m, _nr - 1);
      p -= m;
      p = MIN(p, 1.0);
      auto coeff = rhor_spline.Get(m);
      EAM_rho = ((coeff[3] * p + coeff[4]) * p + coeff[5]) * p + coeff[6];
    }
    return EAM_rho;
  }

  template<typename rhor_splineType>
  VTKM_EXEC Real ComputeEAMrhoRs(const Real& rc,
                                 const Real& rs,
                                 const Vec3f& r_ij,
                                 const rhor_splineType& rhor_spline) const
  {
    Real EAM_rho = 0;
    auto rdr = 1 / _dr;

    auto rsq = r_ij[0] * r_ij[0] + r_ij[1] * r_ij[1] + r_ij[2] * r_ij[2];
    auto cutsq = rc * rc;
    auto rs_2 = rs * rs;
    if (rsq < cutsq && rsq > rs_2)
    {
      auto p = vtkm::Sqrt(rsq) * rdr + 1.0;
      auto m = static_cast<Id>(p);
      // Id m = p;
      m = MIN(m, _nr - 1);
      p -= m;
      p = MIN(p, 1.0);
      auto coeff = rhor_spline.Get(m);
      EAM_rho = ((coeff[3] * p + coeff[4]) * p + coeff[5]) * p + coeff[6];
    }
    return EAM_rho;
  }

  template<typename rhoType, typename frho_splineType>
  VTKM_EXEC Real ComputeEAMfp(const Id& atoms_id,
                              const rhoType& EAM_rho,
                              const frho_splineType& frho_spline) const
  {
    Real fp = 0;
    auto rdrho = 1 / _drho;
    //auto p = EAM_rho.Get(atoms_id) * rdrho + 1.0;
    auto p = EAM_rho * rdrho + 1.0;
    auto m = static_cast<int>(p);
    //Id m = p;
    m = MAX(1, MIN(m, _nrho - 1));
    p -= m;
    p = MIN(p, 1.0);
    auto coeff = frho_spline.Get(m);
    fp = (coeff[0] * p + coeff[1]) * p + coeff[2];
    return fp;
  }

  template<typename rhoType, typename frho_splineType>
  VTKM_EXEC Real ComputeEAMfpOriginal(const Id& atoms_id,
                              const rhoType& EAM_rho,
                              const frho_splineType& frho_spline) const
  {
    Real fp = 0;
    auto rdrho = 1 / _drho;
    auto p = EAM_rho.Get(atoms_id) * rdrho + 1.0;
    //auto p = EAM_rho * rdrho + 1.0;
    auto m = static_cast<int>(p);
    //Id m = p;
    m = MAX(1, MIN(m, _nrho - 1));
    p -= m;
    p = MIN(p, 1.0);
    auto coeff = frho_spline.Get(m);
    fp = (coeff[0] * p + coeff[1]) * p + coeff[2];
    return fp;
  }


  template<typename fpType, typename rhor_splineType, typename z2r_splineType>
  VTKM_EXEC Vec3f ComputeEAMforce(const Real& eam_cut_off,
                                  const Id& atoms_id,
                                  const Id& pts_id_j,
                                  const Vec3f& r_ij,
                                  const fpType fp,
                                  const rhor_splineType& rhor_spline,
                                  const z2r_splineType& z2r_spline) const
  {
    Vec3f eam_force = { 0, 0, 0 };

    auto rdr = 1 / _dr;
    auto rsq = r_ij[0] * r_ij[0] + r_ij[1] * r_ij[1] + r_ij[2] * r_ij[2];
    auto cutsq = eam_cut_off * eam_cut_off;

    if (rsq < cutsq && rsq > 0.01)
    {
      auto r = vtkm::Sqrt(rsq);
      auto p = r * rdr + 1.0;
      auto m = static_cast<int>(p);
      //Id m = p;
      m = MIN(m, _nr - 1);
      p -= m;
      p = MIN(p, 1.0);

      // rhoip = derivative of (density at atom j due to atom i)
      // rhojp = derivative of (density at atom i due to atom j)
      // phi = pair potential energy
      // phip = phi'
      // z2 = phi * r
      // z2p = (phi * r)' = (phi' r) + phi
      // psip needs both fp[i] and fp[j] terms since r_ij appears in two
      //   terms of embed eng: Fi(sum rho_ij) and Fj(sum rho_ji)
      //   hence embed' = Fi(sum rho_ij) rhojp + Fj(sum rho_ji) rhoip
      // scale factor can be applied by thermodynamic integration

      auto coeffi = rhor_spline.Get(m);
      auto rhoip = (coeffi[0] * p + coeffi[1]) * p + coeffi[2];

      auto coeffj = rhor_spline.Get(m);
      auto rhojp = (coeffj[0] * p + coeffj[1]) * p + coeffj[2];

      auto coeff = z2r_spline.Get(m);
      auto z2p = (coeff[0] * p + coeff[1]) * p + coeff[2];
      auto z2 = ((coeff[3] * p + coeff[4]) * p + coeff[5]) * p + coeff[6];

      auto recip = 1.0 / r;
      auto phi = z2 * recip;                 //pair potential energy
      auto phip = z2p * recip - phi * recip; //pair force
      auto psip = fp.Get(atoms_id) * rhojp + fp.Get(pts_id_j) * rhoip + phip;
      auto fpair = -psip * recip;
      //std::cout << fpair << "," << r_ij[0] << "," << r_ij[1] << ","  << r_ij[2] << std::endl;

      //compute f
      eam_force[0] = r_ij[0] * fpair;
      eam_force[1] = r_ij[1] * fpair;
      eam_force[2] = r_ij[2] * fpair;
    }
    return eam_force;
  }

  template<typename fpType, typename rhor_splineType, typename z2r_splineType>
  VTKM_EXEC Vec3f ComputeEAMforceRBL(const Real& eam_cut_off,
                                     const Real& rs,
                                     const Id& atoms_id,
                                     const Id& pts_id_j,
                                     const Vec3f& r_ij,
                                     const fpType fp,
                                     const rhor_splineType& rhor_spline,
                                     const z2r_splineType& z2r_spline) const
  {
    Vec3f eam_force = { 0, 0, 0 };

    auto rdr = 1 / _dr;
    auto rsq = r_ij[0] * r_ij[0] + r_ij[1] * r_ij[1] + r_ij[2] * r_ij[2];
    auto cutsq = eam_cut_off * eam_cut_off;
    auto rs_2 = rs * rs;
    if (rsq < cutsq && rsq > rs_2)
    {
      auto r = vtkm::Sqrt(rsq);
      auto p = r * rdr + 1.0;
      auto m = static_cast<int>(p);
      //Id m = p;
      m = MIN(m, _nr - 1);
      p -= m;
      p = MIN(p, 1.0);

      auto coeffi = rhor_spline.Get(m);
      auto rhoip = (coeffi[0] * p + coeffi[1]) * p + coeffi[2];

      auto coeffj = rhor_spline.Get(m);
      auto rhojp = (coeffj[0] * p + coeffj[1]) * p + coeffj[2];

      auto coeff = z2r_spline.Get(m);
      auto z2p = (coeff[0] * p + coeff[1]) * p + coeff[2];
      auto z2 = ((coeff[3] * p + coeff[4]) * p + coeff[5]) * p + coeff[6];

      auto recip = 1.0 / r;
      auto phi = z2 * recip;                 //pair potential energy
      auto phip = z2p * recip - phi * recip; //pair force
      auto psip = fp.Get(atoms_id) * rhojp + fp.Get(pts_id_j) * rhoip + phip;
      auto fpair = -psip * recip;
      //std::cout << fpair << "," << r_ij[0] << "," << r_ij[1] << ","  << r_ij[2] << std::endl;

      //compute f
      eam_force[0] = r_ij[0] * fpair;
      eam_force[1] = r_ij[1] * fpair;
      eam_force[2] = r_ij[2] * fpair;
    }
    return eam_force;
  }


  //compute EAM Energy = EmbeddingEnergy + PairEnergy;
  template<typename rhor_splineType>
  VTKM_EXEC Real ComputeEAMrhoOUT(const Real& eam_cut_off,
                                  const Vec3f& r_ij,
                                  const rhor_splineType& rhor_spline) const
  {
    Real EAM_rho = 0;
    auto rdr = 1 / _dr;

    auto rsq = r_ij[0] * r_ij[0] + r_ij[1] * r_ij[1] + r_ij[2] * r_ij[2];
    auto cutsq = eam_cut_off * eam_cut_off;

    if (rsq < cutsq && rsq > 0.01)
    {
      auto p = vtkm::Sqrt(rsq) * rdr + 1.0;
      auto m = static_cast<Id>(p);
      //Id m = p;
      m = MIN(m, _nr - 1);
      p -= m;
      p = MIN(p, 1.0);
      auto coeff = rhor_spline.Get(m);
      EAM_rho = ((coeff[3] * p + coeff[4]) * p + coeff[5]) * p + coeff[6];
    }
    return EAM_rho;
  }

  template<typename rhoType, typename frho_splineType>
  VTKM_EXEC Real ComputeEmbeddingEnergy(const Id& atoms_id,
                                        const rhoType& EAM_rho,
                                        const frho_splineType& frho_spline) const
  {
    Real phi = 0;

    auto rdrho = 1 / _drho;
    auto p = EAM_rho.Get(atoms_id) * rdrho + 1.0;
    auto m = static_cast<int>(p);
    //Id m = p;
    m = MAX(1, MIN(m, _nrho - 1));
    p -= m;
    p = MIN(p, 1.0);
    auto coeff = frho_spline.Get(m);
    auto fp = (coeff[0] * p + coeff[1]) * p + coeff[2];
    phi = ((coeff[3] * p + coeff[4]) * p + coeff[5]) * p + coeff[6];
    if (EAM_rho.Get(atoms_id) > _rhomax)
    {
      phi += fp * (EAM_rho.Get(atoms_id) - _rhomax);
    }
    return phi;
  }

  template<typename z2r_splineType>
  VTKM_EXEC Real ComputePairEnergy(const Real& eam_cut_off,
                                   const Vec3f& r_ij,
                                   const z2r_splineType& z2r_spline) const
  {
    Real phi = 0;
    auto rdr = 1 / _dr;
    auto rsq = r_ij[0] * r_ij[0] + r_ij[1] * r_ij[1] + r_ij[2] * r_ij[2];
    auto cutsq = eam_cut_off * eam_cut_off;

    if (rsq < cutsq && rsq > 0.01)
    {
      auto r = vtkm::Sqrt(rsq);
      auto p = r * rdr + 1.0;
      auto m = static_cast<int>(p);
      //Id m = p;
      m = MIN(m, _nr - 1);
      p -= m;
      p = MIN(p, 1.0);

      auto coeff = z2r_spline.Get(m);
      auto z2p = (coeff[0] * p + coeff[1]) * p + coeff[2];
      auto z2 = ((coeff[3] * p + coeff[4]) * p + coeff[5]) * p + coeff[6];

      auto recip = 1.0 / r;
      phi = z2 * recip; //pair potential energy
    }
    return 0.5 * phi;
  }

  private:
  Real _cut_Off;
  Real _alpha;
  Real _volume;
  Real _Vlength;
  Vec3f _box;
  Real _sum_gauss;
  IdComponent Kmax;
  IdComponent RBEP;

  Real _rhomax;
  Id _nr;
  Id _nrho;
  Real _dr;
  Real _drho;

};
