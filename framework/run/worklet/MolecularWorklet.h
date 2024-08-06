#pragma once

#include <vtkm/worklet/WorkletMapField.h>
#include <vtkm/worklet/WorkletMapTopology.h>
#include <vtkm/worklet/WorkletReduceByKey.h>
#include "locator/ExecPointLocator.h"
#include "forceFunction/ExecForceFunction.h"

static constexpr Real SMALL = 0.001;

namespace MolecularWorklet
{
struct UnitRescaleWorklet : vtkm::worklet::WorkletMapField
{
  UnitRescaleWorklet(const Real& unit_factor)
    : _unit_factor(unit_factor)
  {
  }
  using ControlSignature = void(FieldInOut unit_rescale);
  using ExecutionSignature = void(_1);

  template<typename UnitRescaleType>
  VTKM_EXEC void operator()(UnitRescaleType& unit_rescale) const
  {
    unit_rescale = unit_rescale * _unit_factor;
  }

  Real _unit_factor;
};

struct ComputeBondHarmonicWorklet : vtkm::worklet::WorkletMapField
{
  ComputeBondHarmonicWorklet(const Vec3f& box)
    : _box(box)
  {
  }
  using ControlSignature = void(FieldIn bond_type,
                                FieldIn bondlist,
                                FieldIn bond_coeffs_k,
                                FieldIn bond_coeffs_equilibrium,
                                WholeArrayIn _position,
                                FieldOut forcebond,
                                FieldOut bondEnergy,
                                ExecObject locator);
  using ExecutionSignature = void(_1, _2, _3, _4, _5, _6, _7, _8);

  template<typename BondType,
           typename BondListType,
           typename PositionType,
           typename ForceBondType,
           typename BondEnergyType>
  VTKM_EXEC void operator()(const BondType& bond_type,
                            const BondListType& bondlist,
                            const Real& bond_coeffs_k,
                            const Real& bond_coeffs_equilibrium,
                            const PositionType& _position,
                            ForceBondType& forcebond,
                            BondEnergyType& bondEnergy,
                            const ExecPointLocator& locator) const
  {

    vtkm::IdComponent bondi = bondlist[0];
    vtkm::IdComponent bondj = bondlist[1];
    vtkm::IdComponent bondtype = bond_type;
    Vec3f forcebondij;

    ComputeijBondEnergyForce(bondi,
                             bondj,
                             bondtype,
                             bond_coeffs_k,
                             bond_coeffs_equilibrium,
                             _position,
                             forcebondij,
                             bondEnergy,
                             locator);

    //_forcebond.Get(bondi) = _forcebond.Get(bondi) + forcebondij;
    //_forcebond.Get(bondj) = _forcebond.Get(bondj) - forcebondij;

    forcebond[0] = forcebondij;
    forcebond[1] = -forcebondij;
  }
  template<typename PositionType>
  VTKM_EXEC void ComputeijBondEnergyForce(const vtkm::IdComponent& bondi,
                                          const vtkm::IdComponent& bondj,
                                          const vtkm::IdComponent& bondtype,
                                          const Real& k,
                                          const Real& equilibrium,
                                          const PositionType& _position,
                                          Vec3f& forcebondij,
                                          Real& bondEnergy,
                                          const ExecPointLocator& locator) const
  {
    Vec3f p_i = _position.Get(bondi);
    Vec3f p_j = _position.Get(bondj);

    // minimum image distance
    Vec3f r_ij = locator.MinDistanceVec(p_i, p_j, _box);

    Real dis_2 = vtkm::MagnitudeSquared(r_ij);
    Real disij = vtkm::Sqrt(dis_2);

    Real dr = disij - equilibrium;
    Real rk = k * dr;

    if (disij > 0.01)
      forcebondij = -2.0 * rk / disij;
    else
      forcebondij = 0.0;
    forcebondij = r_ij * forcebondij;

    bondEnergy = rk * dr;
  }
  Vec3f _box;
};

struct ReduceForceWorklet : vtkm::worklet::WorkletReduceByKey
{
  using ControlSignature = void(KeysIn atoms, ValuesIn force, ReducedValuesOut reduce_force);

  using ExecutionSignature = _3(_2);
  using InputDomain = _1;
  template<typename ForceVecType>
  VTKM_EXEC typename ForceVecType::ComponentType operator()(const ForceVecType& force) const
  {
    typename ForceVecType::ComponentType sum = 0;
    for (vtkm::IdComponent index = 0; index < force.GetNumberOfComponents(); index++)
    {
      sum = sum + force[index];
    }
    return sum;
  }
};

struct AddForceWorklet : vtkm::worklet::WorkletMapField
{
  using ControlSignature = void(FieldIn force, FieldInOut allForce);
  using ExecutionSignature = void(_1, _2);

  template<typename ForceType, typename allForceType>
  VTKM_EXEC void operator()(const ForceType& force, allForceType& allForce) const
  {
    allForce = allForce + force;
  }
};

struct ComputeAngleHarmonicWorklet : vtkm::worklet::WorkletMapField
{

  ComputeAngleHarmonicWorklet(const Vec3f& box)
    : _box(box)
  {
  }
  using ControlSignature = void(FieldIn angle_type,
                                const FieldIn anglelist,
                                FieldIn angle_coeffs_k,
                                FieldIn angle_coeffs_equilibrium,
                                const WholeArrayIn whole_pts,
                                FieldOut forceangle,
                                FieldOut angleEnergy,
                                ExecObject locator);
  using ExecutionSignature = void(_1, _2, _3, _4, _5, _6, _7, _8);

  template<typename AngleType,
           typename AngleListType,
           typename WholePtsType,
           typename ForceAngleType,
           typename AngleEnergyType>
  VTKM_EXEC void operator()(const AngleType& angle_type,
                            const AngleListType& anglelist,
                            const Real& angle_coeffs_k,
                            const Real& angle_coeffs_equilibrium,
                            const WholePtsType& whole_pts,
                            ForceAngleType& forceangle,
                            AngleEnergyType& angleEnergy,
                            const ExecPointLocator& locator) const
  {

    vtkm::IdComponent anglei = anglelist[0];
    vtkm::IdComponent anglej = anglelist[1];
    vtkm::IdComponent anglek = anglelist[2];
    vtkm::IdComponent angletype = angle_type;
    Vec3f force_anglei;
    Vec3f force_anglek;

    ComputeijAngleEnergyForce(anglei,
                              anglej,
                              anglek,
                              angletype,
                              angle_coeffs_k,
                              angle_coeffs_equilibrium,
                              whole_pts,
                              force_anglei,
                              force_anglek,
                              angleEnergy,
                              locator);
    forceangle[0] = force_anglei;
    forceangle[1] = -force_anglei - force_anglek;
    forceangle[2] = force_anglek;
  }
  template<typename WholePtsType>
  VTKM_EXEC void ComputeijAngleEnergyForce(const vtkm::IdComponent& anglei,
                                           const vtkm::IdComponent& anglej,
                                           const vtkm::IdComponent& anglek,
                                           const vtkm::IdComponent& angletype,
                                           const Real& k,
                                           const Real& equilibrium,
                                           const WholePtsType& whole_pts,
                                           Vec3f& force_anglei,
                                           Vec3f& force_anglek,
                                           Real& angleEnergy,
                                           const ExecPointLocator& locator) const
  {

    Vec3f p_i = whole_pts.Get(anglei);
    Vec3f p_j = whole_pts.Get(anglej);
    Vec3f p_k = whole_pts.Get(anglek);

    // minimum image distance
    Vec3f r_ij = locator.MinDistanceVec(p_i, p_j, _box);
    Vec3f r_kj = locator.MinDistanceVec(p_k, p_j, _box);

    Real disij_2 = vtkm::MagnitudeSquared(r_ij);
    Real disij = vtkm::Sqrt(disij_2);

    Real diskj_2 = vtkm::MagnitudeSquared(r_kj);
    Real diskj = vtkm::Sqrt(diskj_2);

    Real cosangle = vtkm::Dot(r_ij, r_kj);
    cosangle /= disij * diskj;

    if (cosangle > 1.0)
      cosangle = 1.0;
    if (cosangle < -1.0)
      cosangle = -1.0;
    Real s = sqrt(1.0 - cosangle * cosangle);
    if (s < SMALL)
      s = SMALL;
    s = 1.0 / s;

    Real dtheta = acos(cosangle) - (equilibrium * vtkm::Pif()) / 180;
    Real tk = k * dtheta;
    angleEnergy = tk * dtheta;

    Real a = -2.0 * tk * s;
    Real a11 = a * cosangle / disij_2;
    Real a12 = -a / (disij * diskj);
    Real a22 = a * cosangle / diskj_2;

    force_anglei = a11 * r_ij + a12 * r_kj;
    force_anglek = a22 * r_kj + a12 * r_ij;
  }

   Vec3f _box;
};


struct ComputeDihedralOPLSWorklet : vtkm::worklet::WorkletMapField
{

  ComputeDihedralOPLSWorklet(const Vec3f& box)
     : _box(box)
  {
  }
  using ControlSignature = void(FieldIn dihedral_type,
                                const FieldIn dihedrallist,
                                FieldIn dihedral_coeffs_k,
                                FieldIn dihedral_coeffs_equilibrium,
                                const WholeArrayIn whole_pts,
                                FieldOut forcedihedral,
                                FieldOut dihedralEnergy);
  using ExecutionSignature = void(_1, _2, _3, _4, _5, _6, _7);

  template<typename DihedralType,
           typename DihedralListType,
           typename WholePtsType,
           typename ForceDihedralType,
           typename DihedralEnergyType>
  VTKM_EXEC void operator()(const DihedralType& dihedral_type,
                            const DihedralListType& dihedrallist,
                            const Real& dihedral_coeffs_k,
                            const Real& dihedral_coeffs_equilibrium,
                            const WholePtsType& whole_pts,
                            DihedralType& forcedihedral,
                            DihedralEnergyType& dihedralEnergy) const
  {

    vtkm::IdComponent dihedrali = dihedrallist[0];
    vtkm::IdComponent dihedralj = dihedrallist[1];
    vtkm::IdComponent dihedralk = dihedrallist[2];
    vtkm::IdComponent dihedralw = dihedrallist[3];
    vtkm::IdComponent dihedraltype = dihedral_type;
    Vec3f force_dihedrali;
    Vec3f force_dihedralj;
    Vec3f force_dihedralk;
    Vec3f force_dihedralw;

    ComputeijImproperEnergyForce(dihedrali,
                              dihedralj,
                              dihedralk,
                              dihedralw,
                              dihedraltype,
                              dihedral_coeffs_k,
                              dihedral_coeffs_equilibrium,
                              whole_pts,
                              force_dihedrali,
                              force_dihedralj,
                              force_dihedralk,
                              force_dihedralw,
                              dihedralEnergy);
    forcedihedral[0] = force_dihedrali;
    forcedihedral[1] = force_dihedralj;
    forcedihedral[2] = force_dihedralk;
    forcedihedral[3] = force_dihedralw;
  }
  template<typename WholePtsType>
  VTKM_EXEC void ComputeijImproperEnergyForce(const vtkm::IdComponent& dihedrali,
                                           const vtkm::IdComponent& dihedralj,
                                           const vtkm::IdComponent& dihedralk,
                                           const vtkm::IdComponent& dihedralw,
                                           const vtkm::IdComponent& dihedraltype,
                                           const vtkm::Vec4f& k,
                                           const Real& chi,
                                           const WholePtsType& whole_pts,
                                           Vec3f& force_dihedrali,
                                           Vec3f& force_dihedralj,
                                           Vec3f& force_dihedralk,
                                           Vec3f& force_dihedralw,
                                           Real& dihedralEnergy) const
  {

    Vec3f p_i = whole_pts.Get(dihedrali);
    Vec3f p_j = whole_pts.Get(dihedralj);
    Vec3f p_k = whole_pts.Get(dihedralk);
    Vec3f p_w = whole_pts.Get(dihedralw);

    Vec3f vb1 = p_i - p_j; 
    // minimum image distance
    if (vtkm::Abs(vb1[0]) > _box[0] / 2.0)
    {
      if (vb1[0] < 0)
        vb1[0] += _box[0];
      else
        vb1[0] -= _box[0];
    }
    if (vtkm::Abs(vb1[1]) > _box[1] / 2.0)
    {
      if (vb1[1] < 0)
        vb1[1] += _box[1];
      else
        vb1[1] -= _box[1];
    }
    if (vtkm::Abs(vb1[2]) > _box[2] / 2.0)
    {
      if (vb1[2] < 0)
        vb1[2] += _box[2];
      else
        vb1[2] -= _box[2];
    }

    Vec3f vb2 = p_k - p_j;
    //Vec3f r_kj = p_k - p_j;
    // minimum image distance
    if (vtkm::Abs(vb2[0]) > _box[0] / 2.0)
    {
      if (vb2[0] < 0)
        vb2[0] += _box[0];
      else
        vb2[0] -= _box[0];
    }
    if (vtkm::Abs(vb2[1]) > _box[1] / 2.0)
    {
      if (vb2[1] < 0)
        vb2[1] += _box[1];
      else
        vb2[1] -= _box[1];
    }
    if (vtkm::Abs(vb2[2]) > _box[2] / 2.0)
    {
      if (vb2[2] < 0)
        vb2[2] += _box[2];
      else
        vb2[2] -= _box[2];
    }
    
    Vec3f vb2m = -vb2;

    Vec3f vb3 = p_w - p_k;
    // minimum image distance
    if (vtkm::Abs(vb3[0]) > _box[0] / 2.0)
    {
      if (vb3[0] < 0)
        vb3[0] += _box[0];
      else
        vb3[0] -= _box[0];
    }
    if (vtkm::Abs(vb3[1]) > _box[1] / 2.0)
    {
      if (vb3[1] < 0)
        vb3[1] += _box[1];
      else
        vb3[1] -= _box[1];
    }
    if (vtkm::Abs(vb3[2]) > _box[2] / 2.0)
    {
      if (vb3[2] < 0)
        vb3[2] += _box[2];
      else
        vb3[2] -= _box[2];
    }

    /*sb1 = 1.0 / (vb1x * vb1x + vb1y * vb1y + vb1z * vb1z);
    sb2 = 1.0 / (vb2x * vb2x + vb2y * vb2y + vb2z * vb2z);
    sb3 = 1.0 / (vb3x * vb3x + vb3y * vb3y + vb3z * vb3z);

    rb1 = sqrt(sb1);
    rb3 = sqrt(sb3);
    */
    Real sb1 = 1.0 / vtkm::MagnitudeSquared(vb1);
    Real rb1 = vtkm::Sqrt(sb1);
    Real sb2 = 1.0 / vtkm::MagnitudeSquared(vb2);
    Real rb2 = vtkm::Sqrt(sb2);
    Real sb3 = 1.0 / vtkm::MagnitudeSquared(vb3);
    Real rb3 = vtkm::Sqrt(sb3);
   
    Real c0 = vtkm::Dot(vb1,vb3) * rb1 * rb3;
    //Real c1 = vtkm::Dot(vb1,vb2) * r1 * r2;
    //Real c2 = -vtkm::Dot(vb3,vb2) * r3 * r2;

    // 1st and 2nd angle
    Real b1mag2 = vtkm::MagnitudeSquared(vb1);
    Real b1mag = vtkm::Sqrt(b1mag2);
    Real b2mag2 = vtkm::MagnitudeSquared(vb2);
    Real b2mag = vtkm::Sqrt(b2mag2);
    Real b3mag2 = vtkm::MagnitudeSquared(vb3);
    Real b3mag = vtkm::Sqrt(b3mag2);

    Real ctmp = vtkm::Dot(vb1,vb2);
    Real r12c1 = 1.0 / (b1mag * b2mag);
    Real c1mag = ctmp * r12c1;
    ctmp = vtkm::Dot(vb2m,vb3);
    Real r12c2 = 1.0 / (b2mag * b3mag);
    Real c2mag = ctmp * r12c2;

    // cos and sin of 2 angles and final c

    Real sin2 = vtkm::Max(1.0 - c1mag * c1mag, 0.0);
    Real sc1 = vtkm::Sqrt(sin2);
    if (sc1 < SMALL) sc1 = SMALL;
    sc1 = 1.0 / sc1;
    sin2 = vtkm::Max(1.0 - c2mag * c2mag, 0.0);
    Real sc2 = vtkm::Sqrt(sin2);
    if (sc2 < SMALL) sc2 = SMALL;
    sc2 = 1.0 / sc2;

    Real s1 = sc1 * sc1;
    Real s2 = sc2 * sc2;
    Real s12 = sc1 * sc2;
    Real c = (c0 + c1mag * c2mag) * s12;
    Real cx = vb1[1] * vb2[2] - vb1[2] * vb2[1];
    Real cy = vb1[2] * vb2[0] - vb1[0] * vb2[2];
    Real cz = vb1[2] * vb2[1] - vb1[1] * vb2[0];
    Real cmag = vtkm::Sqrt(cx * cx + cy * cy + cz * cz);
    Real dx = (cx * vb3[0] + cy * vb3[1] + cz * vb3[2]) / cmag / b3mag;

    if (c > 1.0) c = 1.0;
    if (c < -1.0) c = -1.0;

    // force & energy
    // p = sum (i=1,4) k_i * (1 + (-1)**(i+1)*cos(i*phi) )
    // pd = dp/dc
    Real phi = acos(c);
    if (dx < 0.0) phi *= -1.0;
    Real si = sin(phi);
    if (vtkm::Abs(si) < 0.00001) si = 0.00001;
    Real siinv = 1.0 / si;

    Real p = k[0] * (1.0 + c) + k[1] * (1.0 - cos(2.0 * phi)) +
        k[2] * (1.0 + cos(3.0 * phi)) + k[3] * (1.0 - cos(4.0 * phi));
    Real pd = k[0] - 2.0 * k[1] * sin(2.0 * phi) * siinv +
        3.0 * k[2] * sin(3.0 * phi) * siinv - 4.0 * k[3] * sin(4.0 * phi) * siinv;

    dihedralEnergy = p;
    Real a = pd;
    c = c * a;
    s12 = s12 * a;
    Real a11 = c * sb1 * s1;
    Real a22 = -sb2 * (2.0 * c0 * s12 - c * (s1 + s2));
    Real a33 = c * sb3 * s2;
    Real a12 = -r12c1 * (c1mag * c * s1 + c2mag * s12);
    Real a13 = -rb1 * rb3 * s12;
    Real a23 = r12c2 * (c2mag * c * s2 + c1mag * s12);
    Real sx2 = a12 * vb1[0] + a22 * vb2[0] + a23 * vb3[0];
    Real sy2 = a12 * vb1[1] + a22 * vb2[1] + a23 * vb3[1];
    Real sz2 = a12 * vb1[2] + a22 * vb2[2] + a23 * vb3[2];

    force_dihedrali[0] = a11 * vb1[0] + a12 * vb2[0] + a13 * vb3[0];
    force_dihedrali[1] = a11 * vb1[1] + a12 * vb2[1] + a13 * vb3[1];
    force_dihedrali[2] = a11 * vb1[2] + a12 * vb2[2] + a13 * vb3[2];
    //f1[0] = a11 * vb1x + a12 * vb2x + a13 * vb3x;
    //f1[1] = a11 * vb1y + a12 * vb2y + a13 * vb3y;
    //f1[2] = a11 * vb1z + a12 * vb2z + a13 * vb3z;

    force_dihedralj[0] =  -sx2 - force_dihedrali[0];
    force_dihedralj[1] =  -sy2 - force_dihedrali[1];
    force_dihedralj[2] =  -sz2 - force_dihedrali[2];
    //f2[0] = -sx2 - f1[0];
    //f2[1] = -sy2 - f1[1];
    //f2[2] = -sz2 - f1[2];

    force_dihedralw[0] = a13 * vb1[0] + a23 * vb2[0] + a33 * vb3[0];
    force_dihedralw[1] = a13 * vb1[1] + a23 * vb2[1] + a33 * vb3[1];
    force_dihedralw[2] = a13 * vb1[2] + a23 * vb2[2] + a33 * vb3[2];
    //f4[0] = a13 * vb1x + a23 * vb2x + a33 * vb3x;
    //f4[1] = a13 * vb1y + a23 * vb2y + a33 * vb3y;
    //f4[2] = a13 * vb1z + a23 * vb2z + a33 * vb3z;

    force_dihedralk[0] =  sx2 - force_dihedralw[0];
    force_dihedralk[1] =  sy2 - force_dihedralw[1];
    force_dihedralk[2] =  sz2 - force_dihedralw[2];
    //f3[0] = sx2 - f4[0];
    //f3[1] = sy2 - f4[1];
    //f3[2] = sz2 - f4[2];
  }

  Vec3f _box;
};

struct ComputeDihedralHarmonicWorklet : vtkm::worklet::WorkletMapField
{

  ComputeDihedralHarmonicWorklet(const Vec3f& box)
    : _box(box)
  {
  }
  using ControlSignature = void(FieldIn dihedral_type,
                                FieldIn dihedrallist,
                                FieldIn dihedral_coeffs_k,
                                FieldIn dihedral_coeffs_sign,
                                FieldIn dihedral_coeffs_multiplicity,
                                WholeArrayIn whole_pts,
                                FieldOut forcedihedral,
                                FieldOut dihedralEnergy,
                                ExecObject locator);
  using ExecutionSignature = void(_1, _2, _3, _4, _5, _6, _7, _8, _9);

  template<typename DihedralType,
           typename DihedralListType,
           typename WholePtsType,
           typename ForceDihedralType,
           typename DihedralEnergyType>
  VTKM_EXEC void operator()(const DihedralType& dihedral_type,
                            const DihedralListType& dihedrallist,
                            const Real& dihedral_coeffs_k,
                            const vtkm::IdComponent& dihedral_coeffs_sign,
                            const vtkm::IdComponent& dihedral_coeffs_multiplicity,
                            const WholePtsType& whole_pts,
                            ForceDihedralType& forcedihedral,
                            DihedralEnergyType& dihedralEnergy,
                            const ExecPointLocator& locator) const
  {

    vtkm::IdComponent dihedrali = dihedrallist[0];
    vtkm::IdComponent dihedralj = dihedrallist[1];
    vtkm::IdComponent dihedralk = dihedrallist[2];
    vtkm::IdComponent dihedralw = dihedrallist[3];
    vtkm::IdComponent dihedraltype = dihedral_type;

    Real cos_shift, sin_shift;
    if (dihedral_coeffs_sign == 1)
    {
      cos_shift = 1.0;
      sin_shift = 0.0;
    }
    else
    {
      cos_shift = -1.0;
      sin_shift = 0.0;
    }
    Vec3f force_dihedrali;
    Vec3f force_dihedralj;
    Vec3f force_dihedralk;
    Vec3f force_dihedralw;

    ComputeijDihedralHarmonicEnergyForce(dihedrali,
                                         dihedralj,
                                         dihedralk,
                                         dihedralw,
                                         dihedraltype,
                                         dihedral_coeffs_k,
                                         dihedral_coeffs_multiplicity,
                                         cos_shift,
                                         sin_shift,
                                         whole_pts,
                                         force_dihedrali,
                                         force_dihedralj,
                                         force_dihedralk,
                                         force_dihedralw,
                                         dihedralEnergy,
                                         locator);
    forcedihedral[0] = force_dihedrali;
    forcedihedral[1] = force_dihedralj;
    forcedihedral[2] = force_dihedralk;
    forcedihedral[3] = force_dihedralw;
  }
  template<typename WholePtsType>
  VTKM_EXEC void ComputeijDihedralHarmonicEnergyForce(
    const vtkm::IdComponent& dihedrali,
    const vtkm::IdComponent& dihedralj,
    const vtkm::IdComponent& dihedralk,
    const vtkm::IdComponent& dihedralw,
    const vtkm::IdComponent& dihedraltype,
    const Real& dihedral_coeffs_k,
    const vtkm::IdComponent& dihedral_coeffs_multiplicity,
    const Real& cos_shift,
    const Real& sin_shift,
    const WholePtsType& whole_pts,
    Vec3f& force_dihedrali,
    Vec3f& force_dihedralj,
    Vec3f& force_dihedralk,
    Vec3f& force_dihedralw,
    Real& dihedralEnergy,
    const ExecPointLocator& locator) const
  {

    Vec3f p_i = whole_pts.Get(dihedrali);
    Vec3f p_j = whole_pts.Get(dihedralj);
    Vec3f p_k = whole_pts.Get(dihedralk);
    Vec3f p_w = whole_pts.Get(dihedralw);

    
    Vec3f vb1 = locator.MinDistanceVec(p_i, p_j, _box);
    Vec3f vb2 = locator.MinDistanceVec(p_k, p_j, _box);

    Vec3f vb2m = -vb2;
    Vec3f vb3 = locator.MinDistanceVec(p_w, p_k, _box);

    // c,s calculation

    Real ax = vb1[1] * vb2m[2] - vb1[2] * vb2m[1];
    Real ay = vb1[2] * vb2m[0] - vb1[0] * vb2m[2];
    Real az = vb1[0] * vb2m[1] - vb1[1] * vb2m[0];
    Real bx = vb3[1] * vb2m[2] - vb3[2] * vb2m[1];
    Real by = vb3[2] * vb2m[0] - vb3[0] * vb2m[2];
    Real bz = vb3[0] * vb2m[1] - vb3[1] * vb2m[0];

    /*ax = vb1y * vb2zm - vb1z * vb2ym;
    ay = vb1z * vb2xm - vb1x * vb2zm;
    az = vb1x * vb2ym - vb1y * vb2xm;
    bx = vb3y * vb2zm - vb3z * vb2ym;
    by = vb3z * vb2xm - vb3x * vb2zm;
    bz = vb3x * vb2ym - vb3y * vb2xm;*/

    Real rasq = ax * ax + ay * ay + az * az;
    Real rbsq = bx * bx + by * by + bz * bz;
    Real rgsq = vb2m[0] * vb2m[0] + vb2m[1] * vb2m[1] + vb2m[2] * vb2m[2];
    Real rg = vtkm::Sqrt(rgsq);
    /*rasq = ax * ax + ay * ay + az * az;
    rbsq = bx * bx + by * by + bz * bz;
    rgsq = vb2xm * vb2xm + vb2ym * vb2ym + vb2zm * vb2zm;
    rg = sqrt(rgsq);*/

    Real rginv, ra2inv, rb2inv;
    rginv = ra2inv = rb2inv = 0.0;
    if (rg > 0)
    {
      rginv = 1.0 / rg;
    }
    if (rasq > 0)
    {
      ra2inv = 1.0 / rasq;
    }
    if (rbsq > 0)
    {
      rb2inv = 1.0 / rbsq;
    }
    Real rabinv = vtkm::Sqrt(ra2inv * rb2inv);

    Real c = (ax * bx + ay * by + az * bz) * rabinv;
    Real s = rg * rabinv * (ax * vb3[0] + ay * vb3[1] + az * vb3[2]);
    //s = rg * rabinv * (ax * vb3x + ay * vb3y + az * vb3z);

    // error check

    //if (c > 1.0 + TOLERANCE || c < (-1.0 - TOLERANCE))
    //  problem(FLERR, i1, i2, i3, i4);

    if (c > 1.0)
    {
      c = 1.0;
    }
    if (c < -1.0)
    {
      c = -1.0;
    }

    vtkm::IdComponent m = dihedral_coeffs_multiplicity;
    Real p = 1.0;
    Real ddf1, df1;
    ddf1 = df1 = 0.0;

    /*m = multiplicity[type];
    p = 1.0;
    ddf1 = df1 = 0.0;*/

    for (vtkm::IdComponent i = 0; i < m; i++)
    {
      ddf1 = p * c - df1 * s;
      df1 = p * s + df1 * c;
      p = ddf1;
    }

    p = p * cos_shift + df1 * sin_shift;
    df1 = df1 * cos_shift - ddf1 * sin_shift;
    /*p = p * cos_shift[type] + df1 * sin_shift[type];
    df1 = df1 * cos_shift[type] - ddf1 * sin_shift[type];*/
    df1 *= -m;
    p += 1.0;

    if (m == 0)
    {
      p = 1.0 + cos_shift;
      //p = 1.0 + cos_shift[type];
      df1 = 0.0;
    }

    dihedralEnergy = dihedral_coeffs_k * p;
    //if (eflag)
    //  edihedral = k[type] * p;

    Real fg = vb1[0] * vb2m[0] + vb1[1] * vb2m[1] + vb1[2] * vb2m[2];
    Real hg = vb3[0] * vb2m[0] + vb3[1] * vb2m[1] + vb3[2] * vb2m[2];
    //fg = vb1x * vb2xm + vb1y * vb2ym + vb1z * vb2zm;
    //hg = vb3x * vb2xm + vb3y * vb2ym + vb3z * vb2zm;
    Real fga = fg * ra2inv * rginv;
    Real hgb = hg * rb2inv * rginv;
    Real gaa = -ra2inv * rg;
    Real gbb = rb2inv * rg;

    Real dtfx = gaa * ax;
    Real dtfy = gaa * ay;
    Real dtfz = gaa * az;
    Real dtgx = fga * ax - hgb * bx;
    Real dtgy = fga * ay - hgb * by;
    Real dtgz = fga * az - hgb * bz;
    Real dthx = gbb * bx;
    Real dthy = gbb * by;
    Real dthz = gbb * bz;

    Real df = -dihedral_coeffs_k * df1;
    //df = -k[type] * df1;

    Real sx2 = df * dtgx;
    Real sy2 = df * dtgy;
    Real sz2 = df * dtgz;


    force_dihedrali[0] = df * dtfx;
    force_dihedrali[1] = df * dtfy;
    force_dihedrali[2] = df * dtfz;

    /*f1[0] = df * dtfx;
    f1[1] = df * dtfy;
    f1[2] = df * dtfz;*/

    force_dihedralj[0] = sx2 - force_dihedrali[0];
    force_dihedralj[1] = sy2 - force_dihedrali[1];
    force_dihedralj[2] = sz2 - force_dihedrali[2];

    //f2[0] = sx2 - f1[0];
    //f2[1] = sy2 - f1[1];
    //f2[2] = sz2 - f1[2];

    force_dihedralw[0] = df * dthx;
    force_dihedralw[1] = df * dthy;
    force_dihedralw[2] = df * dthz;

    //f4[0] = df * dthx;
    //f4[1] = df * dthy;
    //f4[2] = df * dthz;

    force_dihedralk[0] = -sx2 - force_dihedralw[0];
    force_dihedralk[1] = -sy2 - force_dihedralw[1];
    force_dihedralk[2] = -sz2 - force_dihedralw[2];

    //f3[0] = -sx2 - f4[0];
    //f3[1] = -sy2 - f4[1];
    //f3[2] = -sz2 - f4[2];
  }

   Vec3f _box;
};

struct ComputeImproperHarmonicWorklet : vtkm::worklet::WorkletMapField
{

  ComputeImproperHarmonicWorklet(const Vec3f& box)
     : _box(box)
  {
  }
  using ControlSignature = void(FieldIn improper_type,
                                const FieldIn improperlist,
                                FieldIn improper_coeffs_k,
                                FieldIn improper_coeffs_equilibrium,
                                const WholeArrayIn whole_pts,
                                FieldOut forceimproper,
                                FieldOut improperEnergy);
  using ExecutionSignature = void(_1, _2, _3, _4, _5, _6, _7);

  template<typename ImproperType,
           typename ImproperListType,
           typename WholePtsType,
           typename ForceImproperType,
           typename ImproperEnergyType>
  VTKM_EXEC void operator()(const ImproperType& improper_type,
                            const ImproperListType& improperlist,
                            const Real& improper_coeffs_k,
                            const Real& improper_coeffs_equilibrium,
                            const WholePtsType& whole_pts,
                            ForceImproperType& forceimproper,
                            ImproperEnergyType& improperEnergy) const
  {

    vtkm::IdComponent improperi = improperlist[0];
    vtkm::IdComponent improperj = improperlist[1];
    vtkm::IdComponent improperk = improperlist[2];
    vtkm::IdComponent improperw = improperlist[3];
    vtkm::IdComponent impropertype = improper_type;
    Vec3f force_improperi;
    Vec3f force_improperj;
    Vec3f force_improperk;
    Vec3f force_improperw;

    ComputeijImproperEnergyForce(improperi,
                              improperj,
                              improperk,
                              improperw,
                              impropertype,
                              improper_coeffs_k,
                              improper_coeffs_equilibrium,
                              whole_pts,
                              force_improperi,
                              force_improperj,
                              force_improperk,
                              force_improperw,
                              improperEnergy);
    forceimproper[0] = force_improperi;
    forceimproper[1] = force_improperj;
    forceimproper[2] = force_improperk;
    forceimproper[3] = force_improperw;
  }
  template<typename WholePtsType>
  VTKM_EXEC void ComputeijImproperEnergyForce(const vtkm::IdComponent& improperi,
                                           const vtkm::IdComponent& improperj,
                                           const vtkm::IdComponent& improperk,
                                           const vtkm::IdComponent& improperw,
                                           const vtkm::IdComponent& impropertype,
                                           const Real& k,
                                           const Real& chi,
                                           const WholePtsType& whole_pts,
                                           Vec3f& force_improperi,
                                           Vec3f& force_improperj,
                                           Vec3f& force_improperk,
                                           Vec3f& force_improperw,
                                           Real& improperEnergy) const
  {

    Vec3f p_i = whole_pts.Get(improperi);
    Vec3f p_j = whole_pts.Get(improperj);
    Vec3f p_k = whole_pts.Get(improperk);
    Vec3f p_w = whole_pts.Get(improperw);

    Vec3f vb1 = p_i - p_j; 
    // minimum image distance
    if (vtkm::Abs(vb1[0]) > _box[0] / 2.0)
    {
      if (vb1[0] < 0)
        vb1[0] += _box[0];
      else
        vb1[0] -= _box[0];
    }
    if (vtkm::Abs(vb1[1]) > _box[1] / 2.0)
    {
      if (vb1[1] < 0)
        vb1[1] += _box[1];
      else
        vb1[1] -= _box[1];
    }
    if (vtkm::Abs(vb1[2]) > _box[2] / 2.0)
    {
      if (vb1[2] < 0)
        vb1[2] += _box[2];
      else
        vb1[2] -= _box[2];
    }

    Vec3f vb2 = p_k - p_j;
    //Vec3f r_kj = p_k - p_j;
    // minimum image distance
    if (vtkm::Abs(vb2[0]) > _box[0] / 2.0)
    {
      if (vb2[0] < 0)
        vb2[0] += _box[0];
      else
        vb2[0] -= _box[0];
    }
    if (vtkm::Abs(vb2[1]) > _box[1] / 2.0)
    {
      if (vb2[1] < 0)
        vb2[1] += _box[1];
      else
        vb2[1] -= _box[1];
    }
    if (vtkm::Abs(vb2[2]) > _box[2] / 2.0)
    {
      if (vb2[2] < 0)
        vb2[2] += _box[2];
      else
        vb2[2] -= _box[2];
    }


    Vec3f vb3 = p_w - p_k;
    // minimum image distance
    if (vtkm::Abs(vb3[0]) > _box[0] / 2.0)
    {
      if (vb3[0] < 0)
        vb3[0] += _box[0];
      else
        vb3[0] -= _box[0];
    }
    if (vtkm::Abs(vb3[1]) > _box[1] / 2.0)
    {
      if (vb3[1] < 0)
        vb3[1] += _box[1];
      else
        vb3[1] -= _box[1];
    }
    if (vtkm::Abs(vb3[2]) > _box[2] / 2.0)
    {
      if (vb3[2] < 0)
        vb3[2] += _box[2];
      else
        vb3[2] -= _box[2];
    }

    /*ss1 = 1.0 / (vb1x * vb1x + vb1y * vb1y + vb1z * vb1z);
    ss2 = 1.0 / (vb2x * vb2x + vb2y * vb2y + vb2z * vb2z);
    ss3 = 1.0 / (vb3x * vb3x + vb3y * vb3y + vb3z * vb3z);

    r1 = sqrt(ss1);
    r2 = sqrt(ss2);
    r3 = sqrt(ss3);
    */
    Real ss1 = 1.0 / vtkm::MagnitudeSquared(vb1);
    Real r1 = vtkm::Sqrt(ss1);
    Real ss2 = 1.0 / vtkm::MagnitudeSquared(vb2);
    Real r2 = vtkm::Sqrt(ss2);
    Real ss3 = 1.0 / vtkm::MagnitudeSquared(vb3);
    Real r3 = vtkm::Sqrt(ss3);
   
    Real c0 = vtkm::Dot(vb1,vb3) * r1 * r3;
    Real c1 = vtkm::Dot(vb1,vb2) * r1 * r2;
    Real c2 = -vtkm::Dot(vb3,vb2) * r3 * r2;

    Real s1 = 1.0 - c1 * c1;
    if (s1 < SMALL) 
    {
      s1 = SMALL;
    }
    s1 = 1.0 / s1;
    Real s2 = 1.0 - c2 * c2;
    if (s2 < SMALL) 
    {
      s2 = SMALL;
    }
    s2 = 1.0 / s2;
    Real s12 = sqrt(s1 * s2);
    Real c = (c1 * c2 + c0) * s12;

    if (c > 1.0) c = 1.0;
    if (c < -1.0) c = -1.0;
    Real s = sqrt(1.0 - c * c);
    if (s < SMALL) 
    {
      s = SMALL;
    }

    Real domega = acos(c) - (chi * vtkm::Pif()) / 180;
    Real a_improper = k * domega;
    improperEnergy = a_improper * domega;

    a_improper = -a_improper * 2.0 / s;
    c = c * a_improper;
    s12 = s12 * a_improper;
    Real a11 = c * ss1 * s1;
    Real a22 = -ss2 * (2.0 * c0 * s12 - c * (s1 + s2));
    Real a33 = c * ss3 * s2;
    Real a12 = -r1 * r2 * (c1 * c * s1 + c2 * s12);
    Real a13 = -r1 * r3 * s12;
    Real a23 = r2 * r3 * (c2 * c * s2 + c1 * s12);

    Real sx2 = a22 * vb2[0] + a23 * vb3[0] + a12 * vb1[0];
    Real sy2 = a22 * vb2[1] + a23 * vb3[1] + a12 * vb1[1];
    Real sz2 = a22 * vb2[2] + a23 * vb3[2] + a12 * vb1[2];

    //sx2 = a22 * vb2x + a23 * vb3x + a12 * vb1x;
    //sy2 = a22 * vb2y + a23 * vb3y + a12 * vb1y;
    //sz2 = a22 * vb2z + a23 * vb3z + a12 * vb1z;

    force_improperi[0] = a12 * vb2[0] + a13 * vb3[0] + a11 * vb1[0];
    force_improperi[1] = a12 * vb2[1] + a13 * vb3[1] + a11 * vb1[1];
    force_improperi[2] = a12 * vb2[2] + a13 * vb3[2] + a11 * vb1[2];
    //f1[0] = a12 * vb2x + a13 * vb3x + a11 * vb1x;
    //f1[1] = a12 * vb2y + a13 * vb3y + a11 * vb1y;
    //f1[2] = a12 * vb2z + a13 * vb3z + a11 * vb1z;

    force_improperj[0] =  -sx2 - force_improperi[0];
    force_improperj[1] =  -sy2 - force_improperi[1];
    force_improperj[2] =  -sz2 - force_improperi[2];
    //f2[0] = -sx2 - f1[0];
    //f2[1] = -sy2 - f1[1];
    //f2[2] = -sz2 - f1[2];

    force_improperw[0] = a23 * vb2[0] + a33 * vb3[0] + a13 * vb1[0];
    force_improperw[1] = a23 * vb2[1] + a33 * vb3[1] + a13 * vb1[1];
    force_improperw[2] = a23 * vb2[2] + a33 * vb3[2] + a13 * vb1[2];
    //f4[0] = a23 * vb2x + a33 * vb3x + a13 * vb1x;
    //f4[1] = a23 * vb2y + a33 * vb3y + a13 * vb1y;
    //f4[2] = a23 * vb2z + a33 * vb3z + a13 * vb1z;

    force_improperk[0] =  sx2 - force_improperw[0];
    force_improperk[1] =  sy2 - force_improperw[1];
    force_improperk[2] =  sz2 - force_improperw[2];
    //f3[0] = sx2 - f4[0];
    //f3[1] = sy2 - f4[1];
    //f3[2] = sz2 - f4[2];
  }

   Vec3f _box;
};

// "ConstraintWaterBondAngleWorklet"未使用
struct ConstraintWaterBondAngleWorklet : vtkm::worklet::WorkletMapField
{
  ConstraintWaterBondAngleWorklet(const Real& vlength,
                                  const Real& dt,
                                  const Real fmt2v,
                                  const Vec<Vec2f, 3>& range)
    : _vlength(vlength)
    , _dt(dt)
    , _fmt2v(fmt2v)
    , _range(range)
  {
  }
  using ControlSignature = void(const FieldIn anglelist,
                                WholeArrayInOut whole_pts,
                                WholeArrayInOut whole_velocity,
                                const WholeArrayIn whole_all_force,
                                const WholeArrayIn whole_mass,
                                FieldOut test_angle);
  using ExecutionSignature = void(_1, _2, _3, _4, _5, _6);

  template<typename AngleListType,
           typename PositionType,
           typename VelocityType,
           typename AllForceType,
           typename MassType>
  VTKM_EXEC void operator()(const AngleListType& anglelist,
                            PositionType& whole_pts,
                            VelocityType& whole_velocity,
                            const AllForceType& whole_all_force,
                            const MassType& whole_mass,
                            Vec2f& test_angle) const
  {
    IdComponent i0 = anglelist[1];
    IdComponent i1 = anglelist[0];
    IdComponent i2 = anglelist[2];

    vtkm::Vec<Vec3f, 3> shake_position;
    shake_position[0] = whole_pts.Get(i0) + whole_velocity.Get(i0) * _dt +
      0.5 * _dt * _dt * whole_all_force.Get(i0) / whole_mass.Get(i0) * _fmt2v;
    shake_position[1] = whole_pts.Get(i1) + whole_velocity.Get(i1) * _dt +
      0.5 * _dt * _dt * whole_all_force.Get(i1) / whole_mass.Get(i1) * _fmt2v;
    shake_position[2] = whole_pts.Get(i2) + whole_velocity.Get(i2) * _dt +
      0.5 * _dt * _dt * whole_all_force.Get(i2) / whole_mass.Get(i2) * _fmt2v;

    Real bond1 = 1.0;
    Real bond2 = 1.0;
    Real bond12 = vtkm::Sqrt(bond1 * bond1 + bond2 * bond2 -
                             2.0 * bond1 * bond2 * vtkm::Cos((109.4700 / 180.0) * vtkm::Pif()));

    // minimum image
    Vec3f r01;
    r01 = whole_pts.Get(i0) - whole_pts.Get(i1);

    while (vtkm::Abs(r01[0]) > _vlength / 2.0)
    {
      if (r01[0] < 0)
        r01[0] = r01[0] + _vlength;
      else
        r01[0] = r01[0] - _vlength;
    }
    while (vtkm::Abs(r01[1]) > _vlength / 2.0)
    {
      if (r01[1] < 0)
        r01[1] = r01[1] + _vlength;
      else
        r01[1] = r01[1] - _vlength;
    }
    while (vtkm::Abs(r01[2]) > _vlength / 2.0)
    {
      if (r01[2] < 0)
        r01[2] = r01[2] + _vlength;
      else
        r01[2] = r01[2] - _vlength;
    }
    Vec3f r02;
    r02 = whole_pts.Get(i0) - whole_pts.Get(i2);
    //domain->minimum_image(r02);

    while (vtkm::Abs(r02[0]) > _vlength / 2.0)
    {
      if (r02[0] < 0)
        r02[0] = r02[0] + _vlength;
      else
        r02[0] = r02[0] - _vlength;
    }
    while (vtkm::Abs(r02[1]) > _vlength / 2.0)
    {
      if (r02[1] < 0)
        r02[1] = r02[1] + _vlength;
      else
        r02[1] = r02[1] - _vlength;
    }
    while (vtkm::Abs(r02[2]) > _vlength / 2.0)
    {
      if (r02[2] < 0)
        r02[2] = r02[2] + _vlength;
      else
        r02[2] = r02[2] - _vlength;
    }
    Vec3f r12;
    r12 = whole_pts.Get(i1) - whole_pts.Get(i2);
    //domain->minimum_image(r12);

    while (vtkm::Abs(r12[0]) > _vlength / 2.0)
    {
      if (r12[0] < 0)
        r12[0] = r12[0] + _vlength;
      else
        r12[0] = r12[0] - _vlength;
    }
    while (vtkm::Abs(r12[1]) > _vlength / 2.0)
    {
      if (r12[1] < 0)
        r12[1] = r12[1] + _vlength;
      else
        r12[1] = r12[1] - _vlength;
    }
    while (vtkm::Abs(r12[2]) > _vlength / 2.0)
    {
      if (r12[2] < 0)
        r12[2] = r12[2] + _vlength;
      else
        r12[2] = r12[2] - _vlength;
    }
    // s01,s02,s12 = distance vec after unconstrained update, with PBC

    Vec3f s01;
    s01 = shake_position[0] - shake_position[1];
    //domain->minimum_image_once(s01);

    Vec3f s02;
    s02 = shake_position[0] - shake_position[2];
    //domain->minimum_image_once(s02);

    Vec3f s12;
    s12 = shake_position[1] - shake_position[2];

    if (vtkm::Abs(s01[0]) > _vlength / 2.0)
    {
      if (s01[0] < 0)
        s01[0] = s01[0] + _vlength;
      else
        s01[0] = s01[0] - _vlength;
    }
    if (vtkm::Abs(s01[1]) > _vlength / 2.0)
    {
      if (s01[1] < 0)
        s01[1] = s01[1] + _vlength;
      else
        s01[1] = s01[1] - _vlength;
    }
    if (vtkm::Abs(s01[2]) > _vlength / 2.0)
    {
      if (s01[2] < 0)
        s01[2] = s01[2] + _vlength;
      else
        s01[2] = s01[2] - _vlength;
    }

    if (vtkm::Abs(s02[0]) > _vlength / 2.0)
    {
      if (s02[0] < 0)
        s02[0] = s02[0] + _vlength;
      else
        s02[0] = s02[0] - _vlength;
    }
    if (vtkm::Abs(s02[1]) > _vlength / 2.0)
    {
      if (s02[1] < 0)
        s02[1] = s02[1] + _vlength;
      else
        s02[1] = s02[1] - _vlength;
    }
    if (vtkm::Abs(s02[2]) > _vlength / 2.0)
    {
      if (s02[2] < 0)
        s02[2] = s02[2] + _vlength;
      else
        s02[2] = s02[2] - _vlength;
    }

    if (vtkm::Abs(s12[0]) > _vlength / 2.0)
    {
      if (s12[0] < 0)
        s12[0] = s12[0] + _vlength;
      else
        s12[0] = s12[0] - _vlength;
    }
    if (vtkm::Abs(s12[1]) > _vlength / 2.0)
    {
      if (s12[1] < 0)
        s12[1] = s12[1] + _vlength;
      else
        s12[1] = s12[1] - _vlength;
    }
    if (vtkm::Abs(s12[2]) > _vlength / 2.0)
    {
      if (s12[2] < 0)
        s12[2] = s12[2] + _vlength;
      else
        s12[2] = s12[2] - _vlength;
    }

    // scalar distances between atoms

    Real r01sq = vtkm::Dot(r01, r01);
    Real r02sq = vtkm::Dot(r02, r02);
    Real r12sq = vtkm::Dot(r12, r12);
    Real s01sq = vtkm::Dot(s01, s01);
    Real s02sq = vtkm::Dot(s02, s02);
    Real s12sq = vtkm::Dot(s12, s12);

    // matrix coeffs and rhs for lamda equations

    Real invmass0 = 1.0 / whole_mass.Get(i0);
    Real invmass1 = 1.0 / whole_mass.Get(i1);
    Real invmass2 = 1.0 / whole_mass.Get(i2);
    Real a11 = 2.0 * (invmass0 + invmass1) * vtkm::Dot(s01, r01);
    Real a12 = 2.0 * invmass0 * vtkm::Dot(s01, r02);
    Real a13 = -2.0 * invmass1 * vtkm::Dot(s01, r12);
    Real a21 = 2.0 * invmass0 * vtkm::Dot(s02, r01);
    Real a22 = 2.0 * (invmass0 + invmass2) * vtkm::Dot(s02, r02);
    Real a23 = 2.0 * invmass2 * vtkm::Dot(s02, r12);
    Real a31 = -2.0 * invmass1 * vtkm::Dot(s12, r01);
    Real a32 = 2.0 * invmass2 * vtkm::Dot(s12, r02);
    Real a33 = 2.0 * (invmass1 + invmass2) * (vtkm::Dot(s12, r12));

    // inverse of matrix

    Real determ = a11 * a22 * a33 + a12 * a23 * a31 + a13 * a21 * a32 - a11 * a23 * a32 -
      a12 * a21 * a33 - a13 * a22 * a31;
    if (vtkm::Abs(determ) < 0.0001)
      printf("Shake determinant = 0.0");

    Real determinv = 1.0 / determ;

    Real a11inv = determinv * (a22 * a33 - a23 * a32);
    Real a12inv = -determinv * (a12 * a33 - a13 * a32);
    Real a13inv = determinv * (a12 * a23 - a13 * a22);
    Real a21inv = -determinv * (a21 * a33 - a23 * a31);
    Real a22inv = determinv * (a11 * a33 - a13 * a31);
    Real a23inv = -determinv * (a11 * a23 - a13 * a21);
    Real a31inv = determinv * (a21 * a32 - a22 * a31);
    Real a32inv = -determinv * (a11 * a32 - a12 * a31);
    Real a33inv = determinv * (a11 * a22 - a12 * a21);

    // quadratic correction coeffs

    Real r0102 = vtkm::Dot(r01, r02);
    Real r0112 = vtkm::Dot(r01, r12);
    Real r0212 = vtkm::Dot(r02, r12);

    Real quad1_0101 = (invmass0 + invmass1) * (invmass0 + invmass1) * r01sq;
    Real quad1_0202 = invmass0 * invmass0 * r02sq;
    Real quad1_1212 = invmass1 * invmass1 * r12sq;
    Real quad1_0102 = 2.0 * (invmass0 + invmass1) * invmass0 * r0102;
    Real quad1_0112 = -2.0 * (invmass0 + invmass1) * invmass1 * r0112;
    Real quad1_0212 = -2.0 * invmass0 * invmass1 * r0212;
    Real quad2_0101 = invmass0 * invmass0 * r01sq;
    Real quad2_0202 = (invmass0 + invmass2) * (invmass0 + invmass2) * r02sq;
    Real quad2_1212 = invmass2 * invmass2 * r12sq;
    Real quad2_0102 = 2.0 * (invmass0 + invmass2) * invmass0 * r0102;
    Real quad2_0112 = 2.0 * invmass0 * invmass2 * r0112;
    Real quad2_0212 = 2.0 * (invmass0 + invmass2) * invmass2 * r0212;
    Real quad3_0101 = invmass1 * invmass1 * r01sq;
    Real quad3_0202 = invmass2 * invmass2 * r02sq;
    Real quad3_1212 = (invmass1 + invmass2) * (invmass1 + invmass2) * r12sq;
    Real quad3_0102 = -2.0 * invmass1 * invmass2 * r0102;
    Real quad3_0112 = -2.0 * (invmass1 + invmass2) * invmass1 * r0112;
    Real quad3_0212 = 2.0 * (invmass1 + invmass2) * invmass2 * r0212;

    // iterate until converged


    Real tolerance = 0.0001; // original 0.001

    IdComponent max_iter = 100; // original: 100


    Real lamda01 = 0.0;
    Real lamda02 = 0.0;
    Real lamda12 = 0.0;
    IdComponent niter = 0;
    IdComponent done = 0;
    IdComponent flag_overflow = 0;
    Real quad1, quad2, quad3, b1, b2, b3, lamda01_new, lamda02_new, lamda12_new;

    while (!done && niter < max_iter)
    {

      quad1 = quad1_0101 * lamda01 * lamda01 + quad1_0202 * lamda02 * lamda02 +
        quad1_1212 * lamda12 * lamda12 + quad1_0102 * lamda01 * lamda02 +
        quad1_0112 * lamda01 * lamda12 + quad1_0212 * lamda02 * lamda12;

      quad2 = quad2_0101 * lamda01 * lamda01 + quad2_0202 * lamda02 * lamda02 +
        quad2_1212 * lamda12 * lamda12 + quad2_0102 * lamda01 * lamda02 +
        quad2_0112 * lamda01 * lamda12 + quad2_0212 * lamda02 * lamda12;

      quad3 = quad3_0101 * lamda01 * lamda01 + quad3_0202 * lamda02 * lamda02 +
        quad3_1212 * lamda12 * lamda12 + quad3_0102 * lamda01 * lamda02 +
        quad3_0112 * lamda01 * lamda12 + quad3_0212 * lamda02 * lamda12;

      b1 = bond1 * bond1 - s01sq - quad1;
      b2 = bond2 * bond2 - s02sq - quad2;
      b3 = bond12 * bond12 - s12sq - quad3;

      lamda01_new = a11inv * b1 + a12inv * b2 + a13inv * b3;
      lamda02_new = a21inv * b1 + a22inv * b2 + a23inv * b3;
      lamda12_new = a31inv * b1 + a32inv * b2 + a33inv * b3;

      done = 1;
      if (vtkm::Abs(lamda01_new - lamda01) > tolerance)
        done = 0;
      if (vtkm::Abs(lamda02_new - lamda02) > tolerance)
        done = 0;
      if (vtkm::Abs(lamda12_new - lamda12) > tolerance)
        done = 0;

      //std::cout << "In process: lambda01 = " << lamda01 << std::endl;

      lamda01 = lamda01_new;
      lamda02 = lamda02_new;
      lamda12 = lamda12_new;

      // stop iterations before we have a floating point overflow

      // max double is < 1.0e308, so 1e150 is a reasonable cutoff


      if (vtkm::IsNan(lamda01) || vtkm::IsNan(lamda02)||vtkm::IsNan(lamda12) || 
          vtkm::Abs(lamda01) > 1e20 || vtkm::Abs(lamda02) > 1e20 || vtkm::Abs(lamda12) > 1e20)
      {
        done = 1;
        flag_overflow = 1;
      }
      niter++;
    }

    Vec3f position_constraint_i0;
    Vec3f position_constraint_i1;
    Vec3f position_constraint_i2;

    if (!flag_overflow)
    {
      position_constraint_i0 = lamda01 * r01 * invmass0 + lamda02 * r02 * invmass0;
      position_constraint_i1 = -lamda01 * r01 * invmass1 + lamda12 * r12 * invmass1;
      position_constraint_i2 = -lamda02 * r02 * invmass2 - lamda12 * r12 * invmass2;
    }
    else
    {
      position_constraint_i0 = { 0.0, 0.0, 0.0 };
      position_constraint_i1 = { 0.0, 0.0, 0.0 };
      position_constraint_i2 = { 0.0, 0.0, 0.0 };
    }

    Vec3f velocity_constraint_i0 = position_constraint_i0 / _dt;
    Vec3f velocity_constraint_i1 = position_constraint_i1 / _dt;
    Vec3f velocity_constraint_i2 = position_constraint_i2 / _dt;

    vtkm::Vec<Vec3f, 3> shake_velocity;

    shake_velocity[0] = whole_velocity.Get(i0) +
      0.5 * _dt * whole_all_force.Get(i0) / whole_mass.Get(i0) * _fmt2v;
    shake_velocity[1] = whole_velocity.Get(i1) +
      0.5 * _dt * whole_all_force.Get(i1) / whole_mass.Get(i1) * _fmt2v;
    shake_velocity[2] = whole_velocity.Get(i2) +
      0.5 * _dt * whole_all_force.Get(i2) / whole_mass.Get(i2) * _fmt2v;

    whole_velocity.Set(i0, shake_velocity[0] + velocity_constraint_i0);
    whole_velocity.Set(i1, shake_velocity[1] + velocity_constraint_i1);
    whole_velocity.Set(i2, shake_velocity[2] + velocity_constraint_i2);

    // -----------------------------------------------------------------------
    // test angle before shake
    Vec3f r_ij = whole_pts.Get(anglelist[0]) - whole_pts.Get(anglelist[1]);
    if (vtkm::Abs(r_ij[0]) > _vlength / 2.0)
    {
      if (r_ij[0] < 0)
        r_ij[0] += _vlength;
      else
        r_ij[0] -= _vlength;
    }
    if (vtkm::Abs(r_ij[1]) > _vlength / 2.0)
    {
      if (r_ij[1] < 0)
        r_ij[1] += _vlength;
      else
        r_ij[1] -= _vlength;
    }
    if (vtkm::Abs(r_ij[2]) > _vlength / 2.0)
    {
      if (r_ij[2] < 0)
        r_ij[2] += _vlength;
      else
        r_ij[2] -= _vlength;
    }
    Real disij = vtkm::Magnitude(r_ij);

    Vec3f r_kj = whole_pts.Get(anglelist[2]) - whole_pts.Get(anglelist[1]);
    if (vtkm::Abs(r_kj[0]) > _vlength / 2.0)
    {
      if (r_kj[0] < 0)
        r_kj[0] += _vlength;
      else
        r_kj[0] -= _vlength;
    }
    if (vtkm::Abs(r_kj[1]) > _vlength / 2.0)
    {
      if (r_kj[1] < 0)
        r_kj[1] += _vlength;
      else
        r_kj[1] -= _vlength;
    }
    if (vtkm::Abs(r_kj[2]) > _vlength / 2.0)
    {
      if (r_kj[2] < 0)
        r_kj[2] += _vlength;
      else
        r_kj[2] -= _vlength;
    }
    Real diskj = vtkm::Magnitude(r_kj);
    Real cosangle = vtkm::Dot(r_ij, r_kj);
    cosangle /= disij * diskj;
    if (cosangle > 1.0)
      cosangle = 1.0;
    if (cosangle < -1.0)
      cosangle = -1.0;
    test_angle[0] = acos(cosangle) / vtkm::Pif() * 180;
    //--------------------------------------------------------------------------------


    whole_pts.Set(i0, shake_position[0] + position_constraint_i0);
    whole_pts.Set(i1, shake_position[1] + position_constraint_i1);
    whole_pts.Set(i2, shake_position[2] + position_constraint_i2);

    vtkm::Vec<Vec3f, 3> whole_position_pbc;
    whole_position_pbc[0] = whole_pts.Get(anglelist[0]);
    whole_position_pbc[1] = whole_pts.Get(anglelist[1]);
    whole_position_pbc[2] = whole_pts.Get(anglelist[2]);

    for (auto i = 0; i < 3; i++)
    {
      for (auto j = 0; j < 3; j++)
      {
        if (whole_position_pbc[i][j] < _range[j][0])
        {
          whole_position_pbc[i][j] += (_range[j][1] - _range[j][0]);
        }
        if (whole_position_pbc[i][j] > _range[j][1])
        {
          whole_position_pbc[i][j] -= (_range[j][1] - _range[j][0]);
        }
      }
    }

    for (int i = 0; i < 3; i++)
    {
      whole_pts.Set(anglelist[i], whole_position_pbc[i]);
    }

    //check

    Real db1 = vtkm::Abs(vtkm::Sqrt(vtkm::Dot(r01, r01)) - bond1);
    Real db2 = vtkm::Abs(vtkm::Sqrt(vtkm::Dot(r02, r02)) - bond2);
    Real db12 = vtkm::Abs(vtkm::Sqrt(vtkm::Dot(r12, r12)) - bond12);

    // -----------------------------------------------------------------------------------
    // test angle after shake
    r_ij = whole_pts.Get(anglelist[0]) - whole_pts.Get(anglelist[1]);
    if (vtkm::Abs(r_ij[0]) > _vlength / 2.0)
    {
      if (r_ij[0] < 0)
        r_ij[0] += _vlength;
      else
        r_ij[0] -= _vlength;
    }
    if (vtkm::Abs(r_ij[1]) > _vlength / 2.0)
    {
      if (r_ij[1] < 0)
        r_ij[1] += _vlength;
      else
        r_ij[1] -= _vlength;
    }
    if (vtkm::Abs(r_ij[2]) > _vlength / 2.0)
    {
      if (r_ij[2] < 0)
        r_ij[2] += _vlength;
      else
        r_ij[2] -= _vlength;
    }
    disij = vtkm::Magnitude(r_ij);

    r_kj = whole_pts.Get(anglelist[2]) - whole_pts.Get(anglelist[1]);
    if (vtkm::Abs(r_kj[0]) > _vlength / 2.0)
    {
      if (r_kj[0] < 0)
        r_kj[0] += _vlength;
      else
        r_kj[0] -= _vlength;
    }
    if (vtkm::Abs(r_kj[1]) > _vlength / 2.0)
    {
      if (r_kj[1] < 0)
        r_kj[1] += _vlength;
      else
        r_kj[1] -= _vlength;
    }
    if (vtkm::Abs(r_kj[2]) > _vlength / 2.0)
    {
      if (r_kj[2] < 0)
        r_kj[2] += _vlength;
      else
        r_kj[2] -= _vlength;
    }
    diskj = vtkm::Magnitude(r_kj);
    cosangle = vtkm::Dot(r_ij, r_kj);
    cosangle /= disij * diskj;
    if (cosangle > 1.0)
      cosangle = 1.0;
    if (cosangle < -1.0)
      cosangle = -1.0;
    test_angle[1] = acos(cosangle) / vtkm::Pif() * 180;

    // --------------------------------------------------------------------------
  }

  Real _vlength;
  Real _dt;
  Vec<Vec2f, 3> _range;
  Real _fmt2v;
};

struct ConstraintWaterVelocityBondAngleWorklet : vtkm::worklet::WorkletMapField
{
  ConstraintWaterVelocityBondAngleWorklet(const Vec3f& box,
                                          const Real& dt,
                                          const Real fmt2v,
                                          const Vec<Vec2f, 3>& range)
    : _box(box)
    , _dt(dt)
    , _fmt2v(fmt2v)
    , _range(range)
  {
  }
  using ControlSignature = void(const FieldIn anglelist,
                                const WholeArrayIn whole_pts,
                                WholeArrayInOut whole_velocity,
                                const WholeArrayIn whole_all_force,
                                const WholeArrayIn whole_mass,
                                ExecObject locator);
  using ExecutionSignature = void(_1, _2, _3, _4, _5, _6);

  template<typename AngleListType,
           typename PositionType,
           typename VelocityType,
           typename AllForceType,
           typename MassType>
  VTKM_EXEC void operator()(const AngleListType& anglelist,
                            const PositionType& whole_pts,
                            VelocityType& whole_velocity,
                            const AllForceType& whole_all_force,
                            const MassType& whole_mass,
                            const ExecPointLocator& locator) const
  {

    // here shake_velocity is the vp in Lammps

    vtkm::Vec<vtkm::Vec3f, 3> shake_velocity;
    auto i0 = anglelist[1];
    auto i1 = anglelist[0];
    auto i2 = anglelist[2];

    shake_velocity[0] = whole_velocity.Get(i0) +
      0.5 * _dt * whole_all_force.Get(i0) / whole_mass.Get(i0) * _fmt2v;
    shake_velocity[1] = whole_velocity.Get(i1) +
      0.5 * _dt * whole_all_force.Get(i1) / whole_mass.Get(i1) * _fmt2v;
    shake_velocity[2] = whole_velocity.Get(i2) +
      0.5 * _dt * whole_all_force.Get(i2) / whole_mass.Get(i2) * _fmt2v;
    // minimum image

    // r01,r02,r12 = distance vec between atoms, with PBC

    Vec3f r01 = locator.MinDistanceVec(whole_pts.Get(i1), whole_pts.Get(i0), _box);
    Vec3f r02 = locator.MinDistanceVec(whole_pts.Get(i2), whole_pts.Get(i0), _box);
    Vec3f r12 = locator.MinDistanceVec(whole_pts.Get(i2), whole_pts.Get(i1), _box);
    //std::cout << "Before: r01 = " << r01 << std::endl;

    //domain->minimum_image(r01);

    Vec3f sv01;
    sv01 = shake_velocity[1] - shake_velocity[0];
    Vec3f sv02;                              
    sv02 = shake_velocity[2] - shake_velocity[0];
    Vec3f sv12;                              
    sv12 = shake_velocity[2] - shake_velocity[1];

    Real invmass0 = 1.0 / whole_mass.Get(i0);
    Real invmass1 = 1.0 / whole_mass.Get(i1);
    Real invmass2 = 1.0 / whole_mass.Get(i2);

    Vec3f c, l;
    Vec<Vec3f, 3> a;
    //Vec<Vec3f> a[3][3];

    // setup matrix

    a[0][0] = (invmass1 + invmass0) * vtkm::Dot(r01, r01);
    a[0][1] = (invmass0)*vtkm::Dot(r01, r02);
    a[0][2] = (-invmass1) * vtkm::Dot(r01, r12);
    a[1][0] = a[0][1];
    a[1][1] = (invmass0 + invmass2) * vtkm::Dot(r02, r02);
    a[1][2] = (invmass2)*vtkm::Dot(r02, r12);
    a[2][0] = a[0][2];
    a[2][1] = a[1][2];
    a[2][2] = (invmass2 + invmass1) * vtkm::Dot(r12, r12);

    // sestup RHS

    c[0] = -vtkm::Dot(sv01, r01);
    c[1] = -vtkm::Dot(sv02, r02);
    c[2] = -vtkm::Dot(sv12, r12);

    Vec<Vec3f, 3> ai;
    //Vec<Vec3f> ai[3][3];

    Real determ, determinv = 0.0;
    // calculate the determinant of the matrix

    determ = a[0][0] * a[1][1] * a[2][2] + a[0][1] * a[1][2] * a[2][0] +
      a[0][2] * a[1][0] * a[2][1] - a[0][0] * a[1][2] * a[2][1] - a[0][1] * a[1][0] * a[2][2] -
      a[0][2] * a[1][1] * a[2][0];

    // check if matrix is actually invertible

    if (vtkm::Abs(determ) < 0.0001)
      printf(" Error: Rattle determinant = 0.0 ");

    // calculate the inverse 3x3 matrix: A^(-1) = (ai_jk)

    determinv = 1.0 / determ;
    ai[0][0] = determinv * (a[1][1] * a[2][2] - a[1][2] * a[2][1]);
    ai[0][1] = -determinv * (a[0][1] * a[2][2] - a[0][2] * a[2][1]);
    ai[0][2] = determinv * (a[0][1] * a[1][2] - a[0][2] * a[1][1]);
    ai[1][0] = -determinv * (a[1][0] * a[2][2] - a[1][2] * a[2][0]);
    ai[1][1] = determinv * (a[0][0] * a[2][2] - a[0][2] * a[2][0]);
    ai[1][2] = -determinv * (a[0][0] * a[1][2] - a[0][2] * a[1][0]);
    ai[2][0] = determinv * (a[1][0] * a[2][1] - a[1][1] * a[2][0]);
    ai[2][1] = -determinv * (a[0][0] * a[2][1] - a[0][1] * a[2][0]);
    ai[2][2] = determinv * (a[0][0] * a[1][1] - a[0][1] * a[1][0]);

    // calculate the solution:  (l01, l02, l12)^T = A^(-1) * c

    for (int i = 0; i < 3; i++)
    {
      l[i] = 0;
      for (int j = 0; j < 3; j++)
        l[i] += ai[i][j] * c[j];
    }

    Vec<Vec3f, 3> velocity_constraint;
    for (int k = 0; k < 3; k++)
    {
      velocity_constraint[0][k] = -invmass0 * (l[0] * r01[k] + l[1] * r02[k]);
    }

    for (int k = 0; k < 3; k++)
    {
      velocity_constraint[1][k] = -invmass1 * (-l[0] * r01[k] + l[2] * r12[k]);
    }

    for (int k = 0; k < 3; k++)
    {
      velocity_constraint[2][k] = -invmass2 * (-l[1] * r02[k] - l[2] * r12[k]);
    }

    whole_velocity.Set(i0, shake_velocity[0] + velocity_constraint[0]);
    whole_velocity.Set(i1, shake_velocity[1] + velocity_constraint[1]);
    whole_velocity.Set(i2, shake_velocity[2] + velocity_constraint[2]);

    //whole_velocity.Set(i0, whole_velocity.Get(i0) + velocity_constraint[0]);
    //whole_velocity.Set(i1, whole_velocity.Get(i1) + velocity_constraint[1]);
    //whole_velocity.Set(i2, whole_velocity.Get(i2) + velocity_constraint[2]);

    Real dv1 = vtkm::Dot(r01, (whole_velocity.Get(i1) - whole_velocity.Get(i0)));
    Real dv2 = vtkm::Dot(r02, (whole_velocity.Get(i2) - whole_velocity.Get(i0)));
    Real dv12 = vtkm::Dot(r12, (whole_velocity.Get(i2) - whole_velocity.Get(i1)));

    //printf("dv1 = %f, dv2 = %f, dv12 = %f\n",dv1,dv2,dv12);
    //printf("\n");
  }

  Vec3f _box;
  Real _dt;
  Vec<Vec2f, 3> _range;
  Real _fmt2v;
};

// "ConstraintShakeForceStepWaterBondAngleWorklet"未使用
struct ConstraintShakeForceStepWaterBondAngleWorklet : vtkm::worklet::WorkletMapField
{
  ConstraintShakeForceStepWaterBondAngleWorklet(const Real& vlength,
                                                const Real& dt,
                                                const Real fmt2v,
                                                const Vec<Vec2f, 3>& range)
    : _vlength(vlength)
    , _dt(dt)
    , _fmt2v(fmt2v)
    , _range(range)
  {
  }
  using ControlSignature = void(const FieldIn anglelist,
                                const WholeArrayIn whole_pts,
                                const WholeArrayIn whole_velocity,
                                WholeArrayInOut whole_all_force,
                                const WholeArrayIn whole_mass,
                                FieldOut test_angle);
  using ExecutionSignature = void(_1, _2, _3, _4, _5, _6);

  template<typename AngleListType,
           typename PositionType,
           typename VelocityType,
           typename AllForceType,
           typename MassType>
  VTKM_EXEC void operator()(const AngleListType& anglelist,
                            const PositionType& whole_pts,
                            const VelocityType& whole_velocity,
                            AllForceType& whole_all_force,
                            const MassType& whole_mass,
                            Vec2f& test_angle) const
  {

    IdComponent i0 = anglelist[1];
    IdComponent i1 = anglelist[0];
    IdComponent i2 = anglelist[2];

    vtkm::Vec<vtkm::Vec3f, 3> shake_position;
    shake_position[0] = whole_pts.Get(i0) + whole_velocity.Get(i0) * _dt +
      0.5 * _dt * _dt * whole_all_force.Get(i0) / whole_mass.Get(i0) * _fmt2v;
    shake_position[1] = whole_pts.Get(i1) + whole_velocity.Get(i1) * _dt +
      0.5 * _dt * _dt * whole_all_force.Get(i1) / whole_mass.Get(i1) * _fmt2v;
    shake_position[2] = whole_pts.Get(i2) + whole_velocity.Get(i2) * _dt +
      0.5 * _dt * _dt * whole_all_force.Get(i2) / whole_mass.Get(i2) * _fmt2v;

    Real bond1 = 1.0;
    Real bond2 = 1.0;
    Real bond12 = vtkm::Sqrt(bond1 * bond1 + bond2 * bond2 -
                             2.0 * bond1 * bond2 * vtkm::Cos((109.4700 / 180.0) * vtkm::Pif()));
    //Real bond12 = 1.697056274848; // sqrt(2) * 1.2

    // minimum image
    // r01,r02,r12 = distance vec between atoms, with PBC
    Vec3f r01;
    r01 = whole_pts.Get(i0) - whole_pts.Get(i1);
    //std::cout << "Before: r01 = " << r01 << std::endl;
    //domain->minimum_image(r01);
    while (vtkm::Abs(r01[0]) > _vlength / 2.0)
    {
      if (r01[0] < 0)
        r01[0] = r01[0] + _vlength;
      else
        r01[0] = r01[0] - _vlength;
    }
    while (vtkm::Abs(r01[1]) > _vlength / 2.0)
    {
      if (r01[1] < 0)
        r01[1] = r01[1] + _vlength;
      else
        r01[1] = r01[1] - _vlength;
    }
    while (vtkm::Abs(r01[2]) > _vlength / 2.0)
    {
      if (r01[2] < 0)
        r01[2] = r01[2] + _vlength;
      else
        r01[2] = r01[2] - _vlength;
    }
    Vec3f r02;
    r02 = whole_pts.Get(i0) - whole_pts.Get(i2);
    //domain->minimum_image(r02);
    while (vtkm::Abs(r02[0]) > _vlength / 2.0)
    {
      if (r02[0] < 0)
        r02[0] = r02[0] + _vlength;
      else
        r02[0] = r02[0] - _vlength;
    }
    while (vtkm::Abs(r02[1]) > _vlength / 2.0)
    {
      if (r02[1] < 0)
        r02[1] = r02[1] + _vlength;
      else
        r02[1] = r02[1] - _vlength;
    }
    while (vtkm::Abs(r02[2]) > _vlength / 2.0)
    {
      if (r02[2] < 0)
        r02[2] = r02[2] + _vlength;
      else
        r02[2] = r02[2] - _vlength;
    }
    Vec3f r12;
    r12 = whole_pts.Get(i1) - whole_pts.Get(i2);
    //domain->minimum_image(r12);
    while (vtkm::Abs(r12[0]) > _vlength / 2.0)
    {
      if (r12[0] < 0)
        r12[0] = r12[0] + _vlength;
      else
        r12[0] = r12[0] - _vlength;
    }
    while (vtkm::Abs(r12[1]) > _vlength / 2.0)
    {
      if (r12[1] < 0)
        r12[1] = r12[1] + _vlength;
      else
        r12[1] = r12[1] - _vlength;
    }
    while (vtkm::Abs(r12[2]) > _vlength / 2.0)
    {
      if (r12[2] < 0)
        r12[2] = r12[2] + _vlength;
      else
        r12[2] = r12[2] - _vlength;
    }
    //std::cout << "After: r01 = " << r01 << " r02 = " << r02 << " r12 = " << r12 << std::endl;
    // s01,s02,s12 = distance vec after unconstrained update, with PBC
    Vec3f s01;
    s01 = shake_position[0] - shake_position[1];
                                           
    Vec3f s02;                             
    s02 = shake_position[0] - shake_position[2];
                                            
    Vec3f s12;                             
    s12 = shake_position[1] - shake_position[2];
    //std::cout << "Before: s01 = " << s01 << " s02 = " << s02 << " s12 = " << s12 << std::endl;
    //domain->minimum_image_once(s12);
    //std::cout << "s01 = " << s01 << std::endl;
    if (vtkm::Abs(s01[0]) > _vlength / 2.0)
    {
      if (s01[0] < 0)
        s01[0] = s01[0] + _vlength;
      else
        s01[0] = s01[0] - _vlength;
    }
    if (vtkm::Abs(s01[1]) > _vlength / 2.0)
    {
      if (s01[1] < 0)
        s01[1] = s01[1] + _vlength;
      else
        s01[1] = s01[1] - _vlength;
    }
    if (vtkm::Abs(s01[2]) > _vlength / 2.0)
    {
      if (s01[2] < 0)
        s01[2] = s01[2] + _vlength;
      else
        s01[2] = s01[2] - _vlength;
    }

    if (vtkm::Abs(s02[0]) > _vlength / 2.0)
    {
      if (s02[0] < 0)
        s02[0] = s02[0] + _vlength;
      else
        s02[0] = s02[0] - _vlength;
    }
    if (vtkm::Abs(s02[1]) > _vlength / 2.0)
    {
      if (s02[1] < 0)
        s02[1] = s02[1] + _vlength;
      else
        s02[1] = s02[1] - _vlength;
    }
    if (vtkm::Abs(s02[2]) > _vlength / 2.0)
    {
      if (s02[2] < 0)
        s02[2] = s02[2] + _vlength;
      else
        s02[2] = s02[2] - _vlength;
    }

    if (vtkm::Abs(s12[0]) > _vlength / 2.0)
    {
      if (s12[0] < 0)
        s12[0] = s12[0] + _vlength;
      else
        s12[0] = s12[0] - _vlength;
    }
    if (vtkm::Abs(s12[1]) > _vlength / 2.0)
    {
      if (s12[1] < 0)
        s12[1] = s12[1] + _vlength;
      else
        s12[1] = s12[1] - _vlength;
    }
    if (vtkm::Abs(s12[2]) > _vlength / 2.0)
    {
      if (s12[2] < 0)
        s12[2] = s12[2] + _vlength;
      else
        s12[2] = s12[2] - _vlength;
    }
    //std::cout << "After: s01 = " << s01 << " s02 = " << s02 << " s12 = " << s12 << std::endl;

    // scalar distances between atoms
    Real r01sq = vtkm::Dot(r01, r01);
    Real r02sq = vtkm::Dot(r02, r02);
    Real r12sq = vtkm::Dot(r12, r12);
    Real s01sq = vtkm::Dot(s01, s01);
    Real s02sq = vtkm::Dot(s02, s02);
    Real s12sq = vtkm::Dot(s12, s12);

    // matrix coeffs and rhs for lamda equations
    Real invmass0 = 1.0 / whole_mass.Get(i0);
    Real invmass1 = 1.0 / whole_mass.Get(i1);
    Real invmass2 = 1.0 / whole_mass.Get(i2);
    Real a11 = 2.0 * (invmass0 + invmass1) * vtkm::Dot(s01, r01);
    Real a12 = 2.0 * invmass0 * vtkm::Dot(s01, r02);
    Real a13 = -2.0 * invmass1 * vtkm::Dot(s01, r12);
    Real a21 = 2.0 * invmass0 * vtkm::Dot(s02, r01);
    Real a22 = 2.0 * (invmass0 + invmass2) * vtkm::Dot(s02, r02);
    Real a23 = 2.0 * invmass2 * vtkm::Dot(s02, r12);
    Real a31 = -2.0 * invmass1 * vtkm::Dot(s12, r01);
    Real a32 = 2.0 * invmass2 * vtkm::Dot(s12, r02);
    Real a33 = 2.0 * (invmass1 + invmass2) * (vtkm::Dot(s12, r12));

    // inverse of matrix
    Real determ = a11 * a22 * a33 + a12 * a23 * a31 + a13 * a21 * a32 - a11 * a23 * a32 -
      a12 * a21 * a33 - a13 * a22 * a31;
    if (vtkm::Abs(determ) < 0.0001)
    {
      printf("Shake determinant = 0.0");
    }

    Real determinv = 1.0 / determ;

    Real a11inv = determinv * (a22 * a33 - a23 * a32);
    Real a12inv = -determinv * (a12 * a33 - a13 * a32);
    Real a13inv = determinv * (a12 * a23 - a13 * a22);
    Real a21inv = -determinv * (a21 * a33 - a23 * a31);
    Real a22inv = determinv * (a11 * a33 - a13 * a31);
    Real a23inv = -determinv * (a11 * a23 - a13 * a21);
    Real a31inv = determinv * (a21 * a32 - a22 * a31);
    Real a32inv = -determinv * (a11 * a32 - a12 * a31);
    Real a33inv = determinv * (a11 * a22 - a12 * a21);

    // quadratic correction coeffs
    Real r0102 = vtkm::Dot(r01, r02);
    Real r0112 = vtkm::Dot(r01, r12);
    Real r0212 = vtkm::Dot(r02, r12);

    Real quad1_0101 = (invmass0 + invmass1) * (invmass0 + invmass1) * r01sq;
    Real quad1_0202 = invmass0 * invmass0 * r02sq;
    Real quad1_1212 = invmass1 * invmass1 * r12sq;
    Real quad1_0102 = 2.0 * (invmass0 + invmass1) * invmass0 * r0102;
    Real quad1_0112 = -2.0 * (invmass0 + invmass1) * invmass1 * r0112;
    Real quad1_0212 = -2.0 * invmass0 * invmass1 * r0212;
    Real quad2_0101 = invmass0 * invmass0 * r01sq;
    Real quad2_0202 = (invmass0 + invmass2) * (invmass0 + invmass2) * r02sq;
    Real quad2_1212 = invmass2 * invmass2 * r12sq;
    Real quad2_0102 = 2.0 * (invmass0 + invmass2) * invmass0 * r0102;
    Real quad2_0112 = 2.0 * invmass0 * invmass2 * r0112;
    Real quad2_0212 = 2.0 * (invmass0 + invmass2) * invmass2 * r0212;
    Real quad3_0101 = invmass1 * invmass1 * r01sq;
    Real quad3_0202 = invmass2 * invmass2 * r02sq;
    Real quad3_1212 = (invmass1 + invmass2) * (invmass1 + invmass2) * r12sq;
    Real quad3_0102 = -2.0 * invmass1 * invmass2 * r0102;
    Real quad3_0112 = -2.0 * (invmass1 + invmass2) * invmass1 * r0112;
    Real quad3_0212 = 2.0 * (invmass1 + invmass2) * invmass2 * r0212;

    // iterate until converged

    Real tolerance = 0.0001;    // original 0.001
    IdComponent max_iter = 100; // original: 100

    Real lamda01 = 0.0;
    Real lamda02 = 0.0;
    Real lamda12 = 0.0;
    IdComponent niter = 0;
    IdComponent done = 0;
    IdComponent flag_overflow = 0;
    Real quad1, quad2, quad3, b1, b2, b3, lamda01_new, lamda02_new, lamda12_new;

    while (!done && niter < max_iter)
    {

      quad1 = quad1_0101 * lamda01 * lamda01 + quad1_0202 * lamda02 * lamda02 +
        quad1_1212 * lamda12 * lamda12 + quad1_0102 * lamda01 * lamda02 +
        quad1_0112 * lamda01 * lamda12 + quad1_0212 * lamda02 * lamda12;

      quad2 = quad2_0101 * lamda01 * lamda01 + quad2_0202 * lamda02 * lamda02 +
        quad2_1212 * lamda12 * lamda12 + quad2_0102 * lamda01 * lamda02 +
        quad2_0112 * lamda01 * lamda12 + quad2_0212 * lamda02 * lamda12;

      quad3 = quad3_0101 * lamda01 * lamda01 + quad3_0202 * lamda02 * lamda02 +
        quad3_1212 * lamda12 * lamda12 + quad3_0102 * lamda01 * lamda02 +
        quad3_0112 * lamda01 * lamda12 + quad3_0212 * lamda02 * lamda12;

      b1 = bond1 * bond1 - s01sq - quad1;
      b2 = bond2 * bond2 - s02sq - quad2;
      b3 = bond12 * bond12 - s12sq - quad3;

      //b1 = bond1*bond1 - s01sq;
      //b2 = bond2*bond2 - s02sq;
      //b3 = bond12*bond12 - s12sq;

      lamda01_new = a11inv * b1 + a12inv * b2 + a13inv * b3;
      lamda02_new = a21inv * b1 + a22inv * b2 + a23inv * b3;
      lamda12_new = a31inv * b1 + a32inv * b2 + a33inv * b3;

      //std::cout << "vtkm::Abs(lamda01_new-lamda01) = " << vtkm::Abs(lamda01_new-lamda01) << std::endl;
      done = 1;
      if (vtkm::Abs(lamda01_new - lamda01) > tolerance)
        done = 0;
      if (vtkm::Abs(lamda02_new - lamda02) > tolerance)
        done = 0;
      if (vtkm::Abs(lamda12_new - lamda12) > tolerance)
        done = 0;

      //std::cout << "In process: lambda01 = " << lamda01 << std::endl;
      lamda01 = lamda01_new;
      lamda02 = lamda02_new;
      lamda12 = lamda12_new;

      // stop iterations before we have a floating point overflow
      // max double is < 1.0e308, so 1e150 is a reasonable cutoff

      if (vtkm::IsNan(lamda01) || vtkm::IsNan(lamda02) || vtkm::IsNan(lamda12) ||
          vtkm::Abs(lamda01) > 1e20 || vtkm::Abs(lamda02) > 1e20 || vtkm::Abs(lamda12) > 1e20)
      {
        done = 1;
        flag_overflow = 1;
      }

      //std::cout << "niter = " << niter << std::endl;
      niter++;
    }

    Real dtfsq = 0.5 * _dt * _dt * _fmt2v;
    lamda01 = lamda01 / dtfsq;
    lamda02 = lamda02 / dtfsq;
    lamda12 = lamda12 / dtfsq;

    Vec3f shake_force_i0;
    Vec3f shake_force_i1;
    Vec3f shake_force_i2;

    if (!flag_overflow)
    {
      shake_force_i0 = { lamda01 * r01[0] + lamda02 * r02[0],
                         lamda01 * r01[1] + lamda02 * r02[1],
                         lamda01 * r01[2] + lamda02 * r02[2] };
      shake_force_i1 = { lamda01 * r01[0] - lamda12 * r12[0],
                         lamda01 * r01[1] - lamda12 * r12[1],
                         lamda01 * r01[2] - lamda12 * r12[2] };
      shake_force_i2 = { lamda02 * r02[0] + lamda12 * r12[0],
                         lamda02 * r02[1] + lamda12 * r12[1],
                         lamda02 * r02[2] + lamda12 * r12[2] };
    }
    else
    {
      shake_force_i0 = { 0.0, 0.0, 0.0 };
      shake_force_i1 = { 0.0, 0.0, 0.0 };
      shake_force_i2 = { 0.0, 0.0, 0.0 };
      printf("flag_overflow = %d\n", flag_overflow);
    }


    whole_all_force.Set(i0, whole_all_force.Get(i0) + shake_force_i0);
    whole_all_force.Set(i1, whole_all_force.Get(i1) - shake_force_i1);
    whole_all_force.Set(i2, whole_all_force.Get(i2) - shake_force_i2);

    //printf("i0: force[0]=%f, force[1]=%f, force[2]=%f\n",shake_force_i0[0],shake_force_i0[1],shake_force_i0[2] );
    //printf("i1: force[0]=%f, force[1]=%f, force[2]=%f\n",shake_force_i1[0],shake_force_i1[1],shake_force_i1[2] );
    //printf("i2: force[0]=%f, force[1]=%f, force[2]=%f\n",shake_force_i2[0],shake_force_i2[1],shake_force_i2[2] );
    //printf("\n");

    Vec3f r_ij = whole_pts.Get(anglelist[0]) - whole_pts.Get(anglelist[1]);
    if (vtkm::Abs(r_ij[0]) > _vlength / 2.0)
    {
      if (r_ij[0] < 0)
        r_ij[0] += _vlength;
      else
        r_ij[0] -= _vlength;
    }
    if (vtkm::Abs(r_ij[1]) > _vlength / 2.0)
    {
      if (r_ij[1] < 0)
        r_ij[1] += _vlength;
      else
        r_ij[1] -= _vlength;
    }
    if (vtkm::Abs(r_ij[2]) > _vlength / 2.0)
    {
      if (r_ij[2] < 0)
        r_ij[2] += _vlength;
      else
        r_ij[2] -= _vlength;
    }
    Real disij = vtkm::Magnitude(r_ij);

    Vec3f r_kj = whole_pts.Get(anglelist[2]) - whole_pts.Get(anglelist[1]);
    if (vtkm::Abs(r_kj[0]) > _vlength / 2.0)
    {
      if (r_kj[0] < 0)
        r_kj[0] += _vlength;
      else
        r_kj[0] -= _vlength;
    }
    if (vtkm::Abs(r_kj[1]) > _vlength / 2.0)
    {
      if (r_kj[1] < 0)
        r_kj[1] += _vlength;
      else
        r_kj[1] -= _vlength;
    }
    if (vtkm::Abs(r_kj[2]) > _vlength / 2.0)
    {
      if (r_kj[2] < 0)
        r_kj[2] += _vlength;
      else
        r_kj[2] -= _vlength;
    }
    Real diskj = vtkm::Magnitude(r_kj);
    Real cosangle = vtkm::Dot(r_ij, r_kj);
    cosangle /= disij * diskj;
    if (cosangle > 1.0)
      cosangle = 1.0;
    if (cosangle < -1.0)
      cosangle = -1.0;
    test_angle[1] = acos(cosangle) / vtkm::Pif() * 180;
  }

  Real _vlength;
  Real _dt;
  Vec<Vec2f, 3> _range;
  Real _fmt2v;
};

struct NewConstraintAWaterBondAngleWorklet : vtkm::worklet::WorkletMapField
{
  NewConstraintAWaterBondAngleWorklet(const Vec3f& box,
                                      const Real& dt,
                                      const Real fmt2v,
                                      const Vec<Vec2f, 3>& range)
    : _box(box)
    , _dt(dt)
    , _fmt2v(fmt2v)
    , _range(range)
  {
  }
  using ControlSignature = void(const FieldIn anglelist,
                                WholeArrayInOut whole_pts,
                                WholeArrayInOut whole_velocity,
                                const WholeArrayIn whole_all_force,

                                const WholeArrayIn whole_mass , 
                                WholeArrayInOut whole_pts_flag,
                                ExecObject locator);
  using ExecutionSignature = void(_1, _2, _3, _4, _5 ,_6,_7);

  template<typename AngleListType,
           typename PositionType,
           typename VelocityType,
           typename AllForceType,
           typename MassType,
           typename PtsFlagType>
  VTKM_EXEC void operator()(const AngleListType& anglelist,
                            PositionType& whole_pts,
                            VelocityType& whole_velocity,
                            const AllForceType& whole_all_force,
                            const MassType& whole_mass,
                            PtsFlagType& whole_pts_flag,
                            const ExecPointLocator& locator) const
  {
    IdComponent i0 = anglelist[0]; //H
    IdComponent i1 = anglelist[1]; //O
    IdComponent i2 = anglelist[2]; //H

    vtkm::Vec<Vec3f, 3> shake_position;
    shake_position[0] = whole_pts.Get(i0) + whole_velocity.Get(i0) * _dt +
      0.5 * _dt * _dt * whole_all_force.Get(i0) / whole_mass.Get(i0) * _fmt2v;
    shake_position[1] = whole_pts.Get(i1) + whole_velocity.Get(i1) * _dt +
      0.5 * _dt * _dt * whole_all_force.Get(i1) / whole_mass.Get(i1) * _fmt2v;
    shake_position[2] = whole_pts.Get(i2) + whole_velocity.Get(i2) * _dt +
      0.5 * _dt * _dt * whole_all_force.Get(i2) / whole_mass.Get(i2) * _fmt2v;


    Real bond1 = 1.0;
    Real bond2 = 1.0;
    Real bond12 = vtkm::Sqrt(bond1 * bond1 + bond2 * bond2 -
                             2.0 * bond1 * bond2 * vtkm::Cos((109.4700 / 180.0) * vtkm::Pif()));

    // minimum image
    Vec3f r01 = locator.MinDistanceVec(whole_pts.Get(i1), whole_pts.Get(i0), _box);
    Vec3f r12 = locator.MinDistanceVec(whole_pts.Get(i2), whole_pts.Get(i1), _box);
    Vec3f r20 = locator.MinDistanceVec(whole_pts.Get(i0), whole_pts.Get(i2), _box);

    // s01,s02,s12 = distance vec after unconstrained update, with PBC

    Vec3f s10 = locator.MinDistanceVec(shake_position[0], shake_position[1], _box);
    Vec3f s21 = locator.MinDistanceVec(shake_position[1], shake_position[2], _box);
    Vec3f s02 = locator.MinDistanceVec(shake_position[2], shake_position[0], _box);

    // scalar distances between atoms

    Real r01sq = vtkm::Dot(r01, r01);
    Real r02sq = vtkm::Dot(r20, r20);
    Real r12sq = vtkm::Dot(r12, r12);
    Real s01sq = vtkm::Dot(s10, s10);
    Real s02sq = vtkm::Dot(s02, s02);
    Real s12sq = vtkm::Dot(s21, s21);

    // matrix coeffs and rhs for lamda equations

    Real invmass0 = 1.0 / whole_mass.Get(i0);
    Real invmass1 = 1.0 / whole_mass.Get(i1);
    Real invmass2 = 1.0 / whole_mass.Get(i2);
    Real a11 = 2.0 * (invmass0 + invmass1) * vtkm::Dot(s10, r01);
    Real a12 = -2.0 * invmass1 * vtkm::Dot(s10, r12);
    Real a13 = -2.0 * invmass0 * vtkm::Dot(s10, r20);
    Real a21 = -2.0 * invmass1 * vtkm::Dot(s21, r01);
    Real a22 = 2.0 * (invmass1 + invmass2) * vtkm::Dot(s21, r12);
    Real a23 = -2.0 * invmass2 * vtkm::Dot(s21, r20);
    Real a31 = -2.0 * invmass0 * vtkm::Dot(s02, r01);
    Real a32 = -2.0 * invmass2 * vtkm::Dot(s02, r12);
    Real a33 = 2.0 * (invmass0 + invmass2) * (vtkm::Dot(s02, r20));

    // inverse of matrix

    Real determ = a11 * a22 * a33 + a12 * a23 * a31 + a13 * a21 * a32 - a11 * a23 * a32 -
      a12 * a21 * a33 - a13 * a22 * a31;
    if (vtkm::Abs(determ) < 0.0001)
    {
      printf("Shake determinant = 0.0");
    }


    Real determinv = 1.0 / determ;

    Real a11inv = determinv * (a22 * a33 - a23 * a32);
    Real a12inv = -determinv * (a12 * a33 - a13 * a32);
    Real a13inv = determinv * (a12 * a23 - a13 * a22);
    Real a21inv = -determinv * (a21 * a33 - a23 * a31);
    Real a22inv = determinv * (a11 * a33 - a13 * a31);
    Real a23inv = -determinv * (a11 * a23 - a13 * a21);
    Real a31inv = determinv * (a21 * a32 - a22 * a31);
    Real a32inv = -determinv * (a11 * a32 - a12 * a31);
    Real a33inv = determinv * (a11 * a22 - a12 * a21);


    Real r0120 = vtkm::Dot(r01, r20);
    Real r0112 = vtkm::Dot(r01, r12);
    Real r2012 = vtkm::Dot(r20, r12);

    Real quad1_0101 = (invmass0 + invmass1) * (invmass0 + invmass1) * r01sq;
    Real quad1_1212 = invmass1 * invmass1 * r12sq;
    Real quad1_2020 = invmass0 * invmass0 * r02sq;
    Real quad1_0120 = -2.0 * (invmass0 + invmass1) * invmass0 * r0120;
    Real quad1_0112 = -2.0 * (invmass0 + invmass1) * invmass1 * r0112;
    Real quad1_2012 = 2.0 * invmass0 * invmass1 * r2012;

    Real quad2_0101 = invmass1 * invmass1 * r01sq;
    Real quad2_1212 = (invmass1 + invmass2) * (invmass1 + invmass2) * r12sq;
    Real quad2_2020 = invmass2 * invmass2 * r02sq;
    Real quad2_0120 = 2.0 * invmass1 * invmass2 * r0120;
    Real quad2_0112 = -2.0 * (invmass1 + invmass2) * invmass1 * r0112;
    Real quad2_2012 = -2.0 * (invmass1 + invmass2) * invmass2 * r2012;

    Real quad3_0101 = invmass0 * invmass0 * r01sq;
    Real quad3_1212 = invmass2 * invmass2 * r12sq;
    Real quad3_2020 = (invmass0 + invmass2) * (invmass0 + invmass2) * r02sq;
    Real quad3_0120 = -2.0 * (invmass0 + invmass2) * invmass0 * r0120;
    Real quad3_0112 = 2.0 * invmass0 * invmass2 * r0112;
    Real quad3_2012 = -2.0 * (invmass0 + invmass2) * invmass2 * r2012;

    // iterate until converged
    Real tolerance = 0.00001;     // original 0.001
    IdComponent max_iter = 5000; // original: 100

    Real lamda01 = 0.0;
    Real lamda20 = 0.0;
    Real lamda12 = 0.0;
    IdComponent niter = 0;
    IdComponent done = 0;
    IdComponent flag_overflow = 0;
    Real quad1, quad2, quad3, b1, b2, b3, lamda01_new, lamda20_new, lamda12_new;

    while (!done && niter < max_iter)
    {

      quad1 = quad1_0101 * lamda01 * lamda01 + quad1_2020 * lamda20 * lamda20 +
        quad1_1212 * lamda12 * lamda12 + quad1_0120 * lamda01 * lamda20 +
        quad1_0112 * lamda01 * lamda12 + quad1_2012 * lamda20 * lamda12;

      quad2 = quad2_0101 * lamda01 * lamda01 + quad2_2020 * lamda20 * lamda20 +
        quad2_1212 * lamda12 * lamda12 + quad2_0120 * lamda01 * lamda20 +
        quad2_0112 * lamda01 * lamda12 + quad2_2012 * lamda20 * lamda12;

      quad3 = quad3_0101 * lamda01 * lamda01 + quad3_2020 * lamda20 * lamda20 +
        quad3_1212 * lamda12 * lamda12 + quad3_0120 * lamda01 * lamda20 +
        quad3_0112 * lamda01 * lamda12 + quad3_2012 * lamda20 * lamda12;

      b1 = bond1 * bond1 - s01sq - quad1;
      b2 = bond2 * bond2 - s12sq - quad2;
      b3 = bond12 * bond12 - s02sq - quad3;

      lamda01_new = a11inv * b1 + a12inv * b2 + a13inv * b3;
      lamda12_new = a21inv * b1 + a22inv * b2 + a23inv * b3;
      lamda20_new = a31inv * b1 + a32inv * b2 + a33inv * b3;

      done = 1;
      if (vtkm::Abs(lamda01_new - lamda01) > tolerance)
        done = 0;
      if (vtkm::Abs(lamda20_new - lamda20) > tolerance)
        done = 0;
      if (vtkm::Abs(lamda12_new - lamda12) > tolerance)
        done = 0;

      lamda01 = lamda01_new;
      lamda20 = lamda20_new;
      lamda12 = lamda12_new;


      if (vtkm::IsNan(lamda01) || vtkm::IsNan(lamda20) || vtkm::IsNan(lamda12) ||
          vtkm::Abs(lamda01) > 1e20 || vtkm::Abs(lamda20) > 1e20 || vtkm::Abs(lamda12) > 1e20)
      {
        done = 1;
        flag_overflow = 1;
      }
      niter++;
    }

    Vec3f position_constraint_i0;
    Vec3f position_constraint_i1;
    Vec3f position_constraint_i2;

    position_constraint_i0 = lamda01 * r01 * invmass0 - lamda20 * r20 * invmass0;
    position_constraint_i1 = lamda12 * r12 * invmass1 - lamda01 * r01 * invmass1;
    position_constraint_i2 = lamda20 * r20 * invmass2 - lamda12 * r12 * invmass2;


    Vec3f velocity_constraint_i0 = position_constraint_i0 / _dt;
    Vec3f velocity_constraint_i1 = position_constraint_i1 / _dt;
    Vec3f velocity_constraint_i2 = position_constraint_i2 / _dt;

    vtkm::Vec<Vec3f, 3> shake_velocity;
    shake_velocity[0] = whole_velocity.Get(i0) + 0.5 * _dt * whole_all_force.Get(i0) / whole_mass.Get(i0) * _fmt2v;
    shake_velocity[1] = whole_velocity.Get(i1) + 0.5 * _dt * whole_all_force.Get(i1) / whole_mass.Get(i1) * _fmt2v;
    shake_velocity[2] = whole_velocity.Get(i2) + 0.5 * _dt * whole_all_force.Get(i2) / whole_mass.Get(i2) * _fmt2v;

    whole_velocity.Set(i0, shake_velocity[0] + velocity_constraint_i0);
    whole_velocity.Set(i1, shake_velocity[1] + velocity_constraint_i1);
    whole_velocity.Set(i2, shake_velocity[2] + velocity_constraint_i2);

    whole_pts.Set(i0, shake_position[0] + position_constraint_i0);
    whole_pts.Set(i1, shake_position[1] + position_constraint_i1);
    whole_pts.Set(i2, shake_position[2] + position_constraint_i2);

    // pbc
    vtkm::Vec<vtkm::Vec3f, 3> whole_position_pbc;
    whole_position_pbc[0] = whole_pts.Get(anglelist[0]);
    whole_position_pbc[1] = whole_pts.Get(anglelist[1]);
    whole_position_pbc[2] = whole_pts.Get(anglelist[2]);

    vtkm::Vec<Id3,3> whole_pts_flag_pbc;
    whole_pts_flag_pbc[0] = whole_pts_flag.Get(anglelist[0]);
    whole_pts_flag_pbc[1] = whole_pts_flag.Get(anglelist[1]);
    whole_pts_flag_pbc[2] = whole_pts_flag.Get(anglelist[2]);
    
    for (auto i = 0; i < 3; i++)
    {
      for (auto j = 0; j < 3; j++)
      {
        if (whole_position_pbc[i][j] < _range[j][0])
        {
          whole_position_pbc[i][j] += (_range[j][1] - _range[j][0]);
          whole_pts_flag_pbc[i][j] -= 1;
        }
        if (whole_position_pbc[i][j] > _range[j][1])
        {
          whole_position_pbc[i][j] -= (_range[j][1] - _range[j][0]);
          whole_pts_flag_pbc[i][j] += 1;
        }
      }
    }

    for (int i = 0; i < 3; i++)
    {
      whole_pts.Set(anglelist[i], whole_position_pbc[i]);
      whole_pts_flag.Set(anglelist[i], whole_pts_flag_pbc[i]);
    }
  }

  Vec3f _box;
  Real _dt;
  Vec<Vec2f, 3> _range;
  Real _fmt2v;
};

struct NewConstraintBWaterBondAngleWorklet : vtkm::worklet::WorkletMapField
{
  NewConstraintBWaterBondAngleWorklet(const Vec3f& box,
                                      const Real& dt,
                                      const Real fmt2v,
                                      const Vec<Vec2f, 3>& range)
    : _box(box)
    , _dt(dt)
    , _fmt2v(fmt2v)
    , _range(range)
  {
  }
  using ControlSignature = void(const FieldIn anglelist,
                                WholeArrayInOut whole_pts,
                                WholeArrayInOut whole_velocity,
                                const WholeArrayIn whole_all_force,
                                const WholeArrayIn whole_mass,
                                ExecObject locator);
  using ExecutionSignature = void(_1, _2, _3, _4, _5, _6);

  template<typename AngleListType,
           typename PositionType,
           typename VelocityType,
           typename AllForceType,
           typename MassType>
  VTKM_EXEC void operator()(const AngleListType& anglelist,
                            PositionType& whole_pts,
                            VelocityType& whole_velocity,
                            const AllForceType& whole_all_force,
                            const MassType& whole_mass,
                            const ExecPointLocator& locator) const
  {
    IdComponent i0 = anglelist[0]; //H
    IdComponent i1 = anglelist[1]; //O
    IdComponent i2 = anglelist[2]; //H

    vtkm::Vec<vtkm::Vec3f, 3> shake_velocity;
    shake_velocity[0] = whole_velocity.Get(i0) +
      0.5 * _dt * whole_all_force.Get(i0) / whole_mass.Get(i0) * _fmt2v;
    shake_velocity[1] = whole_velocity.Get(i1) +
      0.5 * _dt * whole_all_force.Get(i1) / whole_mass.Get(i1) * _fmt2v;
    shake_velocity[2] = whole_velocity.Get(i2) +
      0.5 * _dt * whole_all_force.Get(i2) / whole_mass.Get(i2) * _fmt2v;

    // minimum image

    Vec3f r01 = locator.MinDistanceVec(whole_pts.Get(i1), whole_pts.Get(i0), _box);
    Vec3f r12 = locator.MinDistanceVec(whole_pts.Get(i2), whole_pts.Get(i1), _box);
    Vec3f r20 = locator.MinDistanceVec(whole_pts.Get(i0), whole_pts.Get(i2), _box);

    // s01,s02,s12 = distance vec after unconstrained update, with PBC


    Vec3f sv10;
    sv10 = shake_velocity[0] - shake_velocity[1];
    Vec3f sv21;                            
    sv21 = shake_velocity[1] - shake_velocity[2];
    Vec3f sv02;                              
    sv02 = shake_velocity[2] - shake_velocity[0];


    Real invmass0 = 1.0 / whole_mass.Get(i0);
    Real invmass1 = 1.0 / whole_mass.Get(i1);
    Real invmass2 = 1.0 / whole_mass.Get(i2);

    Vec3f c, l;
    Vec<Vec3f, 3> a;
    //Vec<Vec3f> a[3][3];

    // setup matrix

    a[0][0] = (invmass1 + invmass0) * vtkm::Dot(r01, r01);
    a[0][1] = -invmass1 * vtkm::Dot(r01, r12);
    a[0][2] = (-invmass0) * vtkm::Dot(r01, r20);
    a[1][0] = a[0][1];
    a[1][1] = (invmass1 + invmass2) * vtkm::Dot(r12, r12);
    a[1][2] = -(invmass2)*vtkm::Dot(r20, r12);
    a[2][0] = a[0][2];
    a[2][1] = a[1][2];
    a[2][2] = (invmass0 + invmass2) * vtkm::Dot(r20, r20);

    // sestup RHS

    c[0] = -vtkm::Dot(sv10, r01);
    c[1] = -vtkm::Dot(sv21, r12);
    c[2] = -vtkm::Dot(sv02, r20);

    Vec<Vec3f, 3> ai;
    //Vec<Vec3f> ai[3][3];

    Real determ, determinv = 0.0;
    // calculate the determinant of the matrix

    determ = a[0][0] * a[1][1] * a[2][2] + a[0][1] * a[1][2] * a[2][0] +
      a[0][2] * a[1][0] * a[2][1] - a[0][0] * a[1][2] * a[2][1] - a[0][1] * a[1][0] * a[2][2] -
      a[0][2] * a[1][1] * a[2][0];

    // check if matrix is actually invertible

    if (vtkm::Abs(determ) < 0.0001)
      printf(" Error: Rattle determinant = 0.0 ");

    // calculate the inverse 3x3 matrix: A^(-1) = (ai_jk)

    determinv = 1.0 / determ;
    ai[0][0] = determinv * (a[1][1] * a[2][2] - a[1][2] * a[2][1]);
    ai[0][1] = -determinv * (a[0][1] * a[2][2] - a[0][2] * a[2][1]);
    ai[0][2] = determinv * (a[0][1] * a[1][2] - a[0][2] * a[1][1]);
    ai[1][0] = -determinv * (a[1][0] * a[2][2] - a[1][2] * a[2][0]);
    ai[1][1] = determinv * (a[0][0] * a[2][2] - a[0][2] * a[2][0]);
    ai[1][2] = -determinv * (a[0][0] * a[1][2] - a[0][2] * a[1][0]);
    ai[2][0] = determinv * (a[1][0] * a[2][1] - a[1][1] * a[2][0]);
    ai[2][1] = -determinv * (a[0][0] * a[2][1] - a[0][1] * a[2][0]);
    ai[2][2] = determinv * (a[0][0] * a[1][1] - a[0][1] * a[1][0]);

    // calculate the solution:  (l01, l02, l12)^T = A^(-1) * c

    for (int i = 0; i < 3; i++)
    {
      l[i] = 0;
      for (int j = 0; j < 3; j++)
        l[i] += ai[i][j] * c[j];
    }

    // [l01,l02,l12]^T = [lamda12,lamda23,lamda31]^T

    Vec3f velocity_constraint_i0 = l[0] * r01 * invmass0 - l[2] * r20 * invmass0;
    Vec3f velocity_constraint_i1 = l[1] * r12 * invmass1 - l[0] * r01 * invmass1;
    Vec3f velocity_constraint_i2 = l[2] * r20 * invmass2 - l[1] * r12 * invmass2;

    whole_velocity.Set(i0, shake_velocity[0] + velocity_constraint_i0);
    whole_velocity.Set(i1, shake_velocity[1] + velocity_constraint_i1);
    whole_velocity.Set(i2, shake_velocity[2] + velocity_constraint_i2);

    //whole_velocity.Set(i0, whole_velocity.Get(i0) + velocity_constraint[0]);
    //whole_velocity.Set(i1, whole_velocity.Get(i1) + velocity_constraint[1]);
    //whole_velocity.Set(i2, whole_velocity.Get(i2) + velocity_constraint[2]);

    Vec3f a00 = (invmass1 + invmass0) * r01;
    Vec3f a01 = -invmass1 * r12;
    Vec3f a02 = (-invmass0) * r20;
    Vec3f a10 = -invmass1 * r01;
    Vec3f a11 = (invmass1 + invmass2) * r12;
    Vec3f a12 = -(invmass2)*r20;
    Vec3f a20 = (-invmass0) * r01;
    Vec3f a21 = -(invmass2)*r12;
    Vec3f a22 = (invmass0 + invmass2) * r20;

    Vec3f cal_v12t = sv10 + (a00 * l[0] + a01 * l[1] + a02 * l[2]);
    Vec3f cal_v23t = sv21 + (a10 * l[0] + a11 * l[1] + a12 * l[2]);
    Vec3f cal_v31t = sv02 + (a20 * l[0] + a21 * l[1] + a22 * l[2]);

    Real cal_dv1 = vtkm::Dot(r01, cal_v12t);
    Real cal_dv2 = vtkm::Dot(r20, cal_v31t);
    Real cal_dv12 = vtkm::Dot(r12, cal_v23t);

    Real dv1 = vtkm::Dot(r01, (whole_velocity.Get(i1) - whole_velocity.Get(i0)));
    Real dv2 = vtkm::Dot(r20, (whole_velocity.Get(i0) - whole_velocity.Get(i2)));
    Real dv12 = vtkm::Dot(r12, (whole_velocity.Get(i2) - whole_velocity.Get(i1)));

    if (vtkm::Abs(dv1) > 0.1 || vtkm::Abs(dv2) > 0.1 || vtkm::Abs(dv12) > 0.1)
    {
      printf("i0 = %d, i1 = %d, i2 = %d\n", i0, i1, i2);
      printf("dv1 = %f, dv2 = %f, dv12 = %f\n", dv1, dv2, dv12);
      printf("cal_dv1 = %f, cal_dv2 = %f, cal_dv12 = %f\n", cal_dv1, cal_dv2, cal_dv12);
      //printf("position_i0 = [%f,%f,%f], position_i1 = [%f,%f,%f], position_i2 = [%f,%f,%f]\n", whole_pts.Get(i0)[0],whole_pts.Get(i0)[1],whole_pts.Get(i0)[2],
      //whole_pts.Get(i1)[0],whole_pts.Get(i1)[1],whole_pts.Get(i1)[2],whole_pts.Get(i2)[0],whole_pts.Get(i2)[1],whole_pts.Get(i2)[2]);
      printf("velocity_i0 = [%f,%f,%f], velocity_i1 = [%f,%f,%f], velocity_i2 = [%f,%f,%f]\n",
             whole_velocity.Get(i0)[0],
             whole_velocity.Get(i0)[1],
             whole_velocity.Get(i0)[2],
             whole_velocity.Get(i1)[0],
             whole_velocity.Get(i1)[1],
             whole_velocity.Get(i1)[2],
             whole_velocity.Get(i2)[0],
             whole_velocity.Get(i2)[1],
             whole_velocity.Get(i2)[2]);
      printf("\n");
    }
  }

  Vec3f _box;
  Real _dt;
  Vec<Vec2f, 3> _range;
  Real _fmt2v;
};

struct GetPositionByTypeWorklet : vtkm::worklet::WorkletMapField
{
  using ControlSignature = void(FieldIn atom_id, WholeArrayIn whole_pts, FieldOut pts_by_type);
  using ExecutionSignature = void(_1, _2, _3);

  template<typename WholePtsType, typename PtsType, typename IdType>
  VTKM_EXEC void operator()(const IdType& atom_id,
                            const WholePtsType& whole_pts,
                            PtsType& pts_by_type) const
  {
    pts_by_type = whole_pts.Get(atom_id);
  }
};

}