﻿//==================================================================================
//  RBMD 2.2.0 is developed for random batch molecular dynamics calculation.
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
#include "locator/ContPointLocator.h"
#include "forceFunction/ContForceFunction.h"
#include "topology/ContTopology.h"
#include <vtkm/cont/ArrayHandleGroupVecVariable.h>
#include "staticTable/ContStaticTable.h"

using GroupIdIdType = vtkm::cont::ArrayHandleGroupVecVariable<vtkm::cont::ArrayHandle<vtkm::Id>,
                                                              vtkm::cont::ArrayHandle<vtkm::Id>>;

using GroupRealIdType = vtkm::cont::ArrayHandleGroupVecVariable<vtkm::cont::ArrayHandle<Real>,
                                                                vtkm::cont::ArrayHandle<vtkm::Id>>;

namespace RunWorklet 
{
    void ComputeNeighbours(const Real& cut_off,
                           const Vec3f& box,
                           const vtkm::cont::ArrayHandle<vtkm::Id>& atoms_id,
                           const ContPointLocator& locator,
                           GroupVecType& id_verletlist_group,
                           vtkm::cont::ArrayHandle<vtkm::Id>& num_verletlist,
                           CoordOffsetType& offset_verletlist_group);

    void ComputeRBLNeighboursOnce(const Id& rs_num,
                                  const Id& pice_num,
                                  const vtkm::cont::ArrayHandle<vtkm::Id>& atoms_id,
                                  const ContPointLocator& locator,
                                  GroupVecType& id_verletlist_group,
                                  GroupNumType& num_verletlist,
                                  CoordOffsetType& offset_verletlist_group);

    void NearForceVerletERF(const Real& cut_off,
                            const vtkm::cont::ArrayHandle<vtkm::Id>& atoms_id,
                            const ContPointLocator& locator,
                            const ContTopology& topology,
                            const ContForceFunction& force_function,
                            const ContStaticTable& static_table,
                            const GroupVecType& Group_j,
                            const vtkm::cont::ArrayHandle<vtkm::Id>& num_j,
                            const CoordOffsetType& coord_offset_j,
                            vtkm::cont::ArrayHandle<Vec3f>& LJforce);

   void NearForceVerlet(const Real& cut_off,
                         const Vec3f& box,
                         const vtkm::cont::ArrayHandle<vtkm::Id>& atoms_id,
                         const ContPointLocator& locator,
                         const ContTopology& topology,
                         const ContForceFunction& force_function,
                         const GroupVecType& Group_j,
                         const vtkm::cont::ArrayHandle<vtkm::Id>& num_j,
                         const CoordOffsetType& coord_offset_j,
                         vtkm::cont::ArrayHandle<Vec3f>& nearforce);

   void NearForceVerletWeightVirial(const Real& cut_off,
       const Vec3f& box,
       const Real& qqr2e,
       const vtkm::cont::ArrayHandle<vtkm::Id>& atoms_id,
       const ContPointLocator& locator,
       const ContTopology& topology,
       const ContForceFunction& force_function,
       const GroupVecType& Group_j,
       const vtkm::cont::ArrayHandle<vtkm::Id>& num_j,
       const CoordOffsetType& coord_offset_j,
       const GroupIdIdType& group_ids,
       const GroupRealIdType& group_weights,
       vtkm::cont::ArrayHandle<Vec3f>& nearforce,
       vtkm::cont::ArrayHandle<Vec6f>& nearvirial,
       vtkm::cont::ArrayHandle<Real>& energy_lj,
       vtkm::cont::ArrayHandle<Real>& Energy_coul);

   void NearForceVerletVirial(const Real& cut_off,
       const Vec3f& box,
       const Real& qqr2e,
       const vtkm::cont::ArrayHandle<vtkm::Id>& atoms_id,
       const ContPointLocator& locator,
       const ContTopology& topology,
       const ContForceFunction& force_function,
       const GroupVecType& Group_j,
       const vtkm::cont::ArrayHandle<vtkm::Id>& num_j,
       const CoordOffsetType& coord_offset_j,
       vtkm::cont::ArrayHandle<Vec3f>& nearforce,
       vtkm::cont::ArrayHandle<Vec6f>& nearvirial,
       vtkm::cont::ArrayHandle<Real>& energy_lj,
       vtkm::cont::ArrayHandle<Real>& Energy_coul);

   void LJForceVerlet(const Real& cut_off,
                      const Vec3f& box,
                      const vtkm::cont::ArrayHandle<vtkm::Id>& atoms_id,
                      const ContPointLocator& locator,
                      const ContTopology& topology,
                      const ContForceFunction& force_function,
                      const GroupVecType& Group_j,
                      const vtkm::cont::ArrayHandle<vtkm::Id>& num_j,
                      const CoordOffsetType& coord_offset_j,
                      vtkm::cont::ArrayHandle<Vec3f>& LJforce,
                      vtkm::cont::ArrayHandle<Vec6f>& LJvirial,
                      vtkm::cont::ArrayHandle<Real>& LJenergy);

   void LJEnergyVerlet(const Real& cut_off,
       const Vec3f& box,
       const vtkm::cont::ArrayHandle<vtkm::Id>& atoms_id,
       const ContPointLocator& locator,
       const ContTopology& topology,
       const ContForceFunction& force_function,
       const GroupVecType& Group_j,
       const vtkm::cont::ArrayHandle<vtkm::Id>& num_j,
       const CoordOffsetType& coord_offset_j,
       vtkm::cont::ArrayHandle<Real>& LJPE);

   void SpecialLJCoulEnergyVerletWeightVirial(const Real& cut_off,
       const Vec3f& box,
       const Real& qqr2e,
       const vtkm::cont::ArrayHandle<vtkm::Id>& atoms_id,
       const ContPointLocator& locator,
       const ContTopology& topology,
       const ContForceFunction& force_function,
       const GroupVecType& Group_j,
       const vtkm::cont::ArrayHandle<vtkm::Id>& num_j,
       const CoordOffsetType& coord_offset_j,
       const GroupIdIdType& group_ids,
       const GroupRealIdType& group_weights,
       vtkm::cont::ArrayHandle<Real>& energy_lj,
       vtkm::cont::ArrayHandle<Real>& Energy_coul);

   void LJCoulEnergyVerletVirial(const Real& cut_off,
       const Vec3f& box,
       const Real& qqr2e,
       const vtkm::cont::ArrayHandle<vtkm::Id>& atoms_id,
       const ContPointLocator& locator,
       const ContTopology& topology,
       const ContForceFunction& force_function,
       const GroupVecType& Group_j,
       const vtkm::cont::ArrayHandle<vtkm::Id>& num_j,
       const CoordOffsetType& coord_offset_j,
       vtkm::cont::ArrayHandle<Real>& energy_lj,
       vtkm::cont::ArrayHandle<Real>& Energy_coul);

   void EAMfpVerlet(const Real& cut_off,
                    const Vec3f& box,
                    const vtkm::cont::ArrayHandle<vtkm::Id>& atoms_id,
                    const vtkm::cont::ArrayHandle<Vec7f>& rhor_spline,
                    const vtkm::cont::ArrayHandle<Vec7f>& frho_spline,
                    const ContPointLocator& locator,
                    const ContForceFunction& force_function,
                    const GroupVecType& Group_j,
                    const vtkm::cont::ArrayHandle<vtkm::Id>& num_j,
                    const CoordOffsetType& coord_offset_j,
                    vtkm::cont::ArrayHandle<Real>& EAM_fp);

   void EAMForceVerlet(const Real& cut_off,
                       const Vec3f& box,
                       const vtkm::cont::ArrayHandle<vtkm::Id>& atoms_id,
                       const vtkm::cont::ArrayHandle<Vec7f>& rhor_spline,
                       const vtkm::cont::ArrayHandle<Vec7f>& z2r_spline,
                       const vtkm::cont::ArrayHandle<Real>& EAM_fp,
                       const ContPointLocator& locator,
                       const ContForceFunction& force_function,
                       const GroupVecType& Group_j,
                       const vtkm::cont::ArrayHandle<vtkm::Id>& num_j,
                       const CoordOffsetType& coord_offset_j,
                       vtkm::cont::ArrayHandle<Vec3f>& force);

   void EAM_rho(const Real& eam_cut_off,
                const Vec3f& box,
                const vtkm::cont::ArrayHandle<vtkm::Id>& atoms_id,
                const vtkm::cont::ArrayHandle<Vec7f>& rhor_spline,
                const ContPointLocator& locator,
                const ContTopology& topology,
                const ContForceFunction& force_function,
                vtkm::cont::ArrayHandle<Real>& EAM_rho);

   void EAM_fp(const vtkm::cont::ArrayHandle<vtkm::Id>& atoms_id,
               const vtkm::cont::ArrayHandle<Real>& EAM_rho,
               const vtkm::cont::ArrayHandle<Vec7f>& frho_spline,
               const ContPointLocator& locator,
               const ContTopology& topology,
               const ContForceFunction& force_function,
               vtkm::cont::ArrayHandle<Real>& fp);

   void EAM_force(const Real& eam_cut_off,
                  const Vec3f& box,
                  const vtkm::cont::ArrayHandle<vtkm::Id>& atoms_id,
                  const vtkm::cont::ArrayHandle<Real>& fp,
                  const vtkm::cont::ArrayHandle<Vec7f>& rhor_spline,
                  const vtkm::cont::ArrayHandle<Vec7f>& z2r_spline,
                  const ContPointLocator& locator,
                  const ContTopology& topology,
                  const ContForceFunction& force_function,
                  vtkm::cont::ArrayHandle<Vec3f>& eam_force);
    
     void NearForceRBLERF(const Id& rs_num,
                        const Id& pice_num,
                        const Real& qqr2e,
                        const Vec3f& box,
                        const vtkm::cont::ArrayHandle<vtkm::Id>& atoms_id,
                        const ContPointLocator& locator,
                        const ContTopology& topology,
                        const ContForceFunction& force_function,
                        const ContStaticTable& static_table,
                        const GroupVecType& id_verletlist_group,
                        const GroupNumType& num_verletlist,
                        const CoordOffsetType& offset_verletlist_group,
                        vtkm::cont::ArrayHandle<Vec3f>& corr_force,
                         vtkm::cont::ArrayHandle<Vec6f>& corr_virial);

    void NearForceRBLERFSpecialBonds(const Id& rs_num,
                                      const Id& pice_num,
                                      const Real& qqr2e,
                                      const Vec3f& box,
                                      const vtkm::cont::ArrayHandle<vtkm::Id>& atoms_id,
                                      const ContPointLocator& locator,
                                      const ContTopology& topology,
                                      const ContForceFunction& force_function,
                                      const ContStaticTable& static_table,
                                      const GroupVecType& id_verletlist_group,
                                      const GroupNumType& num_verletlist,
                                      const CoordOffsetType& offset_verletlist_group,
                                      const GroupIdIdType& group_ids,
                                      const GroupRealIdType& group_weights,
                                      vtkm::cont::ArrayHandle<Vec3f>& corr_force,
                                     vtkm::cont::ArrayHandle<Vec6f>& corr_virial);

    void NearForceRBLERFSpecialBonds_fix(const Id& rs_num,
        const Id& pice_num,
        const Real& qqr2e,
        const Vec3f& box,
        const vtkm::cont::ArrayHandle<vtkm::Id>& atoms_id,
        const ContPointLocator& locator,
        const ContTopology& topology,
        const ContForceFunction& force_function,
        const ContStaticTable& static_table,
        const GroupVecType& id_verletlist_group,
        const GroupNumType& num_verletlist,
        const CoordOffsetType& offset_verletlist_group,
        const GroupIdIdType& group_ids,
        const GroupRealIdType& group_weights,
        vtkm::cont::ArrayHandle<Vec3f>& corr_force);

    void NearForceRBL(const Id& rs_num,
                      const Id& pice_num,
                      const Real& qqr2e,
                      const vtkm::cont::ArrayHandle<vtkm::Id>& atoms_id,
                      const ContPointLocator& locator,
                      const ContTopology& topology,
                      const ContForceFunction& force_function,
                      const GroupVecType& id_verletlist_group,
                      const GroupNumType& num_verletlist,
                      const CoordOffsetType& offset_verletlist_group,
                      vtkm::cont::ArrayHandle<Vec3f>& corr_force);

    void LJForceRBL(const Id& rs_num,
                    const Id& pice_num,
                    const Vec3f& box,
                    const vtkm::cont::ArrayHandle<vtkm::Id>& atoms_id,
                    const ContPointLocator& locator,
                    const ContTopology& topology,
                    const ContForceFunction& force_function,
                    const GroupVecType& id_verletlist_group,
                    const GroupNumType& num_verletlist,
                    const CoordOffsetType& offset_verletlist_group,
                    vtkm::cont::ArrayHandle<Vec3f>& corr_ljforce,
                    vtkm::cont::ArrayHandle<Vec6f>& corr_ljvirial);

    void EAMfp(const Real& rc,
               const Vec3f& box,
               const Real& rs,
               const Id& rs_num,
               const Id& pice_num,
               const vtkm::cont::ArrayHandle<vtkm::Id>& atoms_id,
               const vtkm::cont::ArrayHandle<Vec7f>& rhor_spline,
               const vtkm::cont::ArrayHandle<Vec7f>& frho_spline,
               const ContPointLocator& locator,
               const ContForceFunction& force_function,
               const GroupVecType& id_verletlist_group,
               const GroupNumType& num_verletlist_group,
               const CoordOffsetType& offset_verletlist_group,
               vtkm::cont::ArrayHandle<Real>& EAM_fp);

    void EAMRBLForce(const Real& rc,
                     const Vec3f& box,
                     const Real& rs,
                     const Id& rs_num,
                     const Id& pice_num,
                     const vtkm::cont::ArrayHandle<vtkm::Id>& atoms_id,
                     const vtkm::cont::ArrayHandle<Vec7f>& rhor_spline,
                     const vtkm::cont::ArrayHandle<Vec7f>& z2r_spline,
                     const vtkm::cont::ArrayHandle<Real>& EAM_fp,
                     const ContPointLocator& locator,
                     const ContForceFunction& force_function,
                     const GroupVecType& id_verletlist_group,
                     const GroupNumType& num_verletlist_group,
                     const CoordOffsetType& offset_verletlist_group,
                     vtkm::cont::ArrayHandle<Vec3f>& corr_force);
    
    void SumRBLCorrForce(const vtkm::Vec3f corr_value,
                         const vtkm::cont::ArrayHandle<vtkm::Vec3f>& corr_force,
                         vtkm::cont::ArrayHandle<vtkm::Vec3f>& LJforce);

    void SumRBLCorrVirial(const Vec6f corr_value,
        const vtkm::cont::ArrayHandle<Vec6f>& corr_virial,
        vtkm::cont::ArrayHandle<Vec6f>& LJvirial);

    void LJForceWithPeriodicBC(const Real& cut_off,
                               const vtkm::cont::ArrayHandle<vtkm::Id>& atoms_id,
                               const ContPointLocator& locator,
                               const ContTopology& topology,
                               const ContForceFunction& force_function,
                               vtkm::cont::ArrayHandle<Vec3f>& LJforce);
    void SpecicalBondsLJForce(const Real& cut_off,
                              const vtkm::cont::ArrayHandle<vtkm::Id>& atoms_id,
                              const ContPointLocator& locator,
                              const ContTopology& topology,
                              const ContForceFunction& force_function,
                              const GroupIdIdType& group_ids,
                              const GroupRealIdType& group_weights,
                              vtkm::cont::ArrayHandle<Vec3f>& LJforce);
    void ComputeNearElectrostatics(const vtkm::cont::ArrayHandle<vtkm::Id>& atoms_id,
                                   const ContPointLocator& locator,
                                   const ContTopology& topology,
                                   const ContForceFunction& force_function,
                                   vtkm::cont::ArrayHandle<vtkm::Vec3f>& eleNearforce);

    void ComputeNearElectrostaticsERF(const vtkm::cont::ArrayHandle<vtkm::Id>& atoms_id,
                                      const ContStaticTable& static_table,
                                      const ContPointLocator& locator,
                                      const ContTopology& topology,
                                      const ContForceFunction& force_function,
                                      vtkm::cont::ArrayHandle<vtkm::Vec3f>& eleNearforce);

    void ComputeSpecialCoul(const Vec3f& box,
                            const vtkm::cont::ArrayHandle<vtkm::Id>& atoms_id,
                            const GroupVecType& group_vec,
                            const ContForceFunction& force_function,
                            const ContTopology& topology,
                            const ContPointLocator& locator,
                            vtkm::cont::ArrayHandle<vtkm::Vec3f>& SpecCoulforce,
                            vtkm::cont::ArrayHandle<Vec6f>& SpecCoulVirial);
    void ComputeSpecialCoulGeneral(const Vec3f& box,
                                   const vtkm::cont::ArrayHandle<vtkm::Id>& atoms_id,
                                   const GroupVecType& group_vec,
                                   const ContForceFunction& force_function,
                                   const ContTopology& topology,
                                   const ContPointLocator& locator,
                                   const GroupIdIdType& group_ids,
                                   const GroupRealIdType& group_weights,
                                   vtkm::cont::ArrayHandle<vtkm::Vec3f>& SpecCoulforce,
                                     vtkm::cont::ArrayHandle<Vec6f>& SpecCoulVirial);

    void ComputeNewRBEForce(const IdComponent& rbeP, 
                            const Vec3f& box,
                                    const vtkm::cont::ArrayHandle<vtkm::Id>& atoms_id,
                                    const vtkm::cont::ArrayHandle<vtkm::Vec3f>& psample,
                                    const vtkm::cont::ArrayHandle<vtkm::Vec2f>& whole_rhokRBE,
                                    const ContForceFunction& force_function,
                                    const ContTopology& topology,
                                    const ContPointLocator& locator,
                                    vtkm::cont::ArrayHandle<vtkm::Vec3f>& eleRBENewforce,
                                    vtkm::cont::ArrayHandle<Vec6f>& KspaceVirial);

    void InitPosition(const vtkm::cont::UnknownCellSet& cellset,
                      const vtkm::cont::CoordinateSystem& coord,
                      vtkm::cont::ArrayHandle<vtkm::Vec3f>& position);

    void InitCondtion(const vtkm::cont::ArrayHandle<vtkm::Vec3f>& position,
                      const vtkm::cont::ArrayHandle<vtkm::Vec3f>& v_array,
                      vtkm::cont::ArrayHandle<Real>& mass,
                      vtkm::cont::ArrayHandle<vtkm::Vec3f>& velocity);

    void SumFarNearForce(const vtkm::cont::ArrayHandle<vtkm::Vec3f>& eleFarNewforce,
                         const vtkm::cont::ArrayHandle<vtkm::Vec3f>& nearforce,
                         vtkm::cont::ArrayHandle<vtkm::Vec3f>& Allforce);

    void SumFarNearLJForce(const vtkm::cont::ArrayHandle<vtkm::Vec3f>& eleFarNewforce,
                           const vtkm::cont::ArrayHandle<vtkm::Vec3f>& eleNearforce,
                           const vtkm::cont::ArrayHandle<vtkm::Vec3f>& LJforce,
                           vtkm::cont::ArrayHandle<vtkm::Vec3f>& Allforce);

    void UnderdampedLangevin(const Vec3f& gaussian,
                             const Real& kBT,
                             const Real& gamma,
                             const Real& dt,
                             const vtkm::cont::ArrayHandle<Real>& mass,
                             const vtkm::cont::ArrayHandle<vtkm::Vec3f>& velocity,
                             vtkm::cont::ArrayHandle<vtkm::Vec3f>& Allforce);
    void UpdateVelocity(const Real& dt,
                        const Real& unit_factor,
                        const vtkm::cont::ArrayHandle<vtkm::Vec3f>& Allforce,
                        const vtkm::cont::ArrayHandle<Real>& mass,
                        vtkm::cont::ArrayHandle<vtkm::Vec3f>& velocity);

    void ComputerKineticEnergy(const vtkm::cont::ArrayHandle<vtkm::Vec3f>& velocity,
                               const vtkm::cont::ArrayHandle<Real>& mass,
                               vtkm::cont::ArrayHandle<Real>& sq_velocity);

    void UpdateVelocityRescale(const Real& coeff_Berendsen,
                               vtkm::cont::ArrayHandle<vtkm::Vec3f>& velocity);

    void UpdateVelocityRescale_press(const Vec3f& factor,
        vtkm::cont::ArrayHandle<vtkm::Vec3f>& velocity);

    void UpdateVelocityNoseHoover(const Real& dt,
                                  const Real& unit_factor,
                                  const Real& nosehooverxi,
                                  const vtkm::cont::ArrayHandle<vtkm::Vec3f>& Allforce,
                                  const vtkm::cont::ArrayHandle<Real>& mass,
                                  vtkm::cont::ArrayHandle<vtkm::Vec3f>& velocity);
    void ComputeChargeStructureFactorComponent(const vtkm::Vec3f& kl,
                                               const vtkm::cont::ArrayHandle<vtkm::Vec3f>& position,
                                               const vtkm::cont::ArrayHandle<Real>& charge,
                                               vtkm::cont::ArrayHandle<Real>& Density_Real,
                                               vtkm::cont::ArrayHandle<Real>& Density_Image);
    void ComputePnumberChargeStructureFactor(const Vec3f& _box,
                                            const Id& pnumber,
                                            const Id& pos_number,
                                             const vtkm::cont::ArrayHandle<vtkm::Id>& atoms_id,
                                             const vtkm::cont::ArrayHandle<vtkm::Vec3f>& position,
                                             const vtkm::cont::ArrayHandle<Real>& charge,
                                             const vtkm::cont::ArrayHandle<vtkm::Vec3f>& psample,
                                            // vtkm::cont::ArrayHandle<vtkm::Id>& psamplekey_group,
                                             vtkm::cont::ArrayHandle<Real>& rhok_Re_group,
                                             vtkm::cont::ArrayHandle<Real>& rhok_Im_group
        /*                                     GroupVecType& psamplekey_group,
                                             GroupRealType& rhok_Re_group,
                                             GroupRealType& rhok_Im_group*/);
    void ChangePnumberChargeStructureFactor(const vtkm::cont::ArrayHandle<Real>& rhok_Re,
                                             const vtkm::cont::ArrayHandle<Real>& rhok_Im,
                                             vtkm::cont::ArrayHandle<vtkm::Vec2f>& new_rhok);
    void UpdatePosition(const Real& dt,  
                        const vtkm::cont::ArrayHandle<vtkm::Vec3f>& velocity,
                        const ContPointLocator& locator,
                        vtkm::cont::ArrayHandle<vtkm::Vec3f>& position);

    void UpdatePositionFlag(const Real& dt,
                            const vtkm::cont::ArrayHandle<vtkm::Vec3f>& velocity,
                            const ContPointLocator& locator,
                            vtkm::cont::ArrayHandle<vtkm::Vec3f>& position,
                            vtkm::cont::ArrayHandle<vtkm::Id3>& position_flag);

    void ComputeNewFarElectrostatics(const Id3& Kmax,
                                     const Vec3f& box,
                                      const vtkm::cont::ArrayHandle<vtkm::Id>& atoms_id,
                                      const vtkm::cont::ArrayHandle<vtkm::Vec2f>& whole_rhok,
                                      const ContForceFunction& force_function,
                                      const ContTopology& topology,
                                      const ContPointLocator& locator,
                                      vtkm::cont::ArrayHandle<vtkm::Vec3f>& eleFarNewforce,
                                      vtkm::cont::ArrayHandle<Vec6f>& KspaceVirial);

    void ComputeSqCharge(const vtkm::cont::ArrayHandle<Real>& charge,
        vtkm::cont::ArrayHandle<Real>& SelfEnergy);

    void fix_press_berendsen(const Real& scale_factor,
        vtkm::cont::ArrayHandle<vtkm::Vec3f>& position,
        const ContPointLocator& locator);

    void ApplyPbc(const Vec3f& box,
        const Vec<Vec2f, 3>& range,
        vtkm::cont::ArrayHandle<vtkm::Vec3f>& position,
        const ContPointLocator& locator);

    void ApplyPbcFlag(const Vec3f& box,
        const Vec<Vec2f, 3>& range,
        vtkm::cont::ArrayHandle<vtkm::Vec3f>& position,
        const ContPointLocator& locator,
        vtkm::cont::ArrayHandle<vtkm::Id3>& position_flag);
    void X2Lamda(const Vec6f& h_inv,
        const vtkm::Vec<vtkm::Range, 3>& range,
        vtkm::cont::ArrayHandle<vtkm::Vec3f>& position);

    void Lamda2X(const Vec6f& h,
        const vtkm::Vec<vtkm::Range, 3>& range,
        vtkm::cont::ArrayHandle<vtkm::Vec3f>& position);

    void UnWarpPostion(const Vec3f& _box,
        const vtkm::cont::ArrayHandle<vtkm::Id3>& position_flag,
        vtkm::cont::ArrayHandle<vtkm::Vec3f>& unwarp_position);

    void ComputeCOM(const vtkm::cont::ArrayHandle<vtkm::Vec3f>& unwarp_position,
        const vtkm::cont::ArrayHandle<Real>& mass,
        vtkm::cont::ArrayHandle<vtkm::Vec3f>& com);

    void ComputeVCOM(const vtkm::cont::ArrayHandle<vtkm::Vec3f>& velocity,
        const vtkm::cont::ArrayHandle<Real>& mass,
        vtkm::cont::ArrayHandle<vtkm::Vec3f>& vcom);

    void ComputeOmega(const vtkm::Vec3f& com,
        const vtkm::cont::ArrayHandle<vtkm::Vec3f>& unwarp_position,
        const vtkm::cont::ArrayHandle<vtkm::Vec3f>& velocity,
        const vtkm::cont::ArrayHandle<Real>& mass,
        const ContForceFunction& force_function,
        vtkm::cont::ArrayHandle<vtkm::Vec3f>& omega);

    void FixMomentum(const vtkm::Vec3f& pcom,
        const vtkm::Vec3f& vcom,
        const vtkm::Vec3f& omega,
        const vtkm::cont::ArrayHandle<vtkm::Vec3f>& unwarp_position,
        const vtkm::cont::ArrayHandle<Real>& mass,
        const ContForceFunction& force_function,
        vtkm::cont::ArrayHandle<vtkm::Vec3f>& velocity);

}
