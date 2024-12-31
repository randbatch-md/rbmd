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
#include "locator/ContPointLocator.h"
#include "forceFunction/ContForceFunction.h"
#include "topology/ContTopology.h"
#include <vtkm/cont/ArrayHandleGroupVecVariable.h>

using GroupIdIdType = vtkm::cont::ArrayHandleGroupVecVariable<vtkm::cont::ArrayHandle<vtkm::Id>,
                                                              vtkm::cont::ArrayHandle<vtkm::Id>>;

using GroupRealIdType = vtkm::cont::ArrayHandleGroupVecVariable<vtkm::cont::ArrayHandle<Real>,
                                                                vtkm::cont::ArrayHandle<vtkm::Id>>;

namespace OutPut
{
    void ComputeNeighbours(const Real& cut_off,
        const Vec3f& box,
        const vtkm::cont::ArrayHandle<vtkm::Id>& atoms_id,
        const ContPointLocator& locator,
        GroupVecType& id_verletlist_group,
        vtkm::cont::ArrayHandle<vtkm::Id>& num_verletlist,
        CoordOffsetType& offset_verletlist_group);

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

    void LJCoulVerlet(const Real& cut_off,
        const Real& alpha,
        const Vec3f& box,
        const vtkm::cont::ArrayHandle<vtkm::Id>& atoms_id,
        const ContPointLocator& locator,
        const ContTopology& topology,
        const ContForceFunction& force_function,
        const GroupVecType& Group_j,
        const vtkm::cont::ArrayHandle<vtkm::Id>& num_j,
        const CoordOffsetType& coord_offset_j,
        vtkm::cont::ArrayHandle<Real>& LJCoul);

void EAM_rho(const Real& eam_cut_off,
             const Vec3f& box,
             const vtkm::cont::ArrayHandle<vtkm::Id>& atoms_id,
             const vtkm::cont::ArrayHandle<Vec7f>& rhor_spline,
             const ContPointLocator& locator,
             const ContTopology& topology,
             const ContForceFunction& force_function,
             vtkm::cont::ArrayHandle<Real>& EAM_rho);

void EAM_EmbeddingEnergy(const vtkm::cont::ArrayHandle<vtkm::Id>& atoms_id,
                         const vtkm::cont::ArrayHandle<Real>& EAM_rho,
                         const vtkm::cont::ArrayHandle<Vec7f>& frho_spline,
                         const ContPointLocator& locator,
                         const ContTopology& topology,
                         const ContForceFunction& force_function,
                         vtkm::cont::ArrayHandle<Real>& embedding_energy);

void EAM_PairEnergy(const Real& eam_cut_off,
                    const Vec3f& box,
                    const vtkm::cont::ArrayHandle<vtkm::Id>& atoms_id,
                    const vtkm::cont::ArrayHandle<Vec7f>& z2r_spline,
                    const ContPointLocator& locator,
                    const ContTopology& topology,
                    const ContForceFunction& force_function,
                    vtkm::cont::ArrayHandle<Real>& pair_energy);

void ComputePotentialEnergy(const Real& cutoff,
                            const Vec3f& box,
                            const vtkm::cont::ArrayHandle<vtkm::Id>& atoms_id,
                            const ContPointLocator& locator,
                            const ContTopology& topology,
                            const ContForceFunction& force_function,
                            vtkm::cont::ArrayHandle<Real>& potential_energy);

void ComputeSpecialBondsLJPotential(const Real& cutoff,
                                    const Vec3f& box,
                                    const vtkm::cont::ArrayHandle<vtkm::Id>& atoms_id,
                                    const ContPointLocator& locator,
                                    const ContTopology& topology,
                                    const ContForceFunction& force_function,
                                    const GroupIdIdType& group_ids,
                                    const GroupRealIdType& group_weights,
                                    vtkm::cont::ArrayHandle<Real>& potential_energy);

void ComputeSpecialBondsLJPotential_fix(const Real& cutoff,
    const Real nearalpha,
    const Vec3f& box,
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
    vtkm::cont::ArrayHandle<Real>& energy_coul);

void ComputeNearElePotential(const Real& cutoff,
                             const Real& alpha,
                             const Vec3f& box,
                             const vtkm::cont::ArrayHandle<vtkm::Id>& atoms_id,
                             const ContPointLocator& locator,
                             const ContTopology& topology,
                             const ContForceFunction& force_function,
                             vtkm::cont::ArrayHandle<Real>& near_potential_energy);

void ComputeDensity(const vtkm::Vec3f& K,
                    const vtkm::cont::ArrayHandle<vtkm::Vec3f>& position,
                    const vtkm::cont::ArrayHandle<Real>& charge,
                    vtkm::cont::ArrayHandle<Real>& DensityReal,
                    vtkm::cont::ArrayHandle<Real>& DensityImage);
    
void ComputeSqCharge(const vtkm::cont::ArrayHandle<Real>& charge,
                     vtkm::cont::ArrayHandle<Real>& SelfEnergy);

using GroupVecType = vtkm::cont::ArrayHandleGroupVecVariable<vtkm::cont::ArrayHandle<vtkm::Id>,
                                                             vtkm::cont::ArrayHandle<vtkm::Id>>;
void ComputeSpecialFarCoul(const Vec3f& box,
                           const vtkm::cont::ArrayHandle<Id>& atoms_id,
                           const GroupVecType& group_vec,
                           const ContPointLocator& locator,
                           const ContTopology& topology,
                           const ContForceFunction& force_function,
                           vtkm::cont::ArrayHandle<Real>& SpecFarEnergy);
void ComputeSpecialBondsCoul(const Vec3f& box,
                             const vtkm::cont::ArrayHandle<Id>& atoms_id,
                             const GroupVecType& group_vec,
                             const ContPointLocator& locator,
                             const ContTopology& topology,
                             const ContForceFunction& force_function,
                             const GroupIdIdType& group_ids,
                             const GroupRealIdType& group_weights,
                             vtkm::cont::ArrayHandle<Real>& SpecFarEnergy);

 void ComputeRDF(const Id& num_center_pos,
                const Real& _rho,
                const vtkm::cont::ArrayHandle<Real>& radius,
                const vtkm::cont::ArrayHandle<vtkm::Vec3f>& center_position,
                const vtkm::cont::ArrayHandle<vtkm::Vec3f>& target_position,
                const vtkm::cont::ArrayHandle<vtkm::Id>& center_ids,
                const vtkm::cont::ArrayHandle<vtkm::Id>& molecule_id,
                const ContPointLocator& locator,
                vtkm::cont::ArrayHandle<Real>& rdf);

void ComputerKineticEnergy(const vtkm::cont::ArrayHandle<vtkm::Vec3f>& velocity,
                           const vtkm::cont::ArrayHandle<Real>& mass,
                           vtkm::cont::ArrayHandle<Real>& sq_velocity);

void ComputeMSD(const Vec3f& _box,
                const vtkm::cont::ArrayHandle<vtkm::Vec3f>& original_position,
                const vtkm::cont::ArrayHandle<vtkm::Vec3f>& current_pts_position,
                const vtkm::cont::ArrayHandle<vtkm::Id3>& position_flag,
                vtkm::cont::ArrayHandle<vtkm::Vec4f>& MSDoutput);

void ComputePrepareMSD(const Vec3f& _box,
                       const vtkm::cont::ArrayHandle<vtkm::Vec3f>& temp_MSD_position,
                       const vtkm::cont::ArrayHandle<vtkm::Vec3f>& current_pts_position,
                       const vtkm::cont::ArrayHandle<vtkm::Id3>& position_flag,
                       vtkm::cont::ArrayHandle<vtkm::Vec3f>& diff_MSD_position);

void ComputeVACF(const vtkm::cont::ArrayHandle<vtkm::Vec3f>& original_velocity,
                 const vtkm::cont::ArrayHandle<vtkm::Vec3f>& current_velocity,
                 vtkm::cont::ArrayHandle<vtkm::Vec4f>& VACFoutput);

namespace atoms
{
void ComputeRDF(const Id& num_center_pos,
                const Real& _rho,
                const vtkm::cont::ArrayHandle<Real>& radius,
                const vtkm::cont::ArrayHandle<vtkm::Vec3f>& center_position,
                const vtkm::cont::ArrayHandle<vtkm::Vec3f>& target_position,
                const ContPointLocator& locator,
                vtkm::cont::ArrayHandle<Real>& rdf);
}
}