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
void ComputePotentialEnergy(const Real& cutoff,
                             const vtkm::cont::ArrayHandle<vtkm::Id>& atoms_id,
                             const ContPointLocator& locator,
                             const ContTopology& topology,
                             const ContForceFunction& force_function,
                             vtkm::cont::ArrayHandle<Real>& potential_energy);

void ComputeSpecialBondsLJPotential(const Real& cutoff,
                                    const vtkm::cont::ArrayHandle<vtkm::Id>& atoms_id,
                                    const ContPointLocator& locator,
                                    const ContTopology& topology,
                                    const ContForceFunction& force_function,
                                    const GroupIdIdType& group_ids,
                                    const GroupRealIdType& group_weights,
                                    vtkm::cont::ArrayHandle<Real>& potential_energy);

void ComputeNearElePotential(const Real& cutoff,
                             const Real& alpha,
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
void ComputeSpecialFarCoul(const Real& Vlength,
                           const vtkm::cont::ArrayHandle<Id>& atoms_id,
                           const GroupVecType& group_vec,
                           const ContPointLocator& locator,
                           const ContTopology& topology,
                           const ContForceFunction& force_function,
                           vtkm::cont::ArrayHandle<Real>& SpecFarEnergy); 
void ComputeSpecialBondsCoul(const Real& Vlength,
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
                const vtkm::cont::ArrayHandle<vtkm::Id>& center_ids,
                const vtkm::cont::ArrayHandle<vtkm::Id>& molecule_id,
                const ContPointLocator& locator,
                vtkm::cont::ArrayHandle<Real>& rdf);

void ComputerKineticEnergy(const vtkm::cont::ArrayHandle<vtkm::Vec3f>& velocity,
                           const vtkm::cont::ArrayHandle<Real>& mass,
                           vtkm::cont::ArrayHandle<Real>& sq_velocity);

void ComputeMSD(const Real& _Vlength,
                const vtkm::cont::ArrayHandle<vtkm::Vec3f>& original_position,
                const vtkm::cont::ArrayHandle<vtkm::Vec3f>& current_pts_position,
                const vtkm::cont::ArrayHandle<vtkm::Id3>& position_flag,
                vtkm::cont::ArrayHandle<vtkm::Vec4f>& MSDoutput);

void ComputePrepareMSD(const Real& _Vlength,
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