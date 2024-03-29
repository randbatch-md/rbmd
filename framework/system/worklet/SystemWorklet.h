﻿#pragma once
#include "locator/ContPointLocator.h"
#include "forceFunction/ContForceFunction.h"
#include "topology/ContTopology.h"
#include <vtkm/cont/ArrayHandleGroupVecVariable.h>
#include "staticTable/ContStaticTable.h"

namespace SystemWorklet 
{
    void ComputeNeighbours(const Real& cut_off,
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
                           const vtkm::cont::ArrayHandle<vtkm::Id>& atoms_id,
                           const ContPointLocator& locator,
                           const ContTopology& topology,
                           const ContForceFunction& force_function,
                           const GroupVecType& Group_j,
                           const vtkm::cont::ArrayHandle<vtkm::Id>& num_j,
                           const CoordOffsetType& coord_offset_j,
                           vtkm::cont::ArrayHandle<Vec3f>& nearforce);

   void LJForceVerlet(const Real& cut_off,
                        const vtkm::cont::ArrayHandle<vtkm::Id>& atoms_id,
                        const ContPointLocator& locator,
                        const ContTopology& topology,
                        const ContForceFunction& force_function,
                        const GroupVecType& Group_j,
                        const vtkm::cont::ArrayHandle<vtkm::Id>& num_j,
                        const CoordOffsetType& coord_offset_j,
                        vtkm::cont::ArrayHandle<Vec3f>& LJforce);
    
    void NearForceRBLERF(const Id& rs_num,
                         const Id& pice_num,
                         const Real& qqr2e,
                         const vtkm::cont::ArrayHandle<vtkm::Id>& atoms_id,
                         const ContPointLocator& locator,
                         const ContTopology& topology,
                         const ContForceFunction& force_function,
                         const ContStaticTable& static_table,
                         const GroupVecType& id_verletlist_group,
                         const GroupNumType& num_verletlist,
                         const CoordOffsetType& offset_verletlist_group,
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
                    const vtkm::cont::ArrayHandle<vtkm::Id>& atoms_id,
                    const ContPointLocator& locator,
                    const ContTopology& topology,
                    const ContForceFunction& force_function,
                    const GroupVecType& id_verletlist_group,
                    const GroupNumType& num_verletlist,
                    const CoordOffsetType& offset_verletlist_group,
                    vtkm::cont::ArrayHandle<Vec3f>& corr_ljforce);
    
    void SumRBLCorrForce(const vtkm::Vec3f corr_value,
                         const vtkm::cont::ArrayHandle<vtkm::Vec3f>& corr_force,
                         vtkm::cont::ArrayHandle<vtkm::Vec3f>& LJforce);

    void LJForceWithPeriodicBC(const Real& cut_off,
                               const vtkm::cont::ArrayHandle<vtkm::Id>& atoms_id,
                               const ContPointLocator& locator,
                               const ContTopology& topology,
                               const ContForceFunction& force_function,
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

    void ComputeSpecialCoul(const Real& Vlength,
                            const vtkm::cont::ArrayHandle<vtkm::Id>& atoms_id,
                            const GroupVecType& group_vec,
                            const ContForceFunction& force_function,
                            const ContTopology& topology,
                            const ContPointLocator& locator,
                            vtkm::cont::ArrayHandle<vtkm::Vec3f>& SpecCoulforce);
    void ComputeNewRBEForce(const IdComponent& rbeP, 
                                    const vtkm::cont::ArrayHandle<vtkm::Id>& atoms_id,
                                    const vtkm::cont::ArrayHandle<vtkm::Vec3f>& psample,
                                    const vtkm::cont::ArrayHandle<vtkm::Vec2f>& whole_rhokRBE,
                                    const ContForceFunction& force_function,
                                    const ContTopology& topology,
                                    const ContPointLocator& locator,
                                    vtkm::cont::ArrayHandle<vtkm::Vec3f>& eleRBENewforce);

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
    void UpdatePosition(const Real& dt,  
                        const vtkm::cont::ArrayHandle<vtkm::Vec3f>& velocity,
                        const ContPointLocator& locator,
                        vtkm::cont::ArrayHandle<vtkm::Vec3f>& position);

    void UpdatePositionFlag(const Real& dt,
                            const vtkm::cont::ArrayHandle<vtkm::Vec3f>& velocity,
                            const ContPointLocator& locator,
                            vtkm::cont::ArrayHandle<vtkm::Vec3f>& position,
                            vtkm::cont::ArrayHandle<vtkm::Id3>& position_flag);

    void ComputeNewFarElectrostatics(const IdComponent& Kmax,
                                      const vtkm::cont::ArrayHandle<vtkm::Id>& atoms_id,
                                      const vtkm::cont::ArrayHandle<vtkm::Vec2f>& whole_rhok,
                                      const ContForceFunction& force_function,
                                      const ContTopology& topology,
                                      const ContPointLocator& locator,
                                      vtkm::cont::ArrayHandle<vtkm::Vec3f>& eleFarNewforce);
}