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

#include "output/worklet/OutPutWorklet.h"
#include "Executioner.h"
#include <vtkm/worklet/WorkletMapField.h>
#include <vtkm/cont/Algorithm.h>
#include "Application.h"
#include <fmt/format.h>
#include "math/Math.h"

namespace OutPut
{
    struct ComputeNeighboursWorklet : vtkm::worklet::WorkletMapField
    {
        ComputeNeighboursWorklet(const Real& cut_off, const Vec3f& box)
            : _cut_off(cut_off)
            , _box(box)
        {
        }

        using ControlSignature = void(FieldIn atoms_id,
            ExecObject locator,
            FieldInOut id_verletlist_group,
            FieldOut num_verletlist,
            FieldInOut offset_verletlist_group);
        using ExecutionSignature = void(_1, _2, _3, _4, _5);

        template<typename NeighbourGroupVecType, typename CoordOffsetj>
        VTKM_EXEC void operator()(const Id atoms_id,
            const ExecPointLocator& locator,
            NeighbourGroupVecType& id_verletlist,
            Id& num_verletlist,
            CoordOffsetj& offset_verletlist) const
        {
            Id index = 0;
            auto p_i = locator.GetPtsPosition(atoms_id);
            vtkm::Id3 p_i_cell = locator.InWhichCell(p_i);
            const auto num_cycles = locator.GetNumCycles();
            for (Id i = -num_cycles; i <= num_cycles; i++)
            {
                for (Id j = -num_cycles; j <= num_cycles; j++)
                {
                    for (Id k = -num_cycles; k <= num_cycles; k++)
                    {
                        auto neighbor_ijk = p_i_cell + Id3{ i, j, k };
                        auto ijk = locator.PeriodicIndexOffset(neighbor_ijk);
                        auto coord_offset = locator.PeriodicCoordOffset(ijk - neighbor_ijk);
                        auto num_pts = locator.NumberPtsInCell(ijk);

                        for (Id p = 0; p < num_pts; p++)
                        {
                            auto pts_id_j = locator.PtsInCell(ijk, p);
                            auto p_j = locator.GetPtsPosition(pts_id_j) - coord_offset;
                            auto r_ij = p_j - p_i;
                            //auto r_ij = locator.MinDistanceVec(p_i, p_j, _box);
                            const Real dis_2 = r_ij[0] * r_ij[0] + r_ij[1] * r_ij[1] + r_ij[2] * r_ij[2];
                            const Real _cut_off_2 = _cut_off * _cut_off;
                            if (_cut_off_2 - dis_2 > 0.0001 && dis_2 > 0.0001)
                            {
                                id_verletlist[index] = pts_id_j;
                                offset_verletlist[index] = coord_offset;
                                index++;
                            }
                        }
                    }
                }
            }
            num_verletlist = index;
        }
        Real _cut_off;
        Vec3f _box;
    };

    struct ComputeEAMrhoWorklet : vtkm::worklet::WorkletMapField
    {
      ComputeEAMrhoWorklet(const Real& eam_cut_off, const Vec3f& box)
        : _eam_cut_off(eam_cut_off)
        , _box(box)
      {
      }
    
      using ControlSignature = void(FieldIn atoms_id,
                                    WholeArrayIn rhor_spline,
                                    ExecObject locator,
                                    ExecObject topology,
                                    ExecObject force_function,
                                    FieldOut EAM_rho);
      using ExecutionSignature = void(_1, _2, _3, _4, _5, _6);
    
      template<typename SpileType>
      VTKM_EXEC void operator()(const Id atoms_id,
                                const SpileType& rhor_spline,
                                const ExecPointLocator& locator,
                                const ExecTopology& topology,
                                const ExecForceFunction& force_function,
                                Real& EAM_rho) const
      {
        Real rho = 0;
    
        auto function = [&](const Vec3f& p_i, const Vec3f& p_j, const Id& pts_id_j)
        {
          //auto r_ij = p_j - p_i;
          auto r_ij = locator.MinDistanceVec(p_i, p_j, _box);
          rho += force_function.ComputeEAMrhoOUT(_eam_cut_off, r_ij, rhor_spline);
        };
        locator.ExecuteOnNeighbor(atoms_id, function);
        EAM_rho = rho;
      }
      Real _eam_cut_off;
      Vec3f _box;
    };

    struct ComputeEmbeddingEnergyWorklet : vtkm::worklet::WorkletMapField
    {
      ComputeEmbeddingEnergyWorklet() {}
    
      using ControlSignature = void(FieldIn atoms_id,
                                    WholeArrayIn EAM_rho,
                                    WholeArrayIn frho_spline,
                                    ExecObject locator,
                                    ExecObject topology,
                                    ExecObject force_function,
                                    FieldOut embedding_energy);
      using ExecutionSignature = void(_1, _2, _3, _4, _5, _6, _7);
    
      template<typename EAM_rhoype, typename SpileType>
      VTKM_EXEC void operator()(const Id atoms_id,
                                const EAM_rhoype& EAM_rho,
                                const SpileType& frho_spline,
                                const ExecPointLocator& locator,
                                const ExecTopology& topology,
                                const ExecForceFunction& force_function,
                                Real& embedding_energy) const
      {
        embedding_energy = force_function.ComputeEmbeddingEnergy(atoms_id, EAM_rho, frho_spline);
      }
    };

    struct ComputePairEnergyWorklet : vtkm::worklet::WorkletMapField
    {
      ComputePairEnergyWorklet(const Real& eam_cut_off, const Vec3f& box)
        : _eam_cut_off(eam_cut_off)
        , _box(box)
      {
      }
    
      using ControlSignature = void(FieldIn atoms_id,
                                    WholeArrayIn z2r_spline,
                                    ExecObject locator,
                                    ExecObject topology,
                                    ExecObject force_function,
                                    FieldOut pair_energy);
      using ExecutionSignature = void(_1, _2, _3, _4, _5, _6);
    
      template<typename z2rSpileType>
      VTKM_EXEC void operator()(const Id atoms_id,
                                const z2rSpileType& z2r_spline,
                                const ExecPointLocator& locator,
                                const ExecTopology& topology,
                                const ExecForceFunction& force_function,
                                Real& pair_energy) const
      {
        pair_energy = 0;
        auto function = [&](const Vec3f& p_i, const Vec3f& p_j, const Id& pts_id_j)
        {
          //auto r_ij = p_j - p_i;
          auto r_ij = locator.MinDistanceVec(p_i, p_j, _box);
          pair_energy += force_function.ComputePairEnergy(_eam_cut_off, r_ij, z2r_spline);
        };
        locator.ExecuteOnNeighbor(atoms_id, function);
      }
      Real _eam_cut_off;
      Vec3f _box;
    };

    struct ComputeDensityWorklet : vtkm ::worklet::WorkletMapField
    {
      ComputeDensityWorklet(const Vec3f& k)
        : _K(k)
      {
      }
    
      using ControlSignature = void(FieldIn current_pts,
                                    FieldIn current_charge,
                                    FieldOut density_Re,
                                    FieldOut density_Im);
    
      using ExecutionSignature = void(_1, _2, _3, _4);
    
      template<typename CoordType, typename Densitytype>
      VTKM_EXEC void operator()(const CoordType& current_pts,
                                const Real& current_charge,
                                Densitytype& density_Re,
                                Densitytype& density_Im) const
      {
        Vec3f r_i = current_pts;
        density_Re = current_charge * vtkm::Cos(vtkm::Dot(_K, r_i));
        density_Im = current_charge * vtkm::Sin(vtkm::Dot(_K, r_i));
      }
    
    
      Vec3f _K;
    };

    struct ComputeRDFWorklet : vtkm::worklet::WorkletMapField
    {
      ComputeRDFWorklet(const Id& num_center_pos, const Real& rho)
        : _num_center_pos(num_center_pos)
        , _rho(rho)
      {
      }
      using ControlSignature = void(FieldIn current_radius,
                                    WholeArrayIn center_pts,
                                    WholeArrayIn target_pts,
                                    WholeArrayIn center_ids,
                                    WholeArrayIn molecule_id,
                                    ExecObject locator,
                                    FieldInOut rdf);
      using ExecutionSignature = void(_1, _2, _3, _4, _5, _6, _7);

      template<typename RadiusType, typename WholePtsType, typename RdfType, typename IdsType>
      VTKM_EXEC void operator()(const RadiusType& current_radius,
                                const WholePtsType& center_pts,
                                const WholePtsType& target_pts,
                                const IdsType& center_ids,
                                const IdsType& molecule_id,
                                const ExecPointLocator& locator,
                                RdfType& rdf) const
      {
        auto GetNumber = [=](const Vec3f& p_i, const Id& center_molecule_id) -> Id
        {
          Id number = 0;

          auto current_ijk = locator.InWhichCell(p_i);
          vtkm::Id3 dx = vtkm::Vec3f{ current_radius } / locator.GetDxdydz() + 1;
          auto ijk_upper = current_ijk + dx;
          auto ijk_lower = current_ijk - dx;

          for (Id i = ijk_lower[0]; i <= ijk_upper[0]; i++)
          {
            for (Id j = ijk_lower[1]; j <= ijk_upper[1]; j++)
            {
              for (Id k = ijk_lower[2]; k <= ijk_upper[2]; k++)
              {
                auto neighbor_ijk = Id3{ i, j, k };
                auto ijk = locator.PeriodicIndexOffset(neighbor_ijk);
                auto coord_offset = locator.PeriodicCoordOffset(ijk - neighbor_ijk);
                auto num_pts = locator.NumberPtsInCell(ijk);


                for (Id p = 0; p < num_pts; p++)
                {
                  auto pts_id = locator.PtsInCell(ijk, p);
                  auto p_j = target_pts.Get(pts_id) - coord_offset;
                  Real dis = vtkm::Magnitude(p_j - p_i);
                  const auto& target_molecule_id = molecule_id.Get(pts_id);

                  if (center_molecule_id != target_molecule_id && dis > current_radius &&
                      dis < current_radius + current_radius * 0.01)
                  {
                    number++;
                  }
                }
              }
            }
          }

          return number;
        };

        auto total_num = 0;
        for (auto i = 0; i < _num_center_pos; ++i)
        {
              const auto& p_i = center_pts.Get(i);
              const auto& center_molecule_id = molecule_id.Get(center_ids.Get(i));
              total_num += GetNumber(p_i, center_molecule_id);
        }
        auto dr3 = vtkm::Pow((current_radius + current_radius * 0.01), 3);
        auto r3 = vtkm::Pow(current_radius, 3);
        rdf += total_num /
          (_num_center_pos * _rho *
           ((4 / 3.0) * vtkm::Pif() * dr3 - ((4 / 3.0) * vtkm::Pif() * r3)));
      }

      Id _num_center_pos;
      Real _rho;
    };

    struct ComputerKineticEnergyWorklet : vtkm::worklet::WorkletMapField
    {
      using ControlSignature = void(FieldIn velocity, FieldIn mass, FieldOut sq_velocity);
      using ExecutionSignature = void(_1, _2, _3);

      VTKM_EXEC void operator()(const Vec3f& velocity, const Real& mass, Real& sq_velocity) const
      {
        sq_velocity = mass *
          (velocity[0] * velocity[0] + velocity[1] * velocity[1] + velocity[2] * velocity[2]);
      }
    };

    struct ComputeMSDWorklet : vtkm::worklet::WorkletMapField
    {
      ComputeMSDWorklet(const Vec3f& box)
        : _box(box)
      {
      }

      using ControlSignature = void(FieldIn original_pts,
                                    FieldIn current_pts,
                                    FieldIn position_flag,
                                    FieldOut MSDoutput);
      using ExecutionSignature = void(_1, _2, _3,_4);

      template<typename CoordType, typename PositionFlagType>
      VTKM_EXEC void operator()(const CoordType& original_pts,
                                const CoordType& current_pts,
                                const PositionFlagType& position_flag,
                                vtkm::Vec4f& MSDoutput) const
      {
        Vec3f temp_position,d_position;
        //temp_position = current_pts;
        for (int i = 0; i < 3; i++)
        {
              temp_position[i] = current_pts[i] + position_flag[i] * _box[i];
        }
        //temp_position = current_pts + position_flag * _Vlength;
        d_position = temp_position - original_pts;

        MSDoutput[0] += d_position[0] * d_position[0];
        MSDoutput[1] += d_position[1] * d_position[1];
        MSDoutput[2] += d_position[2] * d_position[2];
        MSDoutput[3] += d_position[0] * d_position[0] + d_position[1] * d_position[1] +
          d_position[2] * d_position[2];
      }

     Vec3f _box;
    };

    struct ComputePrepareMSDWorklet : vtkm::worklet::WorkletMapField
    {
     ComputePrepareMSDWorklet(const Vec3f& box)
       : _box(box)
      {
      }

      using ControlSignature = void(FieldIn temp_MSD_pts,
                                    FieldIn current_pts,
                                    FieldIn position_flag,
                                    FieldOut diff_MSD_pts);
      using ExecutionSignature = void(_1, _2, _3, _4);

      template<typename PositionFlagType>
      VTKM_EXEC void operator()(const vtkm::Vec3f& temp_MSD_pts,
                                const vtkm::Vec3f& current_pts,
                                const PositionFlagType& position_flag,
                                vtkm::Vec3f& diff_MSD_pts) const
      {
        Vec3f d_position;
        for (int i = 0; i < 3; i++)
        {
              d_position[i] = current_pts[i] + position_flag[i] * _box[i] - temp_MSD_pts[i];
        }
        //d_position = current_pts + position_flag * _Vlength - temp_MSD_pts;
       
        for (int i = 0; i < 3; i++)
        {
          if (vtkm::Abs(d_position[i]) < 1.0)
          {
            diff_MSD_pts[i] = d_position[i];
          }
          else
            diff_MSD_pts[i] = d_position[i] - position_flag[i] * _box[i];
        }
        
        
      }

     Vec3f _box;
    };

    struct ComputeVACFWorklet : vtkm::worklet::WorkletMapField
    {
      using ControlSignature = void(FieldIn original_velocity,
                                    FieldIn current_velocity,
                                    FieldOut VACFoutput);
      using ExecutionSignature = void(_1, _2, _3);

      VTKM_EXEC void operator()(const vtkm::Vec3f& original_velocity,
                                const vtkm::Vec3f& current_velocity,
                                vtkm::Vec4f& VACFoutput) const
      {

        vtkm::Float64 vx_sq, vy_sq, vz_sq;
        vx_sq = original_velocity[0] * current_velocity[0];
        vy_sq = original_velocity[1] * current_velocity[1];
        vz_sq = original_velocity[2] * current_velocity[2];

        VACFoutput[0] += vx_sq;
        VACFoutput[1] += vy_sq;
        VACFoutput[2] += vz_sq;
        VACFoutput[3] += vx_sq + vy_sq + vz_sq;
      }
    };

    void EAM_rho(const Real& eam_cut_off,
                 const Vec3f& box,
                 const vtkm::cont::ArrayHandle<vtkm::Id>& atoms_id,
                 const vtkm::cont::ArrayHandle<Vec7f>& rhor_spline,
                 const ContPointLocator& locator,
                 const ContTopology& topology,
                 const ContForceFunction& force_function,
                 vtkm::cont::ArrayHandle<Real>& EAM_rho)
    {
      vtkm::cont::Invoker{}(ComputeEAMrhoWorklet{ eam_cut_off, box },
                            atoms_id,
                            rhor_spline,
                            locator,
                            topology,
                            force_function,
                            EAM_rho);
    }

    void EAM_EmbeddingEnergy(const vtkm::cont::ArrayHandle<vtkm::Id>& atoms_id,
                             const vtkm::cont::ArrayHandle<Real>& EAM_rho,
                             const vtkm::cont::ArrayHandle<Vec7f>& frho_spline,
                             const ContPointLocator& locator,
                             const ContTopology& topology,
                             const ContForceFunction& force_function,
                             vtkm::cont::ArrayHandle<Real>& embedding_energy)
    {

      vtkm::cont::Invoker{}(ComputeEmbeddingEnergyWorklet{},
                            atoms_id,
                            EAM_rho,
                            frho_spline,
                            locator,
                            topology,
                            force_function,
                            embedding_energy);
    }

    void EAM_PairEnergy(const Real& eam_cut_off,
                        const Vec3f& box,
                        const vtkm::cont::ArrayHandle<vtkm::Id>& atoms_id,
                        const vtkm::cont::ArrayHandle<Vec7f>& z2r_spline,
                        const ContPointLocator& locator,
                        const ContTopology& topology,
                        const ContForceFunction& force_function,
                        vtkm::cont::ArrayHandle<Real>& pair_energy)
    {

      vtkm::cont::Invoker{}(ComputePairEnergyWorklet{ eam_cut_off, box },
                            atoms_id,
                            z2r_spline,
                            locator,
                            topology,
                            force_function,
                            pair_energy);
    }

    //Statistical PotentialEnergy
    void ComputeNeighbours(const Real& cut_off,
        const Vec3f& box,
        const vtkm::cont::ArrayHandle<vtkm::Id>& atoms_id,
        const ContPointLocator& locator,
        GroupVecType& id_verletlist_group,
        vtkm::cont::ArrayHandle<vtkm::Id>& num_verletlist,
        CoordOffsetType& offset_verletlist_group)
    {
        vtkm::cont::Invoker{}(ComputeNeighboursWorklet{ cut_off, box },
            atoms_id,
            locator,
            id_verletlist_group,
            num_verletlist,
            offset_verletlist_group);
    }

    void ComputeDensity(const vtkm::Vec3f& K,
                        const vtkm::cont::ArrayHandle<vtkm::Vec3f>& position,
                        const vtkm::cont::ArrayHandle<Real>& charge,
                        vtkm::cont::ArrayHandle<Real>& DensityReal,
                        vtkm::cont::ArrayHandle<Real>& DensityImage)
    {
      vtkm::cont::Invoker{}(ComputeDensityWorklet{ K }, position, charge, DensityReal, DensityImage);

    }

    void ComputeRDF(const Id& num_center_pos,
                    const Real& _rho,
                    const vtkm::cont::ArrayHandle<Real>& radius,
                    const vtkm::cont::ArrayHandle<vtkm::Vec3f>& center_position,
                    const vtkm::cont::ArrayHandle<vtkm::Vec3f>& target_position,
                    const vtkm::cont::ArrayHandle<vtkm::Id>& center_ids,
                    const vtkm::cont::ArrayHandle<vtkm::Id>& molecule_id,
                    const ContPointLocator& locator,
                    vtkm::cont::ArrayHandle<Real>& rdf)
    {
      vtkm::cont::Invoker{}(ComputeRDFWorklet{ num_center_pos, _rho },
                            radius,
                            center_position,
                            target_position,
                            center_ids,
                            molecule_id,
                            locator,
                            rdf);
    }

    // Statistical Kinetic energy
    void ComputerKineticEnergy(const vtkm::cont::ArrayHandle<vtkm::Vec3f>& velocity,
                               const vtkm::cont::ArrayHandle<Real>& mass,
                               vtkm::cont::ArrayHandle<Real>& sq_velocity)
    {
      vtkm::cont::Invoker{}(ComputerKineticEnergyWorklet{}, velocity, mass, sq_velocity);
    }

    // Statistical MSD
    void ComputeMSD(const Vec3f& _box,
                    const vtkm::cont::ArrayHandle<vtkm::Vec3f>& original_position,
                    const vtkm::cont::ArrayHandle<vtkm::Vec3f>& current_pts_position,
                    const vtkm::cont::ArrayHandle<vtkm::Id3>& position_flag,
                    vtkm::cont::ArrayHandle<vtkm::Vec4f>& MSDoutput)
    {
      vtkm::cont::Invoker{}(ComputeMSDWorklet{ _box },
                            original_position,
                            current_pts_position,
                            position_flag,
                            MSDoutput);
    }

    // Statistical MSD
    void ComputePrepareMSD(const Real& _Vlength,
                           const vtkm::cont::ArrayHandle<vtkm::Vec3f>& temp_MSD_position,
                           const vtkm::cont::ArrayHandle<vtkm::Vec3f>& current_pts_position,
                           const vtkm::cont::ArrayHandle<vtkm::Id3>& position_flag,
                           vtkm::cont::ArrayHandle<vtkm::Vec3f>& diff_MSD_position)
    {
      vtkm::cont::Invoker{}(ComputePrepareMSDWorklet{ _Vlength },
                            temp_MSD_position,
                            current_pts_position,
                            position_flag,
                            diff_MSD_position);
    }

    void ComputeVACF(const vtkm::cont::ArrayHandle<vtkm::Vec3f>& original_velocity,
                     const vtkm::cont::ArrayHandle<vtkm::Vec3f>& current_velocity,
                     vtkm::cont::ArrayHandle<vtkm::Vec4f>& VACFoutput)
    {
      vtkm::cont::Invoker{}(ComputeVACFWorklet{}, original_velocity, current_velocity, VACFoutput);
    }

  namespace atoms
  {
    struct ComputeRDFWorklet : vtkm::worklet::WorkletMapField
    {
      ComputeRDFWorklet(const Id& num_center_pos, const Real& rho)
        : _num_center_pos(num_center_pos)
        , _rho(rho)
      {
      }
      using ControlSignature = void(FieldIn current_radius,
                                    WholeArrayIn center_pts,
                                    WholeArrayIn target_pts,
                                    ExecObject locator,
                                    FieldInOut rdf);
      using ExecutionSignature = void(_1, _2, _3, _4, _5);

      template<typename RadiusType, typename WholePtsType, typename RdfType>
      VTKM_EXEC void operator()(const RadiusType& current_radius,
                                const WholePtsType& center_pts,
                                const WholePtsType& target_pts,
                                const ExecPointLocator& locator,
                                RdfType& rdf) const
      {
        auto GetNumber = [=](const Vec3f& p_i) -> Id
        {
          Id number = 0;

          auto current_ijk = locator.InWhichCell(p_i);
          vtkm::Id3 dx = vtkm::Vec3f{ current_radius } / locator.GetDxdydz() + 1;
          auto ijk_upper = current_ijk + dx;
          auto ijk_lower = current_ijk - dx;

          for (Id i = ijk_lower[0]; i <= ijk_upper[0]; i++)
          {
            for (Id j = ijk_lower[1]; j <= ijk_upper[1]; j++)
            {
              for (Id k = ijk_lower[2]; k <= ijk_upper[2]; k++)
              {
                auto neighbor_ijk = Id3{ i, j, k };
                auto ijk = locator.PeriodicIndexOffset(neighbor_ijk);
                auto coord_offset = locator.PeriodicCoordOffset(ijk - neighbor_ijk);
                auto num_pts = locator.NumberPtsInCell(ijk);


                for (Id p = 0; p < num_pts; p++)
                {
                  auto pts_id = locator.PtsInCell(ijk, p);
                  auto p_j = target_pts.Get(pts_id) - coord_offset;
                  Real dis = vtkm::Magnitude(p_j - p_i);

                  if (dis > current_radius && dis < current_radius + current_radius * 0.01)
                  {
                    number++;
                  }
                }
              }
            }
          }

          return number;
        };

        auto total_num = 0;
        for (auto i = 0; i < _num_center_pos; ++i)
        {
          const auto& p_i = center_pts.Get(i);
          total_num += GetNumber(p_i);
        }
        auto dr3 = vtkm::Pow((current_radius + current_radius * 0.01), 3);
        auto r3 = vtkm::Pow(current_radius, 3);
        rdf += total_num /
          (_num_center_pos * _rho *
           ((4 / 3.0) * vtkm::Pif() * dr3 - ((4 / 3.0) * vtkm::Pif() * r3)));
      }

      Id _num_center_pos;
      Real _rho;
    };

    void ComputeRDF(const Id& num_center_pos,
                    const Real& _rho,
                    const vtkm::cont::ArrayHandle<Real>& radius,
                    const vtkm::cont::ArrayHandle<vtkm::Vec3f>& center_position,
                    const vtkm::cont::ArrayHandle<vtkm::Vec3f>& target_position,
                    const ContPointLocator& locator,
                    vtkm::cont::ArrayHandle<Real>& rdf)
    {
      vtkm::cont::Invoker{}(ComputeRDFWorklet{ num_center_pos, _rho },
                            radius,
                            center_position,
                            target_position,
                            locator,
                            rdf);
    }

  }
}