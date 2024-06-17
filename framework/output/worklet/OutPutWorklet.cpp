#include "output/worklet/OutPutWorklet.h"
#include "System.h"
#include "Executioner.h"
#include <vtkm/worklet/WorkletMapField.h>
#include <vtkm/cont/Algorithm.h>
#include "Application.h"
#include <fmt/format.h>
#include "math/Math.h"

namespace OutPut
{
    struct ComputeEAMrhoWorklet : vtkm::worklet::WorkletMapField
    {
      ComputeEAMrhoWorklet(const Real& eam_cut_off, const Real& Vlength)
        : _eam_cut_off(eam_cut_off)
        , _Vlength(Vlength)
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
          auto r_ij = locator.MinDistance(p_i, p_j, _Vlength);
          rho += force_function.ComputeEAMrhoOUT(_eam_cut_off, r_ij, rhor_spline);
        };
        locator.ExecuteOnNeighbor(atoms_id, function);
        EAM_rho = rho;
      }
      Real _eam_cut_off;
      Real _Vlength;
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
      ComputePairEnergyWorklet(const Real& eam_cut_off, const Real& Vlength)
        : _eam_cut_off(eam_cut_off)
        , _Vlength(Vlength)
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
          auto r_ij = locator.MinDistance(p_i, p_j, _Vlength);
          pair_energy += force_function.ComputePairEnergy(_eam_cut_off, r_ij, z2r_spline);
        };
        locator.ExecuteOnNeighbor(atoms_id, function);
      }
      Real _eam_cut_off;
      Real _Vlength;
    };

    struct ComputePotentialEnWorklet : vtkm::worklet::WorkletMapField
    {
      ComputePotentialEnWorklet(const Real& cut_off)
        : _cut_off(cut_off)
      {
      }

      using ControlSignature = void(FieldIn atoms_id,
                                    ExecObject locator,
                                    ExecObject topology,
                                    ExecObject force_function,
                                    FieldOut PotentialEnergy);
      using ExecutionSignature = void(_1, _2, _3,_4 ,_5);

      VTKM_EXEC void operator()(vtkm::Id atoms_id,
                                const ExecPointLocator& locator,
                                const ExecTopology& topology,
                                const ExecForceFunction& force_function,
                                Real& PotentialEnergy) const
      {
        Real PE_ij = 0;
        const auto& molecular_id_i = topology.GetMolecularId(atoms_id);
        const auto& pts_type_i = topology.GetAtomsType(atoms_id);
        auto eps_i = topology.GetEpsilon(pts_type_i);
        auto sigma_i = topology.GetSigma(pts_type_i);

        auto function = [&](const Vec3f& p_i, const Vec3f& p_j, const Id& pts_id_j)
        {
          auto pts_type_j = topology.GetAtomsType(pts_id_j);
          auto eps_j = topology.GetEpsilon(pts_type_j);
          auto sigma_j = topology.GetSigma(pts_type_j);
          auto r_ij = p_j - p_i;
          auto eps_ij = vtkm::Sqrt(eps_i * eps_j);
          auto sigma_ij = (sigma_i + sigma_j) / 2;

          // special lj part
          auto molecular_id_j = topology.GetMolecularId(pts_id_j);
          IdComponent force_factor_ij = (molecular_id_i == molecular_id_j) ? 0 : 1.0;

          PE_ij += force_factor_ij * force_function.ComputePotentialEn(r_ij, eps_ij, sigma_ij, _cut_off);
        };
        locator.ExecuteOnNeighbor(atoms_id, function);

        PotentialEnergy = PE_ij;
      }
      Real _cut_off;
    };

    struct ComputeSpecialBondsLJPotentialWorklet : vtkm::worklet::WorkletMapField
    {
      ComputeSpecialBondsLJPotentialWorklet(const Real& cut_off)
        : _cut_off(cut_off)
      {
      }

      using ControlSignature = void(FieldIn atoms_id,
                                    ExecObject locator,
                                    ExecObject topology,
                                    ExecObject force_function,
                                    FieldIn special_ids,
                                    FieldIn special_weights,
                                    FieldOut PotentialEnergy);
      using ExecutionSignature = void(_1, _2, _3, _4, _5, _6,_7);

      template<typename idsType, typename weightsType>
      VTKM_EXEC void operator()(vtkm::Id atoms_id,
                                const ExecPointLocator& locator,
                                const ExecTopology& topology,
                                const ExecForceFunction& force_function,
                                const idsType& special_ids,
                                const weightsType& special_weights,
                                Real& PotentialEnergy) const
      {
        Real PE_ij = 0;
        const auto& molecular_id_i = topology.GetMolecularId(atoms_id);
        const auto& pts_type_i = topology.GetAtomsType(atoms_id);
        auto eps_i = topology.GetEpsilon(pts_type_i);
        auto sigma_i = topology.GetSigma(pts_type_i);
        auto num_components = special_ids.GetNumberOfComponents();

        auto function = [&](const Vec3f& p_i, const Vec3f& p_j, const Id& pts_id_j)
        {
          auto pts_type_j = topology.GetAtomsType(pts_id_j);
          auto eps_j = topology.GetEpsilon(pts_type_j);
          auto sigma_j = topology.GetSigma(pts_type_j);
          auto r_ij = p_j - p_i;
          auto eps_ij = vtkm::Sqrt(eps_i * eps_j);
          auto sigma_ij = (sigma_i + sigma_j) / 2;

          auto pe = 
              force_function.ComputePotentialEn(r_ij, eps_ij, sigma_ij, _cut_off);

           // special lj part
          auto weight = 1.0;
          for (auto i = 0; i < num_components; ++i)
          {
            if (special_ids[i] == pts_id_j)
            {

              weight = special_weights[i];
            }
          }

          PE_ij += weight * pe;
        };
        locator.ExecuteOnNeighbor(atoms_id, function);

        PotentialEnergy = PE_ij;
      }
      Real _cut_off;
    };
    
    struct ComputeNearElePotentialWorklet : vtkm::worklet::WorkletMapField
    {
      ComputeNearElePotentialWorklet(const Real& cut_off, const Real& nearalpha)
        : _cut_off(cut_off)
        , _nearalpha(nearalpha)
      {
      }

      using ControlSignature = void(FieldIn atoms_id,
                                    ExecObject locator,
                                    ExecObject topology,
                                    ExecObject force_function,
                                    FieldOut PotentialEnergy); //pointer
      using ExecutionSignature = void(_1, _2, _3, _4, _5);

      VTKM_EXEC void operator()(const Id& atoms_id,
                                const ExecPointLocator& locator,
                                const ExecTopology& topology,
                                const ExecForceFunction& force_function,
                                Real& PotentialNearEleEnergy) const
      {
        Real PE_ij = 0;
        const auto& charge_i = topology.GetCharge(atoms_id);

        auto function = [&](const Vec3f& p_i, const Vec3f& p_j, const Id& pts_id_j)
        {
          auto charge_j = topology.GetCharge(pts_id_j);
          auto r_ij = p_j - p_i;

          PE_ij +=force_function.ComputeNearEleEnergy(r_ij, charge_i, charge_j, _cut_off, _nearalpha);
        };
        locator.ExecuteOnNeighbor(atoms_id, function);

        PotentialNearEleEnergy = PE_ij;
      }

      Real _cut_off; 
      Real _nearalpha;
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

    struct ComputeSpecialFarCoulWorklet : vtkm ::worklet::WorkletMapField
    {
        ComputeSpecialFarCoulWorklet(const Real& Vlength)
      : _Vlength(Vlength)
    {
    }
        using ControlSignature = void(FieldIn current_atoms_id,
                                      FieldIn group_vec,
                                      ExecObject locator,
                                      ExecObject topology,
                                      ExecObject forcefunction,
                                      FieldOut spec_far_coul_energy);
        using ExecutionSignature = void(_1, _2, _3,_4, _5, _6);
        template<typename GroupVecType>

        VTKM_EXEC void operator()(const Id& atoms_id,
                                  const GroupVecType& group_vec,
                                  const ExecPointLocator& locator,
                                  const ExecTopology& topology,
                                  const ExecForceFunction& force_function,
                                  Real& spec_far_coul_energy) const
        {
            auto num = group_vec.GetNumberOfComponents();
            spec_far_coul_energy = 0;

            if (num < 3)
            {
              spec_far_coul_energy = 0;
            }
            else
            {
              auto charge_p_i = topology.GetCharge(atoms_id);
              auto atoms_coord_i = locator.GetPtsPosition(atoms_id);
              for (Id j = 0; j < num; ++j)
              {
                auto id = group_vec[j];
                if (atoms_id == id)
                  continue;

                auto atoms_charge_j = topology.GetCharge(id);
                auto atoms_coord_j = locator.GetPtsPosition(id);
                Vec3f rij = locator.MinDistanceIf(atoms_coord_i, atoms_coord_j, _Vlength);
                Real dis_ij = vtkm::Magnitude(rij);
                Real energy_atom = charge_p_i * atoms_charge_j / dis_ij;
                spec_far_coul_energy += energy_atom;
              }
            }
        }
        Real _Vlength;
    };

    struct ComputeSpecialBondsCoulWorklet : vtkm ::worklet::WorkletMapField
    {
      ComputeSpecialBondsCoulWorklet(const Real& Vlength)
        : _Vlength(Vlength)
      {
      }
      using ControlSignature = void(FieldIn current_atoms_id,
                                    FieldIn group_vec,
                                    ExecObject locator,
                                    ExecObject topology,
                                    ExecObject forcefunction,
                                    FieldIn special_ids,
                                    FieldIn special_weights,
                                    FieldOut spec_far_coul_energy);
      using ExecutionSignature = void(_1, _2, _3, _4, _5, _6,_7,_8);

      template<typename GroupVecType, typename idsType, typename weightsType>
      VTKM_EXEC void operator()(const Id& atoms_id,
                                const GroupVecType& group_vec,
                                const ExecPointLocator& locator,
                                const ExecTopology& topology,
                                const ExecForceFunction& force_function,
                                const idsType& special_ids,
                                const weightsType& special_weights,
                                Real& spec_far_coul_energy) const
      {
        auto num = group_vec.GetNumberOfComponents();
        spec_far_coul_energy = 0;

        if (num < 3)
        {
          spec_far_coul_energy = 0;
        }
        else
        {
          auto charge_p_i = topology.GetCharge(atoms_id);
          auto atoms_coord_i = locator.GetPtsPosition(atoms_id);
          for (Id j = 0; j < num; ++j)
          {
            auto id = group_vec[j];
            if (atoms_id == id)
              continue;

            auto atoms_charge_j = topology.GetCharge(id);
            auto atoms_coord_j = locator.GetPtsPosition(id);
            Vec3f rij = locator.MinDistanceIf(atoms_coord_i, atoms_coord_j, _Vlength);
            Real dis_ij = vtkm::Magnitude(rij);
            Real energy_atom = charge_p_i * atoms_charge_j / dis_ij;
            //spec_far_coul_energy += energy_atom;

            auto num_components = special_ids.GetNumberOfComponents();
            Real weight = 1.0;
            for (auto i = 0; i < num_components; ++i)
            {
              if (special_ids[i] == id)
              {
                weight = special_weights[i];
              }
            }

            spec_far_coul_energy += (1.0 - weight) * energy_atom;
          }
        }
      }
      Real _Vlength;
    };
    
    struct ComputeSqChargeWorklet : vtkm::worklet::WorkletMapField
    {
      using ControlSignature = void(FieldIn charge, FieldOut sq_charge);
      using ExecutionSignature = void(_1, _2);
    
      VTKM_EXEC void operator()(const Real& charge, Real& sq_charge) const
      {
        sq_charge = charge * charge;
      }
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
      ComputeMSDWorklet(const Real& Vlength)
        : _Vlength(Vlength)
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
        temp_position = current_pts + position_flag * _Vlength;
        d_position = temp_position - original_pts;

        MSDoutput[0] += d_position[0] * d_position[0];
        MSDoutput[1] += d_position[1] * d_position[1];
        MSDoutput[2] += d_position[2] * d_position[2];
        MSDoutput[3] += d_position[0] * d_position[0] + d_position[1] * d_position[1] +
          d_position[2] * d_position[2];
      }

      Real _Vlength;
    };

    struct ComputePrepareMSDWorklet : vtkm::worklet::WorkletMapField
    {
      ComputePrepareMSDWorklet(const Real& Vlength)
        : _Vlength(Vlength)
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
        
        d_position = current_pts + position_flag * _Vlength - temp_MSD_pts;
       
        for (int i = 0; i < 3; i++)
        {
          if (vtkm::Abs(d_position[i]) < 1.0)
          {
            diff_MSD_pts[i] = d_position[i];
          }
          else
            diff_MSD_pts[i] = d_position[i] - position_flag[i] * _Vlength;
        }
        
        
      }

      Real _Vlength;
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
                 const Real& Vlength,
                 const vtkm::cont::ArrayHandle<vtkm::Id>& atoms_id,
                 const vtkm::cont::ArrayHandle<Vec7f>& rhor_spline,
                 const ContPointLocator& locator,
                 const ContTopology& topology,
                 const ContForceFunction& force_function,
                 vtkm::cont::ArrayHandle<Real>& EAM_rho)
    {
      vtkm::cont::Invoker{}(ComputeEAMrhoWorklet{ eam_cut_off, Vlength },
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
                        const Real& Vlength,
                        const vtkm::cont::ArrayHandle<vtkm::Id>& atoms_id,
                        const vtkm::cont::ArrayHandle<Vec7f>& z2r_spline,
                        const ContPointLocator& locator,
                        const ContTopology& topology,
                        const ContForceFunction& force_function,
                        vtkm::cont::ArrayHandle<Real>& pair_energy)
    {

      vtkm::cont::Invoker{}(ComputePairEnergyWorklet{ eam_cut_off, Vlength },
                            atoms_id,
                            z2r_spline,
                            locator,
                            topology,
                            force_function,
                            pair_energy);
    }

    //Statistical PotentialEnergy
    void ComputePotentialEnergy(const Real& cutoff,
                                const vtkm::cont::ArrayHandle<vtkm::Id>& atoms_id,
                                const ContPointLocator& locator,
                                const ContTopology& topology,
                                const ContForceFunction& force_function,
                                vtkm::cont::ArrayHandle<Real>& potential_energy)
    {
      vtkm::cont::Invoker{}(ComputePotentialEnWorklet{ cutoff },
                            atoms_id,
                            locator,
                            topology,
                            force_function,
                            potential_energy);
    }

    void ComputeSpecialBondsLJPotential(const Real& cutoff,
                                        const vtkm::cont::ArrayHandle<vtkm::Id>& atoms_id,
                                        const ContPointLocator& locator,
                                        const ContTopology& topology,
                                        const ContForceFunction& force_function,
                                        const GroupIdIdType& group_ids,
                                        const GroupRealIdType& group_weights,
                                        vtkm::cont::ArrayHandle<Real>& potential_energy)
    {
      vtkm::cont::Invoker{}(ComputeSpecialBondsLJPotentialWorklet{ cutoff },
                            atoms_id,
                            locator,
                            topology,
                            force_function,
                            group_ids,
                            group_weights,
                            potential_energy);
    }

    void ComputeNearElePotential(const Real& cutoff,
                                 const Real& alpha,
                                 const vtkm::cont::ArrayHandle<vtkm::Id>& atoms_id,
                                 const ContPointLocator& locator,
                                 const ContTopology& topology,
                                 const ContForceFunction& force_function,
                                 vtkm::cont::ArrayHandle<Real>& near_potential_energy)
    {
      vtkm::cont::Invoker{}(ComputeNearElePotentialWorklet{ cutoff, alpha },
                            atoms_id,
                            locator,
                            topology,
                            force_function,
                            near_potential_energy);
    }

    void ComputeDensity(const vtkm::Vec3f& K,
                        const vtkm::cont::ArrayHandle<vtkm::Vec3f>& position,
                        const vtkm::cont::ArrayHandle<Real>& charge,
                        vtkm::cont::ArrayHandle<Real>& DensityReal,
                        vtkm::cont::ArrayHandle<Real>& DensityImage)
    {
      vtkm::cont::Invoker{}(ComputeDensityWorklet{ K }, position, charge, DensityReal, DensityImage);

    }

    using GroupVecType = vtkm::cont::ArrayHandleGroupVecVariable<vtkm::cont::ArrayHandle<vtkm::Id>, 
                                                                vtkm::cont::ArrayHandle<vtkm::Id>>;
    void ComputeSpecialFarCoul(const Real& Vlength,
                               const vtkm::cont::ArrayHandle<Id>& atoms_id,
                               const GroupVecType& group_vec,
                               const ContPointLocator& locator,
                               const ContTopology& topology,
                               const ContForceFunction& force_function,
                               vtkm::cont::ArrayHandle<Real>& SpecFarEnergy)
    {
      vtkm::cont::Invoker{}(ComputeSpecialFarCoulWorklet{ Vlength }, atoms_id, group_vec, locator, topology, force_function, SpecFarEnergy);
    }

    void ComputeSpecialBondsCoul(const Real& Vlength,
                               const vtkm::cont::ArrayHandle<Id>& atoms_id,
                               const GroupVecType& group_vec,
                               const ContPointLocator& locator,
                               const ContTopology& topology,
                               const ContForceFunction& force_function,
                               const GroupIdIdType& group_ids,
                               const GroupRealIdType& group_weights,
                               vtkm::cont::ArrayHandle<Real>& SpecFarEnergy)
    {
      vtkm::cont::Invoker{}(ComputeSpecialBondsCoulWorklet{ Vlength },
                            atoms_id,
                            group_vec,
                            locator,
                            topology,
                            force_function,
                            group_ids,
                            group_weights,
                            SpecFarEnergy);
    }

    void ComputeSqCharge(const vtkm::cont::ArrayHandle<Real>& charge,
                         vtkm::cont::ArrayHandle<Real>& SelfEnergy)
    {
      vtkm::cont::Invoker{}(ComputeSqChargeWorklet{}, charge, SelfEnergy);
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

     void GetPositionByType(const vtkm::cont::ArrayHandle<Id>& atom_id,
                           const vtkm::cont::ArrayHandle<Vec3f>& whole_pts,
                           vtkm::cont::ArrayHandle<Vec3f>& pts_by_type)
    {
       vtkm::cont::Invoker{}(GetPositionByTypeWorklet{}, atom_id, whole_pts, pts_by_type);
    }

    // Statistical Kinetic energy
    void ComputerKineticEnergy(const vtkm::cont::ArrayHandle<vtkm::Vec3f>& velocity,
                               const vtkm::cont::ArrayHandle<Real>& mass,
                               vtkm::cont::ArrayHandle<Real>& sq_velocity)
    {
      vtkm::cont::Invoker{}(ComputerKineticEnergyWorklet{}, velocity, mass, sq_velocity);
    }

    // Statistical MSD
    void ComputeMSD(const Real& _Vlength,
                    const vtkm::cont::ArrayHandle<vtkm::Vec3f>& original_position,
                    const vtkm::cont::ArrayHandle<vtkm::Vec3f>& current_pts_position,
                    const vtkm::cont::ArrayHandle<vtkm::Id3>& position_flag,
                    vtkm::cont::ArrayHandle<vtkm::Vec4f>& MSDoutput) 
    {
      vtkm::cont::Invoker{}(ComputeMSDWorklet{ _Vlength },
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