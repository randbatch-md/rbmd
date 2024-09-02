#include "DataObject.h"
#include "worklet/RunWorklet.h"
#include "math/Math.h"
#include <vtkm/worklet/WorkletMapField.h>
#include <vtkm/worklet/WorkletMapTopology.h>
#include "locator/ContPointLocator.h"
#include "locator/ExecPointLocator.h"
#include "forceFunction/ContForceFunction.h"
#include "topology/ContTopology.h"
#include "math/Utiles.h"

namespace RunWorklet
{
    struct ComputeLJVirialWorklet : vtkm::worklet::WorkletMapField
    {
  ComputeLJVirialWorklet(const Real& cut_off, const Vec3f& box)
    : _cut_off(cut_off)
    , _box(box)
  {
  }

     using ControlSignature = void(FieldIn atoms_id,
                                ExecObject locator,
                                ExecObject topology,
                                ExecObject force_function,
                                FieldIn group_j,
                                FieldIn num_j,
                                FieldIn coord_offset_j,
                                FieldOut LJVirial);
  using ExecutionSignature = void(_1, _2, _3, _4, _5, _6, _7, _8);

  template<typename NeighbourGroupVecType, typename CoordOffsetj>
  VTKM_EXEC void operator()(const Id atoms_id,
                            const ExecPointLocator& locator,
                            const ExecTopology& topology,
                            const ExecForceFunction& force_function,
                            const NeighbourGroupVecType& group_j,
                            const Id& num_j,
                            const CoordOffsetj& coord_offset_j,
                            Vec6f& LJVirial) const
  {
    Vec6f lj_v = { 0, 0, 0, 0, 0, 0 };

    const auto& molecular_id_i = topology.GetMolecularId(atoms_id);
    const auto& pts_type_i = topology.GetAtomsType(atoms_id);
    auto eps_i = topology.GetEpsilon(pts_type_i);
    auto sigma_i = topology.GetSigma(pts_type_i);

    auto function = [&](const Vec3f& p_i, const Vec3f& p_j, const Id& pts_id_j)
    {
      auto molecular_id_j = topology.GetMolecularId(pts_id_j);
      auto pts_type_j = topology.GetAtomsType(pts_id_j);
      auto eps_j = topology.GetEpsilon(pts_type_j);
      auto sigma_j = topology.GetSigma(pts_type_j);
      //auto r_ij = p_j - p_i;
      auto r_ij = locator.MinDistanceVec(p_j, p_i, _box);

      if (molecular_id_i == molecular_id_j)
        return;
      lj_v += force_function.ComputeLJVirial(r_ij, eps_i, eps_j, sigma_i, sigma_j, _cut_off);
    };

    auto p_i = locator.GetPtsPosition(atoms_id);

    for (Id p = 0; p < num_j; p++)
    {
      auto idj = group_j[p];
      auto p_j = locator.GetPtsPosition(idj) - coord_offset_j[p];
      function(p_i, p_j, idj);
    }

    LJVirial = lj_v;
  }
  Real _cut_off;
  Vec3f _box;
  };

    struct ComputeCoulVirialWorklet : vtkm::worklet::WorkletMapField
  {
  ComputeCoulVirialWorklet(const Real& cut_off, const Vec3f& box)
    : _cut_off(cut_off)
    , _box(box)
  {
  }

  using ControlSignature = void(FieldIn atoms_id,
                                ExecObject locator,
                                ExecObject topology,
                                ExecObject force_function,
                                FieldIn group_j,
                                FieldIn num_j,
                                FieldIn coord_offset_j,
                                FieldOut CoulVirial);
  using ExecutionSignature = void(_1, _2, _3, _4, _5, _6, _7, _8);

  template<typename NeighbourGroupVecType, typename CoordOffsetj>
  VTKM_EXEC void operator()(const Id atoms_id,
                            const ExecPointLocator& locator,
                            const ExecTopology& topology,
                            const ExecForceFunction& force_function,
                            const NeighbourGroupVecType& group_j,
                            const Id& num_j,
                            const CoordOffsetj& coord_offset_j,
                            Vec6f& CoulVirial) const
  {
    Vec6f coul_v = { 0, 0, 0, 0, 0, 0 };

    auto charge_pi = topology.GetCharge(atoms_id);
    auto function = [&](const Vec3f& p_i, const Vec3f& p_j, const Id& pts_id_j)
    {
      auto charge_pj = topology.GetCharge(pts_id_j);

      //auto r_ij = p_j - p_i;
      auto r_ij = locator.MinDistanceVec(p_j, p_i, _box);

      coul_v += force_function.ComputeCoulVirial(r_ij, charge_pi, charge_pj,_cut_off);

    };

    auto p_i = locator.GetPtsPosition(atoms_id);

    for (Id p = 0; p < num_j; p++)
    {
      auto idj = group_j[p];
      auto p_j = locator.GetPtsPosition(idj) - coord_offset_j[p];
      function(p_i, p_j, idj);
    }
    CoulVirial = coul_v;
  }
  Real _cut_off;
  Vec3f _box;
  };

    struct ComputeEwaldVirialWorklet : vtkm ::worklet::WorkletMapField
  {
  ComputeEwaldVirialWorklet(const IdComponent& Kmax, const Real& unit_factor)
    : _Kmax(Kmax)
    , _unit_factor(unit_factor)
  {
  }

  using ControlSignature = void(FieldIn atoms_id,
                                WholeArrayIn whole_rhok,
                                ExecObject force_function,
                                ExecObject topology,
                                ExecObject locator,
                                FieldOut EwaldVirial);
  using ExecutionSignature = void(_1, _2, _3, _4, _5, _6);

  template<typename WholeRhokType>
  VTKM_EXEC void operator()(const Id atoms_id,
                            const WholeRhokType& whole_rhok,
                            ExecForceFunction& force_function,
                            const ExecTopology& topology,
                            const ExecPointLocator& locator,
                            Vec6f& EwaldVirial) const
  {
    Vec6f  Ewald_v = { 0, 0, 0, 0, 0, 0 };
    const auto& charge_p_i = topology.GetCharge(atoms_id);
    const auto& p_i = locator.GetPtsPosition(atoms_id);
    auto function = [&](const vtkm::Vec3f& M, const vtkm::Id& indexEwald)
    {
      auto rhok_ri = whole_rhok.Get(indexEwald - 1);
      Ewald_v += force_function.ComputeLongVirial(M, p_i, charge_p_i, rhok_ri);
    };

    locator.ExecuteOnKNeighbor(_Kmax, atoms_id, function);
    EwaldVirial = Ewald_v;
  }
  IdComponent _Kmax;
  Real _unit_factor;
  };   

    struct fix_press_berendsenWorklet : vtkm::worklet::WorkletMapField
    {
    fix_press_berendsenWorklet(const Real& scale_factor)
    : _scale_factor(scale_factor)
    {
     }

     using ControlSignature = void(FieldInOut position, ExecObject locator);
     using ExecutionSignature = void(_1, _2);

     VTKM_EXEC void operator()(Vec3f& position, const ExecPointLocator locator) const
    {
      position = _scale_factor * position;
      //locator.UpdateOverRangePoint(position);
    }
     Real _scale_factor;
     };

    struct ApplyPbcWorklet : vtkm::worklet::WorkletMapField
     {
     ApplyPbcWorklet(const Vec3f& box, const Vec<Vec2f, 3>& range)
       : _box(box)
       , _range(range)
     {
     }

     using ControlSignature = void(FieldInOut position, ExecObject locator);
     using ExecutionSignature = void(_1, _2);

     template<typename CoordType>
     VTKM_EXEC void operator()(CoordType& position,
                               const ExecPointLocator locator) const
     {
      for (vtkm::Id i = 0; i < 3; ++i)
      {
        if (position[i] < _range[i][0])
        {
          position[i] += _box[i];
        }
        else if (position[i] >= _range[i][1])
        {
          position[i] -= _box[i];
        }
      }
     }
     Vec3f _box;
     Vec<Vec2f, 3> _range;
     };

    struct ApplyPbcFlagWorklet : vtkm::worklet::WorkletMapField
     {
     ApplyPbcFlagWorklet(const Vec3f& box, const Vec<Vec2f, 3>& range)
       : _box(box)
       , _range(range)
     {
     }

     using ControlSignature = void(FieldInOut position,
                                   ExecObject locator,
                                   FieldInOut position_flag);
     using ExecutionSignature = void(_1, _2, _3);

     template<typename CoordType>
     VTKM_EXEC void operator()(CoordType& position,
                               const ExecPointLocator locator,
                               Id3& position_flag) const
     {
      for (vtkm::Id i = 0; i < 3; ++i)
      {
        if (position[i] < _range[i][0])
        {
          position[i] += _box[i];
          position_flag[i] -= 1;
        }
        else if (position[i] >= _range[i][1])
        {
          position[i] -= _box[i];
          position_flag[i] += 1;
        }
      }
     }
     Vec3f _box;
     Vec<Vec2f, 3> _range;
     };

    struct X2LamdaWorklet : vtkm::worklet::WorkletMapField
    {                   
        X2LamdaWorklet(const Vec6f& h_inv,
                       const vtkm::Vec<vtkm::Range, 3>& range)
        : _h_inv(h_inv)
        , _range(range)
        {
        }
        
        using ControlSignature = void(FieldInOut position);
        using ExecutionSignature = void(_1);
        
        template<typename CoordType>
        VTKM_EXEC void operator()(CoordType& position) const
        {
         Vec3f delta;
         delta[0] = position[0] - _range[0].Min;
         delta[1] = position[1] - _range[1].Min;
         delta[2] = position[2] - _range[2].Min;
        
         // 计算 lamda 坐标
         Vec3f lamda_position;
         lamda_position[0] = _h_inv[0] * delta[0] + _h_inv[5] * delta[1] + _h_inv[4] * delta[2];
         lamda_position[1] = _h_inv[1] * delta[1] + _h_inv[3] * delta[2];
         lamda_position[2] = _h_inv[2] * delta[2];
        
         // 更新位置
         position = lamda_position;
        }
        Vec6f _h_inv;
        vtkm::Vec<vtkm::Range, 3> _range;  
    };

    struct Lamda2XWorklet : vtkm::worklet::WorkletMapField
    {       
        Lamda2XWorklet(const Vec6f& h, const vtkm::Vec<vtkm::Range, 3>& range)
          : _h(h)
          , _range(range)
        {
        }
        using ControlSignature = void(FieldInOut position);
        using ExecutionSignature = void(_1);

        template<typename CoordType>
        VTKM_EXEC void operator()(CoordType& position) const
        {
         Vec3f position_base;
         position_base[0] = position[0];
         position_base[1] = position[1];
         position_base[2] = position[2];

         // compute  real position
         Vec3f x_position;
         x_position[0] = _h[0] * position_base[0] + _h[5] * position_base[1] +
           _h[4] * position_base[2] + static_cast<vtkm::FloatDefault>(_range[0].Min);

         x_position[1] = _h[1] * position_base[1] + _h[3] * position_base[2] +
           static_cast<vtkm::FloatDefault>(_range[1].Min);

         x_position[2] = _h[2] * position_base[2] + static_cast<vtkm::FloatDefault>(_range[2].Min);
        
         position = x_position;
        }
        Vec6f _h;                         
        vtkm::Vec<vtkm::Range, 3> _range; 
    };

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
                if ( _cut_off_2 - dis_2 > 0.0001 && dis_2 > 0.0001)
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

    struct ComputeNearForceVerletERFWorklet : vtkm::worklet::WorkletMapField
    {
      ComputeNearForceVerletERFWorklet(const Real& cut_off)
        : _cut_off(cut_off)
      {
      }

      using ControlSignature = void(FieldIn atoms_id,
                                    ExecObject locator,
                                    ExecObject topology,
                                    ExecObject force_function,
                                    ExecObject static_table,
                                    FieldIn group_j,
                                    FieldIn num_j,
                                    FieldIn coord_offset_j,
                                    FieldOut LJforce);
      using ExecutionSignature = void(_1, _2, _3, _4, _5, _6, _7,_8,_9);

      template<typename NeighbourGroupVecType, typename CoordOffsetj>
      VTKM_EXEC void operator()(const Id atoms_id,
                                const ExecPointLocator& locator,
                                const ExecTopology& topology,
                                const ExecForceFunction& force_function,
                                const ExecStaticTable& static_table,
                                const NeighbourGroupVecType& group_j,
                                const Id& num_j,
                                const CoordOffsetj& coord_offset_j,
                                Vec3f& LJforce) const
      {
        Vec3f LJ_force = { 0, 0, 0 };
        Vec3f ele_force = { 0, 0, 0 };

        LJforce = { 0, 0, 0 };
        const auto& molecular_id_i = topology.GetMolecularId(atoms_id);
        const auto& pts_type_i = topology.GetAtomsType(atoms_id);
        auto eps_i = topology.GetEpsilon(pts_type_i);
        auto sigma_i = topology.GetSigma(pts_type_i);
        auto charge_pi = topology.GetCharge(atoms_id);

        auto function = [&](const Vec3f& p_i, const Vec3f& p_j, const Id& pts_id_j)
        {
          auto charge_pj = topology.GetCharge(pts_id_j);
          auto molecular_id_j = topology.GetMolecularId(pts_id_j);
          auto pts_type_j = topology.GetAtomsType(pts_id_j);
          auto eps_j = topology.GetEpsilon(pts_type_j);
          auto sigma_j = topology.GetSigma(pts_type_j);
          auto r_ij = p_j - p_i;

          Real dispij = vtkm::Magnitude(r_ij);
          IdComponent index_table_pij = static_table.Extract(dispij);
          auto table_pij = static_table.TableGnearValue(dispij, index_table_pij);

          ele_force += force_function.ComputeNearEnergyForceERF(p_i, p_j, charge_pi, charge_pj, table_pij);

          if (molecular_id_i == molecular_id_j)
            return;
          LJ_force += force_function.ComputeLJForce(r_ij, eps_i, eps_j, sigma_i, sigma_j, _cut_off);
        };
        
        auto p_i = locator.GetPtsPosition(atoms_id);
        
        for (Id p = 0; p < num_j; p++)
        {
          auto idj = group_j[p];
          auto p_j = locator.GetPtsPosition(idj)-coord_offset_j[p];
          function(p_i, p_j, idj);
        }
        LJforce = LJ_force + ele_force;
      }
      Real _cut_off;
    };

    struct ComputeNearForceVerletWorklet : vtkm::worklet::WorkletMapField
    {
      ComputeNearForceVerletWorklet(const Real& cut_off, const Vec3f& box)
        : _cut_off(cut_off)
        , _box(box)
      {
      }

      using ControlSignature = void(FieldIn atoms_id,
                                    ExecObject locator,
                                    ExecObject topology,
                                    ExecObject force_function,
                                    FieldIn group_j,
                                    FieldIn num_j,
                                    FieldIn coord_offset_j,
                                    FieldIn special_ids,
                                    FieldIn special_weights,
                                    FieldOut nearforce);
      using ExecutionSignature = void(_1, _2, _3, _4, _5, _6, _7, _8,_9,_10);

      template<typename NeighbourGroupVecType,
               typename CoordOffsetj,
               typename idsType,
               typename weightsType>
      VTKM_EXEC void operator()(const Id atoms_id,
                                const ExecPointLocator& locator,
                                const ExecTopology& topology,
                                const ExecForceFunction& force_function,
                                const NeighbourGroupVecType& group_j,
                                const Id& num_j,
                                const CoordOffsetj& coord_offset_j,
                                const idsType& special_ids,
                                const weightsType& special_weights,
                                Vec3f& nearforce) const
      {
        Vec3f LJ_force = { 0, 0, 0 };
        Vec3f ele_force = { 0, 0, 0 };

        nearforce = { 0, 0, 0 };
        const auto& molecular_id_i = topology.GetMolecularId(atoms_id);
        const auto& pts_type_i = topology.GetAtomsType(atoms_id);
        auto eps_i = topology.GetEpsilon(pts_type_i);
        auto sigma_i = topology.GetSigma(pts_type_i);
        auto charge_pi = topology.GetCharge(atoms_id);
        auto num_components = special_ids.GetNumberOfComponents();

        auto function = [&](const Vec3f& p_i, const Vec3f& p_j, const Id& pts_id_j)
        {
          auto charge_pj = topology.GetCharge(pts_id_j);
          auto molecular_id_j = topology.GetMolecularId(pts_id_j);
          auto pts_type_j = topology.GetAtomsType(pts_id_j);
          auto eps_j = topology.GetEpsilon(pts_type_j);
          auto sigma_j = topology.GetSigma(pts_type_j);
          //auto r_ij = p_j - p_i;
          auto r_ij = locator.MinDistanceVec(p_j, p_i, _box);

          ele_force += force_function.ComputeNearEnergyForce1(r_ij, charge_pi, charge_pj, _cut_off);

          //if (molecular_id_i == molecular_id_j)
           // return;
          //LJ_force += force_function.ComputeLJForce(r_ij, eps_i, eps_j, sigma_i, sigma_j, _cut_off);
          
          //weight
          Vec3f forceij = force_function.ComputeLJForce(r_ij, eps_i, eps_j, sigma_i, sigma_j, _cut_off);
          auto weight = 1.0;
          for (auto i = 0; i < num_components; ++i)
          {
            if (special_ids[i] == pts_id_j)
            {
              weight = special_weights[i];
            }
          }
          LJ_force += weight * forceij;
        };

        auto p_i = locator.GetPtsPosition(atoms_id);

        for (Id p = 0; p < num_j; p++)
        {
          auto idj = group_j[p];
          auto p_j = locator.GetPtsPosition(idj) - coord_offset_j[p];
          function(p_i, p_j, idj);
        }
        nearforce = LJ_force + ele_force;
      }
      Real _cut_off;
      Vec3f _box;
    };

    struct ComputeLJForceVerletWorklet : vtkm::worklet::WorkletMapField
    {
      ComputeLJForceVerletWorklet(const Real& cut_off, const Vec3f& box)
        : _cut_off(cut_off)
        , _box(box)
      {
      }

      using ControlSignature = void(FieldIn atoms_id,
                                    ExecObject locator,
                                    ExecObject topology,
                                    ExecObject force_function,
                                    FieldIn group_j,
                                    FieldIn num_j,
                                    FieldIn coord_offset_j,
                                    FieldOut LJforce);
      using ExecutionSignature = void(_1, _2, _3, _4, _5, _6, _7, _8);

      template<typename NeighbourGroupVecType, typename CoordOffsetj>
      VTKM_EXEC void operator()(const Id atoms_id,
                                const ExecPointLocator& locator,
                                const ExecTopology& topology,
                                const ExecForceFunction& force_function,
                                const NeighbourGroupVecType& group_j,
                                const Id& num_j,
                                const CoordOffsetj& coord_offset_j,
                                Vec3f& LJforce) const
      {
        Vec3f LJ_force = { 0, 0, 0 };

        LJforce = { 0, 0, 0 };
        const auto& molecular_id_i = topology.GetMolecularId(atoms_id);
        const auto& pts_type_i = topology.GetAtomsType(atoms_id);
        auto eps_i = topology.GetEpsilon(pts_type_i);
        auto sigma_i = topology.GetSigma(pts_type_i);
        auto charge_pi = topology.GetCharge(atoms_id);

        auto function = [&](const Vec3f& p_i, const Vec3f& p_j, const Id& pts_id_j)
        {
          auto charge_pj = topology.GetCharge(pts_id_j);
          auto molecular_id_j = topology.GetMolecularId(pts_id_j);
          auto pts_type_j = topology.GetAtomsType(pts_id_j);
          auto eps_j = topology.GetEpsilon(pts_type_j);
          auto sigma_j = topology.GetSigma(pts_type_j);
          //auto r_ij = p_j - p_i;
          auto r_ij = locator.MinDistanceVec(p_j, p_i, _box);

          if (molecular_id_i == molecular_id_j)
            return;
          LJ_force += force_function.ComputeLJForce(r_ij, eps_i, eps_j, sigma_i, sigma_j, _cut_off);
        };

        auto p_i = locator.GetPtsPosition(atoms_id);

        for (Id p = 0; p < num_j; p++)
        {
          auto idj = group_j[p];
          auto p_j = locator.GetPtsPosition(idj) - coord_offset_j[p];
          function(p_i, p_j, idj);
        }
        LJforce = LJ_force;
      }
      Real _cut_off;
      Vec3f _box;
    };

    struct ComputeThermoWorklet : vtkm::worklet::WorkletMapField
    {
      ComputeThermoWorklet(const Real& cut_off, const Vec3f& box)
        : _cut_off(cut_off)
        , _box(box)
      {
      }

      using ControlSignature = void(FieldIn atoms_id,
                                    ExecObject locator,
                                    ExecObject topology,
                                    ExecObject force_function,
                                    FieldIn group_j,
                                    FieldIn num_j,
                                    FieldIn coord_offset_j,
                                    FieldOut LJforce,
                                    FieldOut LJvirial,
                                    FieldOut LJPE);
      using ExecutionSignature = void(_1, _2, _3, _4, _5, _6, _7, _8, _9, _10);

      template<typename NeighbourGroupVecType, typename CoordOffsetj>
      VTKM_EXEC void operator()(const Id atoms_id,
                                const ExecPointLocator& locator,
                                const ExecTopology& topology,
                                const ExecForceFunction& force_function,
                                const NeighbourGroupVecType& group_j,
                                const Id& num_j,
                                const CoordOffsetj& coord_offset_j,
                                Vec3f& LJforce,
                                Vec6f& LJvirial,
                                Real& LJPE) const
      {
        Vec3f LJ_force = { 0, 0, 0 };
        Vec6f virial = { 0, 0, 0, 0, 0, 0 };
        Real LJ_PE = 0;

        LJforce = { 0, 0, 0 };
        const auto& molecular_id_i = topology.GetMolecularId(atoms_id);
        const auto& pts_type_i = topology.GetAtomsType(atoms_id);
        auto eps_i = topology.GetEpsilon(pts_type_i);
        auto sigma_i = topology.GetSigma(pts_type_i);
        auto charge_pi = topology.GetCharge(atoms_id);

        auto function = [&](const Vec3f& p_i, const Vec3f& p_j, const Id& pts_id_j)
        {
          auto charge_pj = topology.GetCharge(pts_id_j);
          auto molecular_id_j = topology.GetMolecularId(pts_id_j);
          auto pts_type_j = topology.GetAtomsType(pts_id_j);
          auto eps_j = topology.GetEpsilon(pts_type_j);
          auto sigma_j = topology.GetSigma(pts_type_j);
          //auto r_ij = p_j - p_i;
          auto r_ij = locator.MinDistanceVec(p_j, p_i, _box);

          if (molecular_id_i == molecular_id_j)
            return;
          LJ_force += force_function.ComputeLJForce(r_ij, eps_i, eps_j, sigma_i, sigma_j, _cut_off);
          virial += force_function.ComputeLJVirial(r_ij, eps_i, eps_j, sigma_i, sigma_j, _cut_off);
          LJ_PE +=
            force_function.ComputePotentialEn0(r_ij, eps_i, eps_j, sigma_i, sigma_j, _cut_off);
        };

        auto p_i = locator.GetPtsPosition(atoms_id);

        for (Id p = 0; p < num_j; p++)
        {
          auto idj = group_j[p];
          auto p_j = locator.GetPtsPosition(idj) - coord_offset_j[p];
          function(p_i, p_j, idj);
        }
        LJforce = LJ_force;
        LJvirial = virial;
        LJPE = LJ_PE;
      }
      Real _cut_off;
      Vec3f _box;
    };

    struct ComputeEAMfpVerletWorklet : vtkm::worklet::WorkletMapField
    {
      ComputeEAMfpVerletWorklet(const Real& rc, const Vec3f& box)
        : _rc(rc)
        , _box(box)
      {
      }

      using ControlSignature = void(FieldIn atoms_id,
                                    WholeArrayIn rhor_spline,
                                    WholeArrayIn frho_spline,
                                    ExecObject locator,
                                    ExecObject force_function,
                                    FieldIn group_j,
                                    FieldIn num_j,
                                    FieldIn coord_offset_j,
                                    FieldOut EAM_fp);
      using ExecutionSignature = void(_1, _2, _3, _4, _5, _6, _7, _8,_9);

      template<typename NeighbourGroupVecType,
               typename CoordOffsetj,
               typename SpileType,
               typename frhoType>
      VTKM_EXEC void operator()(const Id atoms_id,
                                const SpileType& rhor_spline,
                                const frhoType& frho_spline,
                                const ExecPointLocator& locator,
                                const ExecForceFunction& force_function,
                                const NeighbourGroupVecType& group_j,
                                const Id& num_j,
                                const CoordOffsetj& coord_offset_j,
                                Real& EAM_fp) const
      {
        Real rho = 0;
        auto function = [&](const Vec3f& p_i, const Vec3f& p_j, const Id& pts_id_j)
        {
          auto r_ij = locator.MinDistanceVec(p_i, p_j, _box);
          rho += force_function.ComputeEAMrho(_rc, r_ij, rhor_spline);
        };

        auto p_i = locator.GetPtsPosition(atoms_id);

        for (Id i = 0; i < num_j; i++)
        {
          auto id_j = group_j[i];
          auto p_j = locator.GetPtsPosition(id_j) - coord_offset_j[i];
          function(p_i, p_j, id_j);
        }

        //EAM_rho = rho;
        EAM_fp = force_function.ComputeEAMfp(atoms_id, rho, frho_spline);
      }
      Real _rc;
      Vec3f _box;
    };

    struct ComputeEAMForceVerletWorklet : vtkm::worklet::WorkletMapField
    {
      ComputeEAMForceVerletWorklet(const Real& cut_off, const Vec3f& box)
        : _rc(cut_off)
        , _box(box)
      {
      }

      using ControlSignature = void(FieldIn atoms_id,
                                    WholeArrayIn rhor_spline,
                                    WholeArrayIn z2r_spline,
                                    WholeArrayIn EAM_fp,
                                    ExecObject locator,
                                    ExecObject force_function,
                                    FieldIn group_j,
                                    FieldIn num_j,
                                    FieldIn coord_offset_j,
                                    FieldOut force);
      using ExecutionSignature = void(_1, _2, _3, _4, _5, _6, _7, _8,_9,_10);

      template<typename NeighbourGroupVecType,
               typename CoordOffsetj,
               typename SpileType,
               typename z2rSpileType,
               typename EAM_fpType>
      VTKM_EXEC void operator()(const Id atoms_id,
                                const SpileType& rhor_spline,
                                const z2rSpileType& z2r_spline,
                                const EAM_fpType& EAM_fp,
                                const ExecPointLocator& locator,
                                const ExecForceFunction& force_function,
                                const NeighbourGroupVecType& group_j,
                                const Id& num_j,
                                const CoordOffsetj& coord_offset_j,
                                Vec3f& force) const
      {
        Vec3f rc_force = { 0, 0, 0 };

        auto function_force =
          [&](const Vec3f& p_i, const Vec3f& p_j, const Id& pts_id_j)
        {
          auto r_ij = locator.MinDistanceVec(p_j, p_i, _box);
          rc_force += force_function.ComputeEAMforce(_rc, atoms_id, pts_id_j, r_ij, EAM_fp, rhor_spline, z2r_spline);
        };

        auto p_i = locator.GetPtsPosition(atoms_id);
        for (Id p = 0; p < num_j; p++)
        {
          Id id_j = group_j[p];
          auto p_j = locator.GetPtsPosition(id_j) - coord_offset_j[p];
          function_force(p_i, p_j, id_j);
        }
        force = rc_force;
      }
      Real _rc;
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
          rho += force_function.ComputeEAMrho(_eam_cut_off, r_ij, rhor_spline);
        };
        locator.ExecuteOnNeighbor(atoms_id, function);
        EAM_rho = rho;
      }
      Real _eam_cut_off;
      Vec3f _box;
    };

    struct ComputeEAMfpWorklet : vtkm::worklet::WorkletMapField
    {
      ComputeEAMfpWorklet() {}
      using ControlSignature = void(FieldIn atoms_id,
                                    WholeArrayIn EAM_rho,
                                    WholeArrayIn frho_spline,
                                    ExecObject locator,
                                    ExecObject topology,
                                    ExecObject force_function,
                                    FieldOut fp);
      using ExecutionSignature = void(_1, _2, _3, _4, _5, _6, _7);

      template<typename EAM_rhoype, typename frho_splineType>
      VTKM_EXEC void operator()(const Id atoms_id,
                                const EAM_rhoype& EAM_rho,
                                const frho_splineType& frho_spline,
                                const ExecPointLocator& locator,
                                const ExecTopology& topology,
                                const ExecForceFunction& force_function,
                                Real& fp) const
      {
        fp = force_function.ComputeEAMfpOriginal(atoms_id, EAM_rho, frho_spline);
      }
    };

    struct ComputeEAMforceWorklet : vtkm::worklet::WorkletMapField
    {
      ComputeEAMforceWorklet(const Real& eam_cut_off, const Vec3f& box)
        : _eam_cut_off(eam_cut_off)
        , _box(box)
      {
      }

      using ControlSignature = void(FieldIn atoms_id,
                                    WholeArrayIn fp,
                                    WholeArrayIn rhor_spline,
                                    WholeArrayIn z2r_spline,
                                    ExecObject locator,
                                    ExecObject topology,
                                    ExecObject force_function,
                                    FieldOut eam_force);
      using ExecutionSignature = void(_1, _2, _3, _4, _5, _6, _7, _8);

      template<typename fpType, typename rhor_splineType, typename z2rSpileType>
      VTKM_EXEC void operator()(const Id atoms_id,
                                const fpType& fp,
                                const rhor_splineType& rhor_spline,
                                const z2rSpileType& z2r_spline,
                                const ExecPointLocator& locator,
                                const ExecTopology& topology,
                                const ExecForceFunction& force_function,
                                Vec3f& eam_force) const
      {
        Vec3f force = { 0, 0, 0 };
        auto function = [&](const Vec3f& p_i, const Vec3f& p_j, const Id& pts_id_j)
        {
          //auto r_ij = p_j - p_i;
          auto r_ij = locator.MinDistanceVec(p_i, p_j, _box);
          force += force_function.ComputeEAMforce(
            _eam_cut_off, atoms_id, pts_id_j, r_ij, fp, rhor_spline, z2r_spline);
        };
        locator.ExecuteOnNeighbor(atoms_id, function);
        eam_force = force;
      }
      Real _eam_cut_off;
      Vec3f _box;
    };

    struct ComputeLJForceWithPeriodicBCWorklet : vtkm::worklet::WorkletMapField
    {
      ComputeLJForceWithPeriodicBCWorklet(const Real& cut_off)
        : _cut_off(cut_off)
      {
      }

      using ControlSignature = void(FieldIn atoms_id,
                                    ExecObject locator,
                                    ExecObject topology,
                                    ExecObject force_function,
                                    FieldOut LJforce);
      using ExecutionSignature = void(_1, _2, _3, _4, _5);

      VTKM_EXEC void operator()(const Id atoms_id,
                                const ExecPointLocator& locator,
                                const ExecTopology& topology,
                                const ExecForceFunction& force_function,
                                Vec3f& LJforce) const
      {
        LJforce = { 0, 0, 0 };
        const auto& molecular_id_i = topology.GetMolecularId(atoms_id);
        const auto& pts_type_i = topology.GetAtomsType(atoms_id);
        auto eps_i = topology.GetEpsilon(pts_type_i);
        auto sigma_i = topology.GetSigma(pts_type_i);

        auto function = [&](const Vec3f& p_i, const Vec3f& p_j, const Id& pts_id_j)
        {
          auto molecular_id_j = topology.GetMolecularId(pts_id_j);
          if (molecular_id_i == molecular_id_j)
            return;

          auto pts_type_j = topology.GetAtomsType(pts_id_j);
          auto eps_j = topology.GetEpsilon(pts_type_j);
          auto sigma_j = topology.GetSigma(pts_type_j);
          auto r_ij = p_j - p_i;

          LJforce += force_function.ComputeLJForce(r_ij, eps_i, eps_j, sigma_i, sigma_j, _cut_off);
        };
        locator.ExecuteOnNeighbor(atoms_id, function);
      }
      Real _cut_off;
    };

    struct ComputeSpecialBondsLJWorklet : vtkm::worklet::WorkletMapField
    {
      ComputeSpecialBondsLJWorklet(const Real& cut_off)
        : _cut_off(cut_off)
      {
      }

      using ControlSignature = void(FieldIn atoms_id,
                                    ExecObject locator,
                                    ExecObject topology,
                                    ExecObject force_function,
                                    FieldIn special_ids,
                                    FieldIn special_weights,
                                    FieldOut LJforce);
      using ExecutionSignature = void(_1, _2, _3, _4, _5,_6,_7);

      template<typename idsType, typename weightsType>
      VTKM_EXEC void operator()(const Id atoms_id,
                                const ExecPointLocator& locator,
                                const ExecTopology& topology,
                                const ExecForceFunction& force_function,
                                const idsType& special_ids,
                                const weightsType& special_weights,
                                Vec3f& LJforce) const
      {
        LJforce = { 0, 0, 0 };
        const auto& molecular_id_i = topology.GetMolecularId(atoms_id);
        const auto& pts_type_i = topology.GetAtomsType(atoms_id);
        auto eps_i = topology.GetEpsilon(pts_type_i);
        auto sigma_i = topology.GetSigma(pts_type_i);
        auto num_components = special_ids.GetNumberOfComponents();

        auto function = [&](const Vec3f& p_i, const Vec3f& p_j, const Id& pts_id_j)
        {
          auto molecular_id_j = topology.GetMolecularId(pts_id_j);

          //if (molecular_id_i == molecular_id_j)
          //  return;

          auto pts_type_j = topology.GetAtomsType(pts_id_j);
          auto eps_j = topology.GetEpsilon(pts_type_j);
          auto sigma_j = topology.GetSigma(pts_type_j);
          auto r_ij = p_j - p_i;

          Vec3f forceij =
            force_function.ComputeLJForce(r_ij, eps_i, eps_j, sigma_i, sigma_j, _cut_off);

          auto weight = 1.0;
          for (auto i = 0; i < num_components; ++i)
          {
            if (special_ids[i] == pts_id_j)
            {
              weight = special_weights[i];
            }
          }

          LJforce += weight * forceij;
        };

        locator.ExecuteOnNeighbor(atoms_id, function);
      }
      Real _cut_off;
    };

    struct ComputeRBLNeighboursOnceWorklet : vtkm::worklet::WorkletMapField
    {
      ComputeRBLNeighboursOnceWorklet(const Id& rs_num, const Id& pice_num)
        : _rs_num(rs_num),_pice_num(pice_num)
      {
      }

      using ControlSignature = void(FieldIn atoms_id,
                                    ExecObject locator,
                                    FieldInOut id_verletlist_group,
                                    FieldInOut num_verletlist,
                                    FieldInOut offset_verletlist_group);
      using ExecutionSignature = void(_1, _2, _3, _4, _5);

      template<typename NeighbourGroupVecType, typename numType, typename CoordOffsetType>
      VTKM_EXEC void operator()(const Id atoms_id,
                                const ExecPointLocator& locator,
                                NeighbourGroupVecType& id_verletlist,
                                numType& num_verletlist,
                                CoordOffsetType& offset_verletlist) const
      {
        Id index_shell = 0;
        Id index_random = 0;
        Id index_rs = 0;
        auto p_i = locator.GetPtsPosition(atoms_id);
        vtkm::Id3 p_i_cell = locator.InWhichCell(p_i);
        auto num_cycles = locator.GetNumCycles();
        auto rc = locator._cut_off;
        auto rs = locator._rs;
        auto rc_2 = rc * rc;
        auto rs_2 = rs * rs;
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
                const Real dis_2 = r_ij[0] * r_ij[0] + r_ij[1] * r_ij[1] + r_ij[2] * r_ij[2];

                if (rc_2 - dis_2 > 0.0001 && dis_2 - rs_2 > 0.0001)
                {
                  if ( index_shell % _pice_num ==0) // 根据壳内粒子总数的百分比搜寻：index_shell % _pice_num ==0；
                  {
                    id_verletlist[_rs_num + index_random] = pts_id_j;
                    offset_verletlist[_rs_num + index_random] = coord_offset;
                    index_random++;
                  }
                  index_shell++;
                }
                else if (rs_2 - dis_2 > 0.0001 && dis_2 > 0.0001)  // 小球内粒子
                {
                    id_verletlist[index_rs] = pts_id_j;
                    offset_verletlist[index_rs] = coord_offset;
                    index_rs++;
                }       
              }
            }
          }
        }
        num_verletlist[0] = index_rs;
        num_verletlist[1] = index_random;
      }
      Id _rs_num;
      Id _pice_num;
    };

    struct ComputeNearForceRBLERFWorklet : vtkm::worklet::WorkletMapField
    {
      ComputeNearForceRBLERFWorklet(const Id& rs_num,
                                    const Id& pice_num,
                                    const Real& qqr2e,
                                    const Vec3f& box)
        : _rs_num(rs_num)
        , _pice_num(pice_num)
        , _qqr2e(qqr2e)
        , _box(box)
      {
      }

      using ControlSignature = void(FieldIn atoms_id,
                                    ExecObject locator,
                                    ExecObject topology,
                                    ExecObject force_function,
                                    ExecObject static_table,
                                    FieldInOut id_verletlist_group,
                                    FieldInOut num_verletlist,
                                    FieldInOut offset_verletlist_group,
                                    FieldOut corr_force);
      using ExecutionSignature = void(_1, _2, _3, _4, _5, _6, _7,_8,_9);

      template<typename NeighbourGroupVecType, typename numType, typename CoordOffsetType>
      VTKM_EXEC void operator()(const Id atoms_id,
                                const ExecPointLocator& locator,
                                const ExecTopology& topology,
                                const ExecForceFunction& force_function,
                                const ExecStaticTable& static_table,
                                const NeighbourGroupVecType& id_verletlist,
                                const numType& num_verletlist,
                                const CoordOffsetType& offset_verletlist,
                                Vec3f& corr_force) const
      {
        vtkm::Vec3f rc_force_lj = { 0, 0, 0 };
        vtkm::Vec3f rcs_force_lj = { 0, 0, 0 };
        vtkm::Vec3f rc_force_eleforce = { 0, 0, 0 };
        vtkm::Vec3f rcs_force_eleforce = { 0, 0, 0 };

        const auto& molecular_id_i = topology.GetMolecularId(atoms_id);
        const auto& pts_type_i = topology.GetAtomsType(atoms_id);
        auto eps_i = topology.GetEpsilon(pts_type_i);
        auto sigma_i = topology.GetSigma(pts_type_i);
        auto rc = locator._cut_off;
        auto rs = locator._rs;
        auto charge_pi = topology.GetCharge(atoms_id);

        auto function = [&](const Vec3f& p_i, const Vec3f& p_j, const Id& pts_id_j, const Id& flag)
        {
          auto molecular_id_j = topology.GetMolecularId(pts_id_j);
          auto pts_type_j = topology.GetAtomsType(pts_id_j);
          auto eps_j = topology.GetEpsilon(pts_type_j);
          auto sigma_j = topology.GetSigma(pts_type_j);
          //auto r_ij = p_j - p_i;
          auto r_ij = locator.MinDistanceVec(p_j, p_i, _box);

          auto charge_pj = topology.GetCharge(pts_id_j);
          Real dispij = vtkm::Magnitude(r_ij);
          IdComponent index_table_pij = static_table.Extract(dispij);
          auto table_pij = static_table.TableGnearValue(dispij, index_table_pij);

          if (flag == 1)
          {
            rc_force_eleforce += force_function.ComputeNearEnergyForceERF(p_i, p_j, charge_pi, charge_pj, table_pij);

            if (molecular_id_i == molecular_id_j)
              return;
            rc_force_lj += force_function.ComputeLJForce(r_ij, eps_i, eps_j, sigma_i, sigma_j, rc); 
          }
          if (flag == 2)
          {
            force_function.ComputeNearEnergyForceERF_box(r_ij, charge_pi, charge_pj, table_pij);

            if (molecular_id_i == molecular_id_j)
              return;
            rcs_force_lj += _pice_num * force_function.ComputeLJForceRcs(r_ij, eps_i, eps_j, sigma_i, sigma_j, rc, rs);
          }
        };

        auto p_i = locator.GetPtsPosition(atoms_id);
        for (Id p = 0; p < num_verletlist[0]; p++)
        {
          Id id_rs = id_verletlist[p];
          auto p_j = locator.GetPtsPosition(id_rs) - offset_verletlist[p];
          function(p_i, p_j, id_rs, 1);
        }

        for (Id p = 0; p < num_verletlist[1]; p++)
        {
          Id id_random = id_verletlist[_rs_num + p];
          auto p_j = locator.GetPtsPosition(id_random) - offset_verletlist[_rs_num + p];
          function(p_i, p_j, id_random, 2);
        }

        // Individual output requirements can be used
        auto LJforce = rc_force_lj + rcs_force_lj;
        auto nearele_force = _qqr2e * (rc_force_eleforce + rcs_force_eleforce);
        corr_force = LJforce + nearele_force;
      }
      Id _rs_num;
      Id _pice_num;
      Real _qqr2e;
      Vec3f _box;
    };

    struct ComputeNearForceRBLERFSpecialBondsWorklet : vtkm::worklet::WorkletMapField
    {
      ComputeNearForceRBLERFSpecialBondsWorklet(const Id& rs_num,
                                                const Id& pice_num,
                                                const Real& qqr2e,
                                                const Vec3f& box)
        : _rs_num(rs_num)
        , _pice_num(pice_num)
        , _qqr2e(qqr2e)
        , _box(box)
      {
      }

      using ControlSignature = void(FieldIn atoms_id,
                                    ExecObject locator,
                                    ExecObject topology,
                                    ExecObject force_function,
                                    ExecObject static_table,
                                    FieldInOut id_verletlist_group,
                                    FieldInOut num_verletlist,
                                    FieldInOut offset_verletlist_group,
                                    FieldIn special_ids,
                                    FieldIn special_weights,
                                    FieldOut corr_force);
      using ExecutionSignature = void(_1, _2, _3, _4, _5, _6, _7, _8, _9, _10, _11);

      template<typename NeighbourGroupVecType,
               typename numType,
               typename CoordOffsetType,
               typename idsType,
               typename weightsType>
      VTKM_EXEC void operator()(const Id atoms_id,
                                const ExecPointLocator& locator,
                                const ExecTopology& topology,
                                const ExecForceFunction& force_function,
                                const ExecStaticTable& static_table,
                                const NeighbourGroupVecType& id_verletlist,
                                const numType& num_verletlist,
                                const CoordOffsetType& offset_verletlist,
                                const idsType& special_ids,
                                const weightsType& special_weights,
                                Vec3f& corr_force) const
      {
        vtkm::Vec3f rc_force_lj = { 0, 0, 0 };
        vtkm::Vec3f rcs_force_lj = { 0, 0, 0 };
        vtkm::Vec3f rc_force_eleforce = { 0, 0, 0 };
        vtkm::Vec3f rcs_force_eleforce = { 0, 0, 0 };

        const auto& molecular_id_i = topology.GetMolecularId(atoms_id);
        const auto& pts_type_i = topology.GetAtomsType(atoms_id);
        auto eps_i = topology.GetEpsilon(pts_type_i);
        auto sigma_i = topology.GetSigma(pts_type_i);
        auto rc = locator._cut_off;
        auto rs = locator._rs;
        auto charge_pi = topology.GetCharge(atoms_id);

        auto num_components = special_ids.GetNumberOfComponents();

        auto function = [&](const Vec3f& p_i, const Vec3f& p_j, const Id& pts_id_j, const Id& flag)
        {
          auto molecular_id_j = topology.GetMolecularId(pts_id_j);
          auto pts_type_j = topology.GetAtomsType(pts_id_j);
          auto eps_j = topology.GetEpsilon(pts_type_j);
          auto sigma_j = topology.GetSigma(pts_type_j);
          //auto r_ij = p_j - p_i;
          auto r_ij = locator.MinDistanceVec(p_j, p_i, _box);

          auto charge_pj = topology.GetCharge(pts_id_j);
          Real dispij = vtkm::Magnitude(r_ij);
          IdComponent index_table_pij = static_table.Extract(dispij);
          auto table_pij = static_table.TableGnearValue(dispij, index_table_pij);

          if (flag == 1)
          {
            rc_force_eleforce +=
              force_function.ComputeNearEnergyForceERF(p_i, p_j, charge_pi, charge_pj, table_pij);

            /*if (molecular_id_i == molecular_id_j)
              return;
            rc_force_lj += force_function.ComputeLJForce(r_ij, eps_i, eps_j, sigma_i, sigma_j, rc);*/
            Vec3f forceij = force_function.ComputeLJForce(r_ij, eps_i, eps_j, sigma_i, sigma_j, rc);
            auto weight = 1.0;
            for (auto i = 0; i < num_components; ++i)
            {
              if (special_ids[i] == pts_id_j)
              {
                weight = special_weights[i];
              }
            }

            rc_force_lj += weight * forceij;
          }
          if (flag == 2)
          {
            rcs_force_eleforce += _pice_num * force_function.ComputeNearEnergyForceRcsERF(p_i, p_j, charge_pi, charge_pj, table_pij, rs);

            /*if (molecular_id_i == molecular_id_j)
              return;
            rcs_force_lj += _pice_num *
              force_function.ComputeLJForceRcs(r_ij, eps_i, eps_j, sigma_i, sigma_j, rc, rs);*/
            Vec3f forceij = _pice_num *
              force_function.ComputeLJForceRcs(r_ij, eps_i, eps_j, sigma_i, sigma_j, rc, rs);
            auto weight = 1.0;
            for (auto i = 0; i < num_components; ++i)
            {
              if (special_ids[i] == pts_id_j)
              {
                weight = special_weights[i];
              }
            }

            rcs_force_lj += weight * forceij;
          }
        };

        auto p_i = locator.GetPtsPosition(atoms_id);
        for (Id p = 0; p < num_verletlist[0]; p++)
        {
          Id id_rs = id_verletlist[p];
          auto p_j = locator.GetPtsPosition(id_rs) - offset_verletlist[p];
          function(p_i, p_j, id_rs, 1);
        }

        for (Id p = 0; p < num_verletlist[1]; p++)
        {
          Id id_random = id_verletlist[_rs_num + p];
          auto p_j = locator.GetPtsPosition(id_random) - offset_verletlist[_rs_num + p];
          function(p_i, p_j, id_random, 2);
        }

        // Individual output requirements can be used
        auto LJforce = rc_force_lj + rcs_force_lj;
        auto nearele_force = _qqr2e * (rc_force_eleforce + rcs_force_eleforce);
        corr_force = LJforce + nearele_force;
      }
      Id _rs_num;
      Id _pice_num;
      Real _qqr2e;
      Vec3f _box;
    };

    struct ComputeNearForceRBLWorklet : vtkm::worklet::WorkletMapField
    {
      ComputeNearForceRBLWorklet(const Id& rs_num, const Id& pice_num, const Real& qqr2e)
        : _rs_num(rs_num)
        , _pice_num(pice_num)
        , _qqr2e(qqr2e)
      {
      }

      using ControlSignature = void(FieldIn atoms_id,
                                    ExecObject locator,
                                    ExecObject topology,
                                    ExecObject force_function,
                                    FieldInOut id_verletlist_group,
                                    FieldInOut num_verletlist,
                                    FieldInOut offset_verletlist_group,
                                    FieldOut corr_force);
      using ExecutionSignature = void(_1, _2, _3, _4, _5, _6, _7, _8);

      template<typename NeighbourGroupVecType, typename numType, typename CoordOffsetType>
      VTKM_EXEC void operator()(const Id atoms_id,
                                const ExecPointLocator& locator,
                                const ExecTopology& topology,
                                const ExecForceFunction& force_function,
                                const NeighbourGroupVecType& id_verletlist,
                                const numType& num_verletlist,
                                const CoordOffsetType& offset_verletlist,
                                Vec3f& corr_force) const
      {
        vtkm::Vec3f rc_force_lj = { 0, 0, 0 };
        vtkm::Vec3f rcs_force_lj = { 0, 0, 0 };
        vtkm::Vec3f rc_force_eleforce = { 0, 0, 0 };
        vtkm::Vec3f rcs_force_eleforce = { 0, 0, 0 };

        const auto& molecular_id_i = topology.GetMolecularId(atoms_id);
        const auto& pts_type_i = topology.GetAtomsType(atoms_id);
        auto eps_i = topology.GetEpsilon(pts_type_i);
        auto sigma_i = topology.GetSigma(pts_type_i);
        auto rc = locator._cut_off;
        auto rs = locator._rs;
        auto charge_pi = topology.GetCharge(atoms_id);

        auto function = [&](const Vec3f& p_i, const Vec3f& p_j, const Id& pts_id_j, const Id& flag)
        {
          auto molecular_id_j = topology.GetMolecularId(pts_id_j);
          auto pts_type_j = topology.GetAtomsType(pts_id_j);
          auto eps_j = topology.GetEpsilon(pts_type_j);
          auto sigma_j = topology.GetSigma(pts_type_j);
          auto r_ij = p_j - p_i;
          auto charge_pj = topology.GetCharge(pts_id_j);

          if (flag == 1)
          {
            rc_force_eleforce += force_function.ComputeNearEnergyForce(p_i,p_j,charge_pi,charge_pj);

            if (molecular_id_i == molecular_id_j)
              return;
            rc_force_lj += force_function.ComputeLJForce(r_ij, eps_i, eps_j, sigma_i, sigma_j, rc);
          }
          if (flag == 2)
          {
            rcs_force_eleforce += _pice_num * force_function.ComputeNearEnergyForceRcs(p_i, p_j, charge_pi, charge_pj, rs);

            if (molecular_id_i == molecular_id_j)
              return;
            rcs_force_lj += _pice_num *
              force_function.ComputeLJForceRcs(r_ij, eps_i, eps_j, sigma_i, sigma_j, rc, rs);
          }
        };

        auto p_i = locator.GetPtsPosition(atoms_id);
        for (Id p = 0; p < num_verletlist[0]; p++)
        {
          Id id_rs = id_verletlist[p];
          auto p_j = locator.GetPtsPosition(id_rs) - offset_verletlist[p];
          function(p_i, p_j, id_rs, 1);
        }

        for (Id p = 0; p < num_verletlist[1]; p++)
        {
          Id id_random = id_verletlist[_rs_num + p];
          auto p_j = locator.GetPtsPosition(id_random) - offset_verletlist[_rs_num + p];
          function(p_i, p_j, id_random, 2);
        }

        // Individual output requirements can be used
        auto LJforce = rc_force_lj + rcs_force_lj;
        auto nearele_force = _qqr2e * (rc_force_eleforce + rcs_force_eleforce);
        corr_force = LJforce + nearele_force;
      }
      Id _rs_num;
      Id _pice_num;
      Real _qqr2e;
    };

   struct ComputeEAMfp_Worklet : vtkm::worklet::WorkletMapField
    {
      ComputeEAMfp_Worklet(const Real& rc,
                           const Vec3f& box,
                           const Real& rs,
                           const Id& rs_num,
                           const Id& pice_num)
        : _rc(rc)
        , _box(box)
        , _rs(rs)
        , _rs_num(rs_num)
        , _pice_num(pice_num)
      {
      }
      
      using ControlSignature = void(FieldIn atoms_id,
                                    WholeArrayIn rhor_spline,
                                    WholeArrayIn frho_spline,
                                    ExecObject locator,
                                    ExecObject force_function,
                                    FieldInOut id_verletlist_group,
                                    FieldInOut num_verletlist,
                                    FieldInOut offset_verletlist_group,
                                    FieldOut EAM_fp);
      using ExecutionSignature = void(_1, _2, _3, _4, _5, _6, _7, _8, _9);

      template<typename NeighbourGroupVecType,
               typename NumType,
               typename CoordOffsetType,
               typename SpileType,
               typename frhoType>
      VTKM_EXEC void operator()(const Id atoms_id,
                                const SpileType& rhor_spline,
                                const frhoType& frho_spline,
                                const ExecPointLocator& locator,
                                const ExecForceFunction& force_function,
                                const NeighbourGroupVecType& id_verletlist,
                                const NumType& num_verletlist,
                                const CoordOffsetType& offset_verletlist,
                                Real& EAM_fp) const
      {
        Real rho = 0;
        auto function = [&](const Vec3f& p_i, const Vec3f& p_j, const Id& pts_id_j, const Id& flag)
        {
          //auto r_ij = p_j - p_i;
          auto r_ij = locator.MinDistanceVec(p_i, p_j, _box);
          if (flag == 1)
          {
            rho += force_function.ComputeEAMrho(_rc, r_ij, rhor_spline);
          }
          if (flag == 2)
          {
            rho += _pice_num * force_function.ComputeEAMrhoRs(_rc, _rs, r_ij, rhor_spline);
          }
        };
        auto p_i = locator.GetPtsPosition(atoms_id);
        for (Id p = 0; p < num_verletlist[0]; p++)
        {
          Id id_rs = id_verletlist[p];
          auto p_j = locator.GetPtsPosition(id_rs) - offset_verletlist[p];
          function(p_i, p_j, id_rs, 1);
        }

        for (Id p = 0; p < num_verletlist[1]; p++)
        {
          Id id_random = id_verletlist[_rs_num + p];
          auto p_j = locator.GetPtsPosition(id_random) - offset_verletlist[_rs_num + p];
          function(p_i, p_j, id_random, 2);
        }
        //EAM_rho = rho;
        EAM_fp = force_function.ComputeEAMfp(atoms_id, rho, frho_spline);
      }
      Real _rc;
      Vec3f _box;
      Real _rs;
      Id _rs_num;
      Id _pice_num;
    };

    struct ComputeEAMRBLWorklet : vtkm::worklet::WorkletMapField
    {
      ComputeEAMRBLWorklet(const Real& rc,
                           const Vec3f& box,
                           const Real& rs,
                           const Id& rs_num,
                           const Id& pice_num)
        : _rc(rc)
        , _box(box)
        , _rs(rs)
        , _rs_num(rs_num)
        , _pice_num(pice_num)
      {
      }

      using ControlSignature = void(FieldIn atoms_id,
                                    WholeArrayIn rhor_spline,
                                    WholeArrayIn z2r_spline,
                                    WholeArrayIn EAM_fp,
                                    ExecObject locator,
                                    ExecObject force_function,
                                    FieldInOut id_verletlist_group,
                                    FieldInOut num_verletlist,
                                    FieldInOut offset_verletlist_group,
                                    FieldOut corr_force);
      using ExecutionSignature = void(_1, _2, _3, _4, _5, _6, _7, _8, _9,_10);

      template<typename NeighbourGroupVecType,
               typename NumType,
               typename CoordOffsetType,
               typename SpileType,
               typename z2rSpileType,
               typename EAM_fpType>
      VTKM_EXEC void operator()(const Id atoms_id,
                                const SpileType& rhor_spline,
                                const z2rSpileType& z2r_spline,
                                const EAM_fpType& EAM_fp,
                                const ExecPointLocator& locator,
                                const ExecForceFunction& force_function,
                                const NeighbourGroupVecType& id_verletlist,
                                const NumType& num_verletlist,
                                const CoordOffsetType& offset_verletlist,
                                Vec3f& corr_force) const
      { 
         Vec3f rc_force = { 0, 0, 0 };
         Vec3f rcs_force = { 0, 0, 0 };

         auto function_force = [&](const Vec3f& p_i, const Vec3f& p_j, const Id& pts_id_j, const Id& flag)
         {
           auto r_ij = locator.MinDistanceVec(p_i, p_j, _box);
           if (flag == 1)
           {
             rc_force += force_function.ComputeEAMforce(
               _rc, atoms_id, pts_id_j, r_ij, EAM_fp, rhor_spline, z2r_spline);
           }
           if (flag == 2)
           {
             rcs_force += _pice_num *
               force_function.ComputeEAMforceRBL(
                 _rc, _rs, atoms_id, pts_id_j, r_ij, EAM_fp, rhor_spline, z2r_spline);
           }
         };

         auto p_i = locator.GetPtsPosition(atoms_id);
         for (Id p = 0; p < num_verletlist[0]; p++)
         {
            Id id_rs = id_verletlist[p];
            auto p_j = locator.GetPtsPosition(id_rs) - offset_verletlist[p];
            function_force(p_i, p_j, id_rs, 1);
         }

         for (Id p = 0; p < num_verletlist[1]; p++)
         {
            Id id_random = id_verletlist[_rs_num + p];
            auto p_j = locator.GetPtsPosition(id_random) - offset_verletlist[_rs_num + p];
            function_force(p_i, p_j, id_random, 2);
         }

         corr_force = rc_force + rcs_force;
      }
      Real _rc;
      Vec3f _box;
      Real _rs;
      Id _rs_num;
      Id _pice_num;
    };

    struct ComputeLJForceRBLWorklet : vtkm::worklet::WorkletMapField
    {
      ComputeLJForceRBLWorklet(const Id& rs_num, const Id& pice_num, const Vec3f& box)
        : _rs_num(rs_num)
        , _pice_num(pice_num)
        , _box(box)
      {
      }

      using ControlSignature = void(FieldIn atoms_id,
                                    ExecObject locator,
                                    ExecObject topology,
                                    ExecObject force_function,
                                    FieldInOut id_verletlist_group,
                                    FieldInOut num_verletlist,
                                    FieldInOut offset_verletlist_group,
                                    FieldOut corr_ljforce);
      using ExecutionSignature = void(_1, _2, _3, _4, _5, _6, _7, _8);

      template<typename NeighbourGroupVecType, typename numType, typename CoordOffsetType>
      VTKM_EXEC void operator()(const Id atoms_id,
                                const ExecPointLocator& locator,
                                const ExecTopology& topology,
                                const ExecForceFunction& force_function,
                                const NeighbourGroupVecType& id_verletlist,
                                const numType& num_verletlist,
                                const CoordOffsetType& offset_verletlist,
                                Vec3f& corr_ljforce) const
      {
        vtkm::Vec3f rc_force_lj = { 0, 0, 0 };
        vtkm::Vec3f rcs_force_lj = { 0, 0, 0 };

        const auto& molecular_id_i = topology.GetMolecularId(atoms_id);
        const auto& pts_type_i = topology.GetAtomsType(atoms_id);
        auto eps_i = topology.GetEpsilon(pts_type_i);
        auto sigma_i = topology.GetSigma(pts_type_i);
        auto rc = locator._cut_off;
        auto rs = locator._rs;
        auto charge_pi = topology.GetCharge(atoms_id);

        auto function = [&](const Vec3f& p_i, const Vec3f& p_j, const Id& pts_id_j, const Id& flag)
        {
          auto molecular_id_j = topology.GetMolecularId(pts_id_j);
          if (molecular_id_i == molecular_id_j)
            return;
          auto pts_type_j = topology.GetAtomsType(pts_id_j);
          auto eps_j = topology.GetEpsilon(pts_type_j);
          auto sigma_j = topology.GetSigma(pts_type_j);
          //auto r_ij = p_j - p_i;
          auto r_ij = locator.MinDistanceVec(p_j, p_i, _box);

          if (flag == 1)
          {
            rc_force_lj += force_function.ComputeLJForce(r_ij, eps_i, eps_j, sigma_i, sigma_j, rc);
          }
          if (flag == 2)
          {
            rcs_force_lj += _pice_num *
              force_function.ComputeLJForceRcs(r_ij, eps_i, eps_j, sigma_i, sigma_j, rc, rs);
          }
        };

        auto p_i = locator.GetPtsPosition(atoms_id);
        for (Id p = 0; p < num_verletlist[0]; p++)
        {
          Id id_rs = id_verletlist[p];
          auto p_j = locator.GetPtsPosition(id_rs) - offset_verletlist[p];
          function(p_i, p_j, id_rs, 1);
        }

        for (Id p = 0; p < num_verletlist[1]; p++)
        {
          Id id_random = id_verletlist[_rs_num + p];
          auto p_j = locator.GetPtsPosition(id_random) - offset_verletlist[_rs_num + p];
          function(p_i, p_j, id_random, 2);
        }

        // Individual output requirements can be used
        auto LJforce = rc_force_lj + rcs_force_lj;
        corr_ljforce = LJforce;
      }
      Id _rs_num;
      Id _pice_num;
      Vec3f _box;
    };

    struct SumRBLCorrForceWorklet : vtkm::worklet::WorkletMapField
    {
      SumRBLCorrForceWorklet(const vtkm::Vec3f corr_value)
        : _corr_value(corr_value)
      {
      }

      using ControlSignature = void(FieldIn corr_force, FieldOut ljforce);
      using ExecutionSignature = void(_1, _2);
    
    
      VTKM_EXEC void operator()(const Vec3f& corr_force,
                                Vec3f& ljforce) const
      {
        ljforce = corr_force - _corr_value;
      }
      vtkm::Vec3f _corr_value;
    };

    struct ComputeNearElectrostaticsWorklet : vtkm::worklet::WorkletMapField
       {
          ComputeNearElectrostaticsWorklet() {}
          using ControlSignature = void(FieldIn atoms_id,
                                        ExecObject locator,
                                        ExecObject topology,
                                        ExecObject force_function,
                                        FieldOut force); 
          using ExecutionSignature = void(_1, _2, _3, _4, _5);
          
          VTKM_EXEC void operator()(const Id& atoms_id,
                                    const ExecPointLocator& locator,
                                    const ExecTopology& topology,
                                    const ExecForceFunction& force_function,
                                    Vec3f& force) const
          {
            force = { 0, 0, 0 };
            const auto& charge_pi = topology.GetCharge(atoms_id);

            auto function = [&](const Vec3f& p_i, const Vec3f& p_j, const Id& pts_id_j)
            {
              auto charge_pj = topology.GetCharge(pts_id_j);
              force += force_function.ComputeNearEnergyForce(p_i, p_j, charge_pi, charge_pj);
            };
            locator.ExecuteOnNeighbor(atoms_id, function);
          }
       };

    struct ComputeNearElectrostaticsERFWorklet : vtkm::worklet::WorkletMapField
       {
          ComputeNearElectrostaticsERFWorklet() {}
          using ControlSignature = void(FieldIn atoms_id,
                                        ExecObject static_table,
                                        ExecObject locator,
                                        ExecObject topology,
                                        ExecObject force_function,
                                        FieldOut force);
          using ExecutionSignature = void(_1, _2, _3, _4, _5,_6);

          VTKM_EXEC void operator()(const Id& atoms_id,
                                    const ExecStaticTable& static_table,
                                    const ExecPointLocator& locator,
                                    const ExecTopology& topology,
                                    const ExecForceFunction& force_function,
                                    Vec3f& force) const
          {
            force = { 0, 0, 0 };
            const auto& charge_pi = topology.GetCharge(atoms_id);

            auto function = [&](const Vec3f& p_i, const Vec3f& p_j, const Id& pts_id_j)
            {
              auto r_ij = p_j - p_i;
              Real dispij = vtkm::Magnitude(r_ij);
              IdComponent index_table_pij = static_table.Extract(dispij);
              auto table_pij = static_table.TableGnearValue(dispij, index_table_pij);
              auto charge_pj = topology.GetCharge(pts_id_j);
              //force += force_function.ComputeNearEnergyForce(p_i, p_j, charge_pi, charge_pj);
              force += force_function.ComputeNearEnergyForceERF(p_i, p_j, charge_pi, charge_pj,table_pij);        
            };
            locator.ExecuteOnNeighbor(atoms_id, function);
          }
       };

    struct ComputeSpecialCoulWorklet : vtkm ::worklet::WorkletMapField
       {
          ComputeSpecialCoulWorklet(const Vec3f& box)
            : _box(box)
         {
         }

         using ControlSignature = void(FieldIn atoms_id,
                                       FieldIn group_vec,
                                       ExecObject force_function,
                                       ExecObject topology,
                                       ExecObject locator,
                                       FieldOut force,
                                       FieldOut SpecialCoulVirial);
         using ExecutionSignature = void(_1, _2, _3, _4, _5, _6,_7);
         template<typename GroupVecType>

         VTKM_EXEC void operator()(const Id& atoms_id,
                                   const GroupVecType& group_vec,
                                   ExecForceFunction& force_function,
                                   const ExecTopology& topology,
                                   const ExecPointLocator& locator,
                                   Vec3f& force,
                                   Vec6f&  SpecialCoulVirial) const
         {
           auto num = group_vec.GetNumberOfComponents();
           force = { 0, 0, 0 };
           SpecialCoulVirial = { 0, 0, 0, 0, 0, 0 };

           if (num < 3)
           {
             force = 0;
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
               Vec3f rij = locator.MinDistanceVec(atoms_coord_i, atoms_coord_j, _box);
               Real dis_ij = vtkm::Magnitude(rij);
               Real dis_ij3 = vtkm::Pow(dis_ij, 3);
               Real force_component = -332.06371 * charge_p_i * atoms_charge_j / dis_ij3;
               force +=  force_component * rij;
               //
               auto f = -0.5 * force_component * rij;

               SpecialCoulVirial[0] += rij[0] * f[0];
               SpecialCoulVirial[1] += rij[1] * f[1];
               SpecialCoulVirial[2] += rij[2] * f[2];
               SpecialCoulVirial[3] += rij[0] * f[1];
               SpecialCoulVirial[4] += rij[0] * f[2];
               SpecialCoulVirial[5] += rij[1] * f[2];
             }
           }
         }
         Vec3f _box;
       };

    struct ComputeSpecialCoulGeneralWorklet : vtkm ::worklet::WorkletMapField
    {
         ComputeSpecialCoulGeneralWorklet(const Vec3f& box)
           : _box(box)
         {
         }
       
         using ControlSignature = void(FieldIn atoms_id,
                                       FieldIn group_vec,
                                       ExecObject force_function,
                                       ExecObject topology,
                                       ExecObject locator,
                                       FieldIn special_ids,
                                       FieldIn special_weights,
                                       FieldOut force,
                                       FieldOut SpecialCoulVirial);
         using ExecutionSignature = void(_1, _2, _3, _4, _5, _6 ,_7,_8,_9);

         template<typename GroupVecType, typename idsType, typename weightsType>
         VTKM_EXEC void operator()(const Id& atoms_id,
                                   const GroupVecType& group_vec,
                                    ExecForceFunction& force_function,
                                   const ExecTopology& topology,
                                   const ExecPointLocator& locator,
                                   const idsType& special_ids,
                                   const weightsType& special_weights,
                                   Vec3f& force,
                                   Vec6f&  SpecialCoulVirial) const
         { 
            auto num = group_vec.GetNumberOfComponents();           
            force = { 0, 0, 0 };
            SpecialCoulVirial = { 0, 0, 0, 0, 0, 0 };

            if (num < 3)
            {
              force = 0;
            }
            else
            {
              auto charge_p_i = topology.GetCharge(atoms_id);
              auto atoms_coord_i = locator.GetPtsPosition(atoms_id);
              auto num_components = special_ids.GetNumberOfComponents();

              for (Id j = 0; j < num; ++j)
              {
                auto id = group_vec[j];
                if (atoms_id == id)
                  continue;
                
                
                auto atoms_charge_j = topology.GetCharge(id);
                auto atoms_coord_j = locator.GetPtsPosition(id);
                Vec3f rij = locator.MinDistanceVec(atoms_coord_i, atoms_coord_j, _box);
                Real dis_ij = vtkm::Magnitude(rij);
                Real dis_ij3 = vtkm::Pow(dis_ij, 3);
                Real force_component = -332.06371 * charge_p_i * atoms_charge_j / dis_ij3;

                Real weight = 1.0;
                for (auto i = 0; i < num_components; ++i)
                {
                  if (special_ids[i] == id)
                  {
                    weight  = special_weights[i];
                  }
                }

                force += (1.0 - weight) * force_component * rij;

                //
                auto f = -0.5 * (1.0 - weight) * force_component * rij;

                SpecialCoulVirial[0] += rij[0] * f[0];
                SpecialCoulVirial[1] += rij[1] * f[1];
                SpecialCoulVirial[2] += rij[2] * f[2];
                SpecialCoulVirial[3] += rij[0] * f[1];
                SpecialCoulVirial[4] += rij[0] * f[2];
                SpecialCoulVirial[5] += rij[1] * f[2];
              }
            } 
         }
         Vec3f _box;
    };

    struct ComputeNewRBEForceWorklet : vtkm ::worklet::WorkletMapField
       {
         ComputeNewRBEForceWorklet(const IdComponent& RBEP)
           : _P(RBEP)
         {
         }
       
         using ControlSignature = void(FieldIn atoms_id,
                                       WholeArrayIn psample,
                                       WholeArrayIn whole_rhok,
                                       ExecObject force_function,
                                       ExecObject topology,
                                       ExecObject locator,
                                       FieldOut force);
         using ExecutionSignature = void(_1, _2, _3, _4, _5, _6, _7);
       
         template<typename WholePType, typename WholeRhokType>
         VTKM_EXEC void operator()(const Id& atoms_id,
                                   const WholePType& psample,
                                   const WholeRhokType& whole_rhok,
                                   ExecForceFunction& force_function,
                                   const ExecTopology& topology,
                                   ExecPointLocator& locator,
                                   Vec3f& force) const
         {
           force = { 0, 0, 0 };
           const auto& current_pts = locator.GetPtsPosition(atoms_id);
           const auto& current_charge = topology.GetCharge(atoms_id);
       
           for (int i = 0; i < _P; i++)
           {
             Vec2f rhok_ri = whole_rhok.Get(i);
             Vec3f kl = psample.Get(i);
       
             force += force_function.ComputeRBEForceSum(kl, current_pts, current_charge, rhok_ri);
           }
           force = force_function.ComputeRBEForce(_P, force);
         }
       
         Id _P;
       };

    struct InitPositionWorklet : vtkm::worklet::WorkletVisitCellsWithPoints
    {
      using ControlSignature = void(CellSetIn, FieldInPoint, FieldOut);
      using ExecutionSignature = void(_2, _3);
    
      template<typename PType>
      VTKM_EXEC void operator()(const PType& pts, vtkm::Vec3f& centroid) const
      {
        if (pts.GetNumberOfComponents() == 4)
        {
          centroid = (pts[0] + pts[1] + pts[2] + pts[3]) / 4;
        }
        else if (pts.GetNumberOfComponents() == 8)
        {
          centroid = (pts[0] + pts[1] + pts[2] + pts[3] + pts[4] + pts[5] + pts[6] + pts[7]) / 8;
        }
        else
        {
          RaiseError("unfinished");
        }
      }
    };
    
    struct InitialCondtionWorklet : vtkm::worklet::WorkletMapField
    {
      using ControlSignature = void(FieldIn pts, FieldIn v, FieldOut mass, FieldOut veclocity);
      using ExecutionSignature = void(_1, _2, _3, _4);
    
      template<typename CoordType>
      VTKM_EXEC void operator()(const CoordType& pts, const Vec3f& v, Real& mass, Vec3f& velocity) const
      {
        mass = 1.;
        velocity = v;
      }
    };
    
    struct ComputeChargeStructureFactorComponentWorklet : vtkm ::worklet::WorkletMapField
    {
      ComputeChargeStructureFactorComponentWorklet(const Vec3f& k)
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

    struct ComputePnumberChargeStructureFactorWorklet : vtkm ::worklet::WorkletMapField
    {
      ComputePnumberChargeStructureFactorWorklet(const Vec3f& box,
                                                 const Id& pnumber,
                                                 const Id& pos_number)
        : _box(box)
        , _pnumber(pnumber)
        , _pos_num(pos_number)
      {
      }

      using ControlSignature = void(FieldIn atoms_id,
                                    FieldIn current_pts,
                                    FieldIn current_charge,
                                    WholeArrayIn psample,
                                   // WholeArrayInOut psamplekey_group,
                                    WholeArrayInOut rhok_Re_group,
                                    WholeArrayInOut rhok_Im_group);

      using ExecutionSignature = void(_1, _2, _3, _4, _5, _6);

      template<typename CoordType,  typename Densitytype, typename SampleType>
      VTKM_EXEC void operator()(const Id atoms_id,
                                const CoordType& current_pts,
                                const Real& current_charge,
                                const SampleType& psample,
                               // Keytype& psamplekey_group,
                                Densitytype& rhok_Re_group,
                                Densitytype& rhok_Im_group) const
      {
        //auto p_i = locator.GetPtsPosition(atoms_id);
        //Vec3f r_i = current_pts;
        /*Real rho_Re_i = 0.0;
        Real rho_Im_i = 0.0;
        density_Re = current_charge * vtkm::Cos(vtkm::Dot(_K, r_i));
        density_Im = current_charge * vtkm::Sin(vtkm::Dot(_K, r_i));*/

        for (Id i = 0; i < _pnumber; i++)
        {
          auto kl = psample.Get(i);
          for (Id j = 0; j < 3; j++)
          {
                kl[j] = 2 * vtkm::Pi() * kl[j] / _box[j];
          }

          auto index = atoms_id + i * _pos_num;
          //psamplekey_group.Set(index,i);
          rhok_Re_group.Set(index,current_charge * vtkm::Cos(vtkm::Dot(kl, current_pts)));
          rhok_Im_group.Set(index,current_charge * vtkm::Sin(vtkm::Dot(kl, current_pts)));
	            //rhok_Re_group.Set(index,current_charge * 1.0);
          //rhok_Im_group.Set(index,current_charge * 1.0);
          //psamplekey_group[index] = i;
          //rhok_Re_group[index] = current_charge * vtkm::Cos(vtkm::Dot(kl, r_i));
          //rhok_Im_group[index] = current_charge * vtkm::Sin(vtkm::Dot(kl, r_i));
        }
      }


      Vec3f _box;
      Id _pnumber;
      Id _pos_num;
    };

    struct ChangePnumberChargeStructureFactorRBEWorklet : vtkm ::worklet::WorkletMapField
    {

      using ControlSignature = void(FieldIn rhok_Re,
                                    FieldIn rhok_Im,
                                    FieldOut new_rhok);

      using ExecutionSignature = void(_1, _2, _3);

      template<typename Densitytype, typename Rhoktype>
      VTKM_EXEC void operator()(const Densitytype& rhok_Re,
                                const Densitytype& rhok_Im,
                                Rhoktype& new_rhok) const
      {
        new_rhok = Vec2f{ rhok_Re, rhok_Im };
        //new_rhok[0] = rhok_Re;
        //new_rhok[1] = rhok_Im;
      }

    };

    struct UnderdampedLangevinWorklet : vtkm::worklet::WorkletMapField
    {
      UnderdampedLangevinWorklet(const Vec3f& gaussian,
                                 const Real& kBT,
                                 const Real& gamma,
                                 const Real& dt)
    : _gaussian(gaussian)
    , _kBT(kBT)
    , _gamma(gamma)
    , _dt(dt)
      {
      }
    
      using ControlSignature = void(const FieldIn mass, const FieldIn velocity, FieldInOut outforce);
      using ExecutionSignature = void(_1, _2,_3);
    
      template<typename Mass, typename VelocityType>
      VTKM_EXEC void operator()(const Mass& mass, const VelocityType& velocity, Vec3f& outforce) const
      {
        outforce = outforce - mass * _gamma * velocity +
          vtkm::Sqrt(2.0 * mass * _gamma * _kBT / _dt) * _gaussian;
      }
    
      // TODO: optional parameters for configuration files
      Real _kBT;
      Real _dt;
      Real _gamma;    
      Vec3f _gaussian;
    };
    
    struct SumFarNearForceWorklet : vtkm::worklet::WorkletMapField
    {
      using ControlSignature = void(FieldIn farforce,
                                    FieldIn nearforce,
                                    FieldOut allforce);
      using ExecutionSignature = void(_1, _2, _3);
    
    
      VTKM_EXEC void operator()(const Vec3f& farforce,
                                const Vec3f& nearforce,
                                Vec3f& allforce) const
      {
        allforce = farforce + nearforce;
      }
    };

    struct SumVirialWorklet : vtkm::worklet::WorkletMapField
    {
      using ControlSignature = void(FieldIn LJVirial, FieldIn EwaldVirial, FieldOut allVirial);
      using ExecutionSignature = void(_1, _2, _3);


      VTKM_EXEC void operator()(const Vec6f& LJVirial,
                                const Vec6f& EwaldVirial,
                                Vec6f& allVirial) const
      {
        allVirial = LJVirial + EwaldVirial;
      }
    };

    struct SumFarNearLJForceWorklet : vtkm::worklet::WorkletMapField
    {
      using ControlSignature = void(FieldIn farforce,
                                    FieldIn nearforce,
                                    FieldIn ljforce,
                                    FieldOut allforce);
      using ExecutionSignature = void(_1, _2, _3 , _4);


      VTKM_EXEC void operator()(const Vec3f& farforce,
                                const Vec3f& nearforce,
                                const Vec3f& ljforce,
                                Vec3f& allforce) const
      {
        allforce = farforce + nearforce + ljforce;
      }
    };
    
    struct UpdateVelocityWorklet : vtkm::worklet::WorkletMapField
    {
      using ControlSignature = void(FieldIn force, FieldIn mass,FieldInOut velocity);
      using ExecutionSignature = void(_1, _2, _3);
    
      UpdateVelocityWorklet(const Real& dt, const Real& unit_factor)
        : _dt(dt)
        , _unit_factor(unit_factor)
      {
      }
    
      VTKM_EXEC void operator()(const Vec3f& force, const Real& mass, Vec3f& velocity) const
      {
        velocity += 0.5 * force / mass * _dt * _unit_factor;
      }   
 
      Real _dt;
      Real _unit_factor;
    };
    
    struct UpdateVelocityNoseHooverWorklet : vtkm::worklet::WorkletMapField
    {
      using ControlSignature = void(FieldIn force, FieldIn mass,FieldInOut velocity);
      using ExecutionSignature = void(_1, _2, _3);
    
      UpdateVelocityNoseHooverWorklet(const Real& dt, const Real unit_factor, const Real& xi)
    : _dt(dt)
    , _unit_factor(unit_factor)
    , _xi(xi)
      {
      }
        
      VTKM_EXEC void operator()(const Vec3f& force, const Real& mass, Vec3f& velocity) const
      {
        velocity += 0.5 * _dt * (force / mass - _xi * velocity) * _unit_factor;
      }
    
      Real _dt;
      Real _unit_factor;
      Real _xi;
    };
    
    struct ComputerKineticEnergyWorklet : vtkm::worklet::WorkletMapField
    {
      using ControlSignature = void(FieldIn velocity, FieldIn mass, FieldOut sq_velocity);
      using ExecutionSignature = void(_1, _2, _3);
    
      VTKM_EXEC void operator()(const Vec3f& velocity, const Real& mass, Real& sq_velocity) const
      {
        sq_velocity =
          mass * (velocity[0] * velocity[0] + velocity[1] * velocity[1] + velocity[2] * velocity[2]);
      }
    };
    
    struct UpdateVelocityRescaleWorklet : vtkm::worklet::WorkletMapField
    {
      using ControlSignature = void(FieldInOut velocity);
      using ExecutionSignature = void(_1);
    
      UpdateVelocityRescaleWorklet(const Real& coeff)
      : _coeff(coeff)
      {
      }
    
      VTKM_EXEC void operator()(Vec3f& velocity) const { velocity = _coeff * velocity; }
    
      Real _coeff;
    };
    
    struct UpdatePositionWorklet : vtkm::worklet::WorkletMapField
    {
      UpdatePositionWorklet(const Real& dt)
        : _dt(dt)
      {
      }
    
      using ControlSignature = void(FieldIn velocity, ExecObject locator,FieldInOut position);
      using ExecutionSignature = void(_1, _2, _3);
    
      template<typename CoordType>
      VTKM_EXEC void operator()(const Vec3f& velocity,
                                const ExecPointLocator locator,
                                CoordType& position) const
      {
        position += velocity * _dt;
        locator.UpdateOverRangePoint(position);
      }
      Real _dt;
    };

    struct UpdatePositionFlagWorklet : vtkm::worklet::WorkletMapField
    {
      UpdatePositionFlagWorklet(const Real& dt)
        : _dt(dt)
      {
      }

      using ControlSignature = void(FieldIn velocity,
                                    ExecObject locator,
                                    FieldInOut position,
                                    FieldInOut position_flag);
      using ExecutionSignature = void(_1, _2, _3, _4);

      template<typename CoordType>
      VTKM_EXEC void operator()(const vtkm::Vec3f& velocity,
                                const ExecPointLocator locator,
                                CoordType& position,
                                Id3& position_flag) const
      {
        position += velocity * _dt;
        //locator.UpdateFlagOverRangePoint(position, position_flag);
      }
      Real _dt;
    };
    
    struct ComputeNewFarElectrostaticsWorklet : vtkm ::worklet::WorkletMapField
    {
      ComputeNewFarElectrostaticsWorklet(const IdComponent& Kmax)
        : _Kmax(Kmax)
      {
      }

      using ControlSignature = void(FieldIn atoms_id,
                                    WholeArrayIn whole_rhok,
                                    ExecObject force_function,
                                    ExecObject topology,
                                    ExecObject locator,
                                    FieldOut force);
      using ExecutionSignature = void(_1, _2, _3, _4, _5, _6);

      template<typename WholeRhokType>
      VTKM_EXEC void operator()(const Id atoms_id,
                                const WholeRhokType& whole_rhok,
                                ExecForceFunction& force_function,
                                const ExecTopology& topology,
                                const ExecPointLocator& locator,
                                Vec3f& force) const
      {
         force = { 0, 0, 0 };
         const auto& charge_p_i = topology.GetCharge(atoms_id);
         const auto& p_i = locator.GetPtsPosition(atoms_id);
         auto function = [&](const vtkm::Vec3f& M, const vtkm::Id& indexEwald)
         {
           auto rhok_ri = whole_rhok.Get(indexEwald - 1);
           force += force_function.ComputeFarEle(M, p_i, charge_p_i, rhok_ri);
         };
         locator.ExecuteOnKNeighbor(_Kmax, atoms_id, function);
      }
      IdComponent _Kmax;
    };   

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

    void ComputeRBLNeighboursOnce(const Id& rs_num,
                                  const Id& pice_num,
                                  const vtkm::cont::ArrayHandle<vtkm::Id>& atoms_id,
                                  const ContPointLocator& locator,
                                  GroupVecType& id_verletlist_group,
                                  GroupNumType& num_verletlist,
                                  CoordOffsetType& offset_verletlist_group)
    {
      vtkm::cont::Invoker{}(ComputeRBLNeighboursOnceWorklet{ rs_num, pice_num },
                            atoms_id,
                            locator,
                            id_verletlist_group,
                            num_verletlist,
                            offset_verletlist_group);
    }

    void NearForceVerletERF(const Real& cut_off,
                            const vtkm::cont::ArrayHandle<vtkm::Id>& atoms_id,
                            const ContPointLocator& locator,
                            const ContTopology& topology,
                            const ContForceFunction& force_function,
                            const ContStaticTable& static_table,
                            const GroupVecType& Group_j,
                            const vtkm::cont::ArrayHandle<vtkm::Id>& num_j,
                            const CoordOffsetType& coord_offset_j,
                            vtkm::cont::ArrayHandle<Vec3f>& LJforce)
    {
      vtkm::cont::Invoker{}(ComputeNearForceVerletERFWorklet{ cut_off },
                            atoms_id,
                            locator,
                            topology,
                            force_function,
                            static_table,
                            Group_j,
                            num_j,
                            coord_offset_j,
                            LJforce);
    }

    void NearForceVerlet(const Real& cut_off,
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
                         vtkm::cont::ArrayHandle<Vec3f>& nearforce)
    {
      vtkm::cont::Invoker{}(ComputeNearForceVerletWorklet{ cut_off, box },
                            atoms_id,
                            locator,
                            topology,
                            force_function,
                            Group_j,
                            num_j,
                            coord_offset_j,
                            group_ids,
                            group_weights,
                            nearforce);
    }
    void LJForceVerlet(const Real& cut_off,
                       const Vec3f& box,
                       const vtkm::cont::ArrayHandle<vtkm::Id>& atoms_id,
                       const ContPointLocator& locator,
                       const ContTopology& topology,
                       const ContForceFunction& force_function,
                       const GroupVecType& Group_j,
                       const vtkm::cont::ArrayHandle<vtkm::Id>& num_j,
                       const CoordOffsetType& coord_offset_j,
                       vtkm::cont::ArrayHandle<Vec3f>& LJforce)
    {
      vtkm::cont::Invoker{}(ComputeLJForceVerletWorklet{ cut_off, box },
                            atoms_id,
                            locator,
                            topology,
                            force_function,
                            Group_j,
                            num_j,
                            coord_offset_j,
                            LJforce);
    }

    void ComputeThermo(const Real& cut_off,
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
                       vtkm::cont::ArrayHandle<Real>& LJPE)
    {
      vtkm::cont::Invoker{}(ComputeThermoWorklet{ cut_off, box },
                            atoms_id,
                            locator,
                            topology,
                            force_function,
                            Group_j,
                            num_j,
                            coord_offset_j,
                            LJforce,
                            LJvirial,
                            LJPE);
    }

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
                     vtkm::cont::ArrayHandle<Real>& EAM_fp)
    {
      vtkm::cont::Invoker{}(ComputeEAMfpVerletWorklet{ cut_off, box },
                            atoms_id,
                            rhor_spline,
                            frho_spline,
                            locator,
                            force_function,
                            Group_j,
                            num_j,
                            coord_offset_j,
                            EAM_fp);
    }

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
                        vtkm::cont::ArrayHandle<Vec3f>& force)
    {
      vtkm::cont::Invoker{}(ComputeEAMForceVerletWorklet{ cut_off, box },
                            atoms_id,
                            rhor_spline,
                            z2r_spline,
                            EAM_fp,
                            locator,
                            force_function,
                            Group_j,
                            num_j,
                            coord_offset_j,
                            force);
    }

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

    void EAM_fp(const vtkm::cont::ArrayHandle<vtkm::Id>& atoms_id,
                const vtkm::cont::ArrayHandle<Real>& EAM_rho,
                const vtkm::cont::ArrayHandle<Vec7f>& frho_spline,
                const ContPointLocator& locator,
                const ContTopology& topology,
                const ContForceFunction& force_function,
                vtkm::cont::ArrayHandle<Real>& fp)
    {
      vtkm::cont::Invoker{}(ComputeEAMfpWorklet{},
                            atoms_id,
                            EAM_rho,
                            frho_spline,
                            locator,
                            topology,
                            force_function,
                            fp);
    }

    void EAM_force(const Real& eam_cut_off,
                   const Vec3f& box,
                   const vtkm::cont::ArrayHandle<vtkm::Id>& atoms_id,
                   const vtkm::cont::ArrayHandle<Real>& fp,
                   const vtkm::cont::ArrayHandle<Vec7f>& rhor_spline,
                   const vtkm::cont::ArrayHandle<Vec7f>& z2r_spline,
                   const ContPointLocator& locator,
                   const ContTopology& topology,
                   const ContForceFunction& force_function,
                   vtkm::cont::ArrayHandle<Vec3f>& eam_force)
    {

      vtkm::cont::Invoker{}(ComputeEAMforceWorklet{ eam_cut_off, box },
                            atoms_id,
                            fp,
                            rhor_spline,
                            z2r_spline,
                            locator,
                            topology,
                            force_function,
                            eam_force);
    }

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
                         vtkm::cont::ArrayHandle<Vec3f>& corr_force)
    {
      vtkm::cont::Invoker{}(ComputeNearForceRBLERFWorklet{ rs_num, pice_num, qqr2e, box },
                            atoms_id,
                            locator,
                            topology,
                            force_function,
                            static_table,
                            id_verletlist_group,
                            num_verletlist,
                            offset_verletlist_group,
                            corr_force);
    }

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
                                     vtkm::cont::ArrayHandle<Vec3f>& corr_force)
    {
      vtkm::cont::Invoker{}(
        ComputeNearForceRBLERFSpecialBondsWorklet{ rs_num, pice_num, qqr2e, box },
        atoms_id,
        locator,
        topology,
        force_function,
        static_table,
        id_verletlist_group,
        num_verletlist,
        offset_verletlist_group,
        group_ids,
        group_weights,
        corr_force);
    }

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
                       vtkm::cont::ArrayHandle<Vec3f>& corr_force)
     {
      vtkm::cont::Invoker{}(ComputeNearForceRBLWorklet{ rs_num, pice_num, qqr2e },
                            atoms_id,
                            locator,
                            topology,
                            force_function,
                            id_verletlist_group,
                            num_verletlist,
                            offset_verletlist_group,
                            corr_force);
     }

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
                     vtkm::cont::ArrayHandle<Vec3f>& corr_ljforce)
     {
      vtkm::cont::Invoker{}(ComputeLJForceRBLWorklet{ rs_num, pice_num, box },
                            atoms_id,
                            locator,
                            topology,
                            force_function,
                            id_verletlist_group,
                            num_verletlist,
                            offset_verletlist_group,
                            corr_ljforce);
     }

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
                vtkm::cont::ArrayHandle<Real>& EAM_fp)
     {
      vtkm::cont::Invoker{}(ComputeEAMfp_Worklet{ rc, box, rs, rs_num, pice_num },
                            atoms_id,
                            rhor_spline,
                            frho_spline,
                            locator,
                            force_function,
                            id_verletlist_group,
                            num_verletlist_group,
                            offset_verletlist_group,
                            EAM_fp);
     }

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
                    vtkm::cont::ArrayHandle<Vec3f>& corr_force)
     {
      vtkm::cont::Invoker{}(ComputeEAMRBLWorklet{ rc, box, rs, rs_num, pice_num },
                            atoms_id,
                            rhor_spline,
                            z2r_spline,
                            EAM_fp,
                            locator,
                            force_function,
                            id_verletlist_group,
                            num_verletlist_group,
                            offset_verletlist_group,
                            corr_force);
     }

    void SumRBLCorrForce(const vtkm::Vec3f corr_value,
                         const vtkm::cont::ArrayHandle<vtkm::Vec3f>& corr_force,
                         vtkm::cont::ArrayHandle<vtkm::Vec3f>& LJforce)
    {
         vtkm::cont::Invoker{}(SumRBLCorrForceWorklet{corr_value}, corr_force, LJforce);
    }

    void LJForceWithPeriodicBC(const Real& cut_off,
                               const vtkm::cont::ArrayHandle<vtkm::Id>& atoms_id,
                               const ContPointLocator& locator,
                               const ContTopology& topology,
                               const ContForceFunction& force_function,
                               vtkm::cont::ArrayHandle<Vec3f>& LJforce)
    {
      vtkm::cont::Invoker{}(ComputeLJForceWithPeriodicBCWorklet{ cut_off },
                            atoms_id,
                            locator,
                            topology,
                            force_function,
                            LJforce);
    }

    void SpecicalBondsLJForce(const Real& cut_off,
                               const vtkm::cont::ArrayHandle<vtkm::Id>& atoms_id,
                               const ContPointLocator& locator,
                               const ContTopology& topology,
                               const ContForceFunction& force_function,
                               const GroupIdIdType& group_ids,
                               const GroupRealIdType& group_weights, 
                               vtkm::cont::ArrayHandle<Vec3f>& LJforce)
    {
      vtkm::cont::Invoker{}(ComputeSpecialBondsLJWorklet{ cut_off },
                            atoms_id,
                            locator,
                            topology,
                            force_function,
                            group_ids,
                            group_weights,
                            LJforce);
    }

    void ComputeNearElectrostatics(const vtkm::cont::ArrayHandle<vtkm::Id>& atoms_id,
                                   const ContPointLocator& locator,
                                   const ContTopology& topology,
                                   const ContForceFunction& force_function,
                                   vtkm::cont::ArrayHandle<vtkm::Vec3f>& eleNearforce)
    {
      vtkm::cont::Invoker{}(ComputeNearElectrostaticsWorklet{},
                            atoms_id,
                            locator,
                            topology,
                            force_function,
                            eleNearforce);
    }

    void ComputeNearElectrostaticsERF(const vtkm::cont::ArrayHandle<vtkm::Id>& atoms_id,
                                      const ContStaticTable& static_table,
                                      const ContPointLocator& locator,
                                      const ContTopology& topology,
                                      const ContForceFunction& force_function,
                                      vtkm::cont::ArrayHandle<vtkm::Vec3f>& eleNearforce)
    {
      vtkm::cont::Invoker{}(ComputeNearElectrostaticsERFWorklet{},
                            atoms_id,
                            static_table,
                            locator,
                            topology,
                            force_function,
                            eleNearforce);
    }

    void ComputeSpecialCoul(const Vec3f& box,
                            const vtkm::cont::ArrayHandle<vtkm::Id>& atoms_id,
                            const GroupVecType& group_vec,
                            const ContForceFunction& force_function,
                            const ContTopology& topology,
                            const ContPointLocator& locator,
                            vtkm::cont::ArrayHandle<vtkm::Vec3f>& SpecCoulforce,
                            vtkm::cont::ArrayHandle<Vec6f>& SpecCoulVirial)
    {
      vtkm::cont::Invoker{}(ComputeSpecialCoulWorklet{ box },
                            atoms_id,
                            group_vec,
                            force_function,
                            topology,
                            locator,
                            SpecCoulforce,
                            SpecCoulVirial);
    }

    void ComputeSpecialCoulGeneral(const Vec3f& box,
                            const vtkm::cont::ArrayHandle<vtkm::Id>& atoms_id,
                            const GroupVecType& group_vec,
                            const ContForceFunction& force_function,
                            const ContTopology& topology,
                            const ContPointLocator& locator,
                            const GroupIdIdType& group_ids,
                            const GroupRealIdType& group_weights,
                            vtkm::cont::ArrayHandle<vtkm::Vec3f>& SpecCoulforce,
                            vtkm::cont::ArrayHandle<Vec6f>& SpecCoulVirial)
    {
      vtkm::cont::Invoker{}(ComputeSpecialCoulGeneralWorklet{ box },
                            atoms_id,
                            group_vec,
                            force_function,
                            topology,
                            locator,
                            group_ids,
                            group_weights,
                            SpecCoulforce,
                            SpecCoulVirial);
    }

    void ComputeNewRBEForce(const IdComponent& rbeP,
                            const vtkm::cont::ArrayHandle<vtkm::Id>& atoms_id,
                            const vtkm::cont::ArrayHandle<vtkm::Vec3f>& psample,
                            const vtkm::cont::ArrayHandle<vtkm::Vec2f>& whole_rhokRBE,
                            const ContForceFunction& force_function,
                            const ContTopology& topology,
                            const ContPointLocator& locator,
                            vtkm::cont::ArrayHandle<vtkm::Vec3f>& eleRBENewforce)
    {
      vtkm::cont::Invoker{}(ComputeNewRBEForceWorklet{ rbeP },
                            atoms_id,
                            psample,
                            whole_rhokRBE,
                            force_function,
                            topology,
                            locator,
                            eleRBENewforce);
    }

    void InitPosition(const vtkm::cont::UnknownCellSet& cellset,
                      const vtkm::cont::CoordinateSystem& coord,
                      vtkm::cont::ArrayHandle<vtkm::Vec3f>& position)
    {
        vtkm::cont::Invoker{}(InitPositionWorklet{}, cellset, coord, position);
    }
    
    void InitCondtion(const vtkm::cont::ArrayHandle<vtkm::Vec3f>& position,
                      const vtkm::cont::ArrayHandle<vtkm::Vec3f>& v_array,
                      vtkm::cont::ArrayHandle<Real>& mass,
                      vtkm::cont::ArrayHandle<vtkm::Vec3f>& velocity)
    {
      vtkm::cont::Invoker{}(InitialCondtionWorklet{}, position, v_array, mass, velocity);
    }
    
    void SumFarNearForce(const vtkm::cont::ArrayHandle<vtkm::Vec3f>& eleFarNewforce,
                         const vtkm::cont::ArrayHandle<vtkm::Vec3f>& nearforce,
                         vtkm::cont::ArrayHandle<vtkm::Vec3f>& Allforce)
    {
      vtkm::cont::Invoker{}(SumFarNearForceWorklet{}, eleFarNewforce, nearforce, Allforce);
    }

    void SumVirial(const vtkm::cont::ArrayHandle<Vec6f>& LJVirial,
                   const vtkm::cont::ArrayHandle<Vec6f>& EwaldVirial,
                   vtkm::cont::ArrayHandle<Vec6f>& AllVirial)
    {
      vtkm::cont::Invoker{}(SumVirialWorklet{}, LJVirial, EwaldVirial, AllVirial);
    }   


    void SumFarNearLJForce(const vtkm::cont::ArrayHandle<vtkm::Vec3f>& eleFarNewforce,
                       const vtkm::cont::ArrayHandle<vtkm::Vec3f>& eleNearforce,
                       const vtkm::cont::ArrayHandle<vtkm::Vec3f>& LJforce,
                       vtkm::cont::ArrayHandle<vtkm::Vec3f>& Allforce)
    {
         vtkm::cont::Invoker{}(SumFarNearLJForceWorklet{}, eleFarNewforce, eleNearforce, LJforce, Allforce);
    }    

    void UnderdampedLangevin(const Vec3f& gaussian,
                             const Real& kBT,
                             const Real& gamma,
                             const Real& dt,
                             const vtkm::cont::ArrayHandle<Real>& mass,
                             const vtkm::cont::ArrayHandle<vtkm::Vec3f>& velocity,
                             vtkm::cont::ArrayHandle<vtkm::Vec3f>& Allforce)
{
      vtkm::cont::Invoker{}(
        UnderdampedLangevinWorklet{ gaussian, kBT, gamma, dt }, mass,velocity, Allforce);
}
    
    void UpdateVelocity(const Real& dt,
                        const Real& unit_factor,
                        const vtkm::cont::ArrayHandle<vtkm::Vec3f>& Allforce,
                        const vtkm::cont::ArrayHandle<Real>& mass,
                        vtkm::cont::ArrayHandle<vtkm::Vec3f>& velocity)
{
      vtkm::cont::Invoker{}(UpdateVelocityWorklet{ dt, unit_factor }, Allforce, mass,velocity);
    }
    
    void ComputerKineticEnergy(const vtkm::cont::ArrayHandle<vtkm::Vec3f>& velocity,
                               const vtkm::cont::ArrayHandle<Real>& mass,
                               vtkm::cont::ArrayHandle<Real>& sq_velocity)
{
  vtkm::cont::Invoker{}(ComputerKineticEnergyWorklet{}, velocity, mass, sq_velocity);
}

    void UpdateVelocityRescale(const Real& coeff_Berendsen,
                               vtkm::cont::ArrayHandle<vtkm::Vec3f>& velocity)
{
      vtkm::cont::Invoker{}(UpdateVelocityRescaleWorklet{ coeff_Berendsen }, velocity);
} 

    void UpdateVelocityNoseHoover(const Real& dt,
                                  const Real& unit_factor,
                                  const Real& nosehooverxi,                                  
                                  const vtkm::cont::ArrayHandle<vtkm::Vec3f>& Allforce,
                                  const vtkm::cont::ArrayHandle<Real>& mass,
                                  vtkm::cont::ArrayHandle<vtkm::Vec3f>& velocity)
{
  vtkm::cont::Invoker{}(
        UpdateVelocityNoseHooverWorklet{ dt, unit_factor, nosehooverxi }, Allforce, mass, velocity);
}
    
    void UpdatePosition(const Real& dt,
                        const vtkm::cont::ArrayHandle<vtkm::Vec3f>& velocity,
                        const ContPointLocator& locator,
                        vtkm::cont::ArrayHandle<vtkm::Vec3f>& position)
    {
         vtkm::cont::Invoker{}(
         UpdatePositionWorklet{ dt}, velocity, locator, position);
    }  

    void UpdatePositionFlag(const Real& dt,
                            const vtkm::cont::ArrayHandle<vtkm::Vec3f>& velocity,
                            const ContPointLocator& locator,
                            vtkm::cont::ArrayHandle<vtkm::Vec3f>& position,
                            vtkm::cont::ArrayHandle<vtkm::Id3>& position_flag)
    {
         vtkm::cont::Invoker{}(UpdatePositionFlagWorklet{ dt }, velocity, locator, position, position_flag);
    }   
    
    void ComputeChargeStructureFactorComponent(const vtkm::Vec3f& kl,
                                               const vtkm::cont::ArrayHandle<vtkm::Vec3f>& position,
                                               const vtkm::cont::ArrayHandle<Real>& charge,
                                               vtkm::cont::ArrayHandle<Real>& Density_Real,
                                               vtkm::cont::ArrayHandle<Real>& Density_Image)
    {
         vtkm::cont::Invoker{}(ComputeChargeStructureFactorComponentWorklet{ kl },
                               position,
                               charge,
                               Density_Real,
                               Density_Image);
    }

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
/*                                             GroupVecType& psamplekey_group,
                                             GroupRealType& rhok_Re_group,
                                             GroupRealType& rhok_Im_group*/)
    {
      vtkm::cont::Invoker{}(ComputePnumberChargeStructureFactorWorklet{ _box, pnumber,pos_number },
                            atoms_id,
                            position,
                            charge,
                            psample,
                            //psamplekey_group,
                            rhok_Re_group,
                            rhok_Im_group);
    }

    void ChangePnumberChargeStructureFactor(const vtkm::cont::ArrayHandle<Real>& rhok_Re,
                                             const vtkm::cont::ArrayHandle<Real>& rhok_Im,
                                             vtkm::cont::ArrayHandle<vtkm::Vec2f>& new_rhok)
    {
      vtkm::cont::Invoker{}(
        ChangePnumberChargeStructureFactorRBEWorklet{}, rhok_Re, rhok_Im, new_rhok);
    }
    

    void ComputeNewFarElectrostatics(const IdComponent& Kmax,
                                     const vtkm::cont::ArrayHandle<vtkm::Id>& atoms_id,
                                     const vtkm::cont::ArrayHandle<vtkm::Vec2f>& whole_rhok,
                                     const ContForceFunction& force_function,
                                     const ContTopology& topology,
                                     const ContPointLocator& locator,
                                     vtkm::cont::ArrayHandle<vtkm::Vec3f>& eleFarNewforce)
    {
         vtkm::cont::Invoker{}(ComputeNewFarElectrostaticsWorklet{ Kmax },
                               atoms_id,
                               whole_rhok,
                               force_function,
                               topology,
                               locator,
                               eleFarNewforce);
    }

      void ComputeLJVirial(const Real& cut_off,
                       const Vec3f& box,
                       const vtkm::cont::ArrayHandle<vtkm::Id>& atoms_id,
                       const ContPointLocator& locator,
                       const ContTopology& topology,
                       const ContForceFunction& force_function,
                       const GroupVecType& Group_j,
                       const vtkm::cont::ArrayHandle<vtkm::Id>& num_j,
                       const CoordOffsetType& coord_offset_j,
                       vtkm::cont::ArrayHandle<Vec6f>& LJVirial)
    {
         vtkm::cont::Invoker{}(ComputeLJVirialWorklet{ cut_off, box },
                               atoms_id,
                               locator,
                               topology,
                               force_function,
                               Group_j,
                               num_j,
                               coord_offset_j,
                               LJVirial);
    }

      
    void ComputeCoulVirial(const Real& cut_off,
                         const Vec3f& box,
                         const vtkm::cont::ArrayHandle<vtkm::Id>& atoms_id,
                         const ContPointLocator& locator,
                         const ContTopology& topology,
                         const ContForceFunction& force_function,
                         const GroupVecType& Group_j,
                         const vtkm::cont::ArrayHandle<vtkm::Id>& num_j,
                         const CoordOffsetType& coord_offset_j,
                          vtkm::cont::ArrayHandle<Vec6f>& CoulVirial)
    {
         vtkm::cont::Invoker{}(ComputeCoulVirialWorklet{ cut_off, box },
                               atoms_id,
                               locator,
                               topology,
                               force_function,
                               Group_j,
                               num_j,
                               coord_offset_j,
                               CoulVirial);
    }

    
     void ComputeEwaldVirial(const IdComponent& Kmax,
                            const Real& unit_factor,
                            const vtkm::cont::ArrayHandle<vtkm::Id>& atoms_id,
                            const vtkm::cont::ArrayHandle<vtkm::Vec2f>& whole_rhok,
                            const ContForceFunction& force_function,
                            const ContTopology& topology,
                            const ContPointLocator& locator,
                            vtkm::cont::ArrayHandle<Vec6f>& EwaldVirial)
    {
         vtkm::cont::Invoker{}(ComputeEwaldVirialWorklet{ Kmax, unit_factor },
                               atoms_id,
                               whole_rhok,
                               force_function,
                               topology,
                               locator,
                               EwaldVirial);
    }

     void fix_press_berendsen(const Real& scale_factor,
                             vtkm::cont::ArrayHandle<vtkm::Vec3f>& position,
                             const ContPointLocator& locator)
     {
         vtkm::cont::Invoker{}(fix_press_berendsenWorklet{ scale_factor }, position, locator);
     }  

      void ApplyPbcFlag(const Vec3f& box,
                    const Vec<Vec2f, 3>& range,
                   vtkm::cont::ArrayHandle<vtkm::Vec3f>& position,
                   const ContPointLocator& locator,
                   vtkm::cont::ArrayHandle<vtkm::Id3>& position_flag)
     {
         vtkm::cont::Invoker{}(ApplyPbcFlagWorklet{ box, range }, position, locator, position_flag);
     }  

      void ApplyPbc(const Vec3f& box,
                   const Vec<Vec2f, 3>& range,
                   vtkm::cont::ArrayHandle<vtkm::Vec3f>& position,
                   const ContPointLocator& locator)
     {
         vtkm::cont::Invoker{}(ApplyPbcWorklet{ box, range }, position, locator);
     }

      void X2Lamda(const Vec6f& h_inv,
                   const vtkm::Vec<vtkm::Range, 3>& range,
                   vtkm::cont::ArrayHandle<vtkm::Vec3f>& position)
     {
         vtkm::cont::Invoker{}(X2LamdaWorklet{ h_inv, range }, position);
     } 

      void Lamda2X(const Vec6f& h,
                  const vtkm::Vec<vtkm::Range, 3>& range,
                  vtkm::cont::ArrayHandle<vtkm::Vec3f>& position)
     {
         vtkm::cont::Invoker{}(Lamda2XWorklet{ h, range }, position);
     } 
 }
