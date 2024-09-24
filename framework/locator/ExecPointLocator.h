#pragma once
#include <vtkm/cont/CoordinateSystem.h>
#include <vtkm/Deprecated.h>
#include <vtkm/Pair.h>
#include <vtkm/VectorAnalysis.h>
#include <vtkm/cont/ArrayHandleGroupVecVariable.h>
#include <vtkm/Math.h>
#include <vtkm/Algorithms.h>
#include "Types.h"
#include "math/Math.h"

class ExecPointLocator
{
public:
  using IdPortalType = typename vtkm::cont::ArrayHandle<vtkm::Id>::ReadPortalType;
  using Vec3fPortalType = typename vtkm::cont::ArrayHandle<vtkm::Vec3f>::ReadPortalType;

  using GroupVecPortalType = typename vtkm::cont::ArrayHandleGroupVecVariable<vtkm::cont::ArrayHandle<vtkm::Id>,
    vtkm::cont::ArrayHandle<vtkm::Id>>::ReadPortalType;
  using CoordOffsetPortalType =typename vtkm::cont::ArrayHandleGroupVecVariable<vtkm::cont::ArrayHandle<vtkm::Vec3f>,
    vtkm::cont::ArrayHandle<vtkm::Id>>::ReadPortalType;

  ExecPointLocator(const vtkm::Vec3f& min,
                   const vtkm::Vec3f& max,
                   const Real& cut_off,
                   const Real& rs,
                   const vtkm::Id3& dims,
                   const vtkm::Id& num_cycles,
                   const vtkm::Id& random_num,
                   const IdPortalType& point_ids,
                   const Vec3fPortalType& position,
                   const IdPortalType& cell_lower,
                   const IdPortalType& cell_upper,
                   const IdPortalType& num_verletlist,
                   const GroupVecPortalType& id_verletlist_group,
                   const CoordOffsetPortalType& offset_verletlist_group)
    : _min(min)
    , _max(max)
    , _cut_off(cut_off)
    , _rs(rs)
    , _dims(dims)
    , _num_cycles(num_cycles)
    , _random_num(random_num)
    , _dxdydz((max - _min) / _dims)
    , _point_ids(point_ids)
    , _position(position)
    , _cell_lower(cell_lower)
    , _cell_upper(cell_upper)
    , _num_verletlist(num_verletlist)
    , _id_verletlist_group(id_verletlist_group)
    , _offset_verletlist_group(offset_verletlist_group)
  {
  }

  template<typename Func>
  VTKM_EXEC void New_ExecuteOnNeighbor(const vtkm::Id& atoms_id, Func& function) const
  {
    auto p_i = _position.Get(atoms_id);
    auto num_pts = _num_verletlist.Get(atoms_id);
    auto id_verletlist = _id_verletlist_group.Get(atoms_id);
    auto offset_verletlist = _offset_verletlist_group.Get(atoms_id);
    //printf("num_pts = %d\n", num_pts);
    for (Id p = 0; p < num_pts; p++)
    {
      auto pts_id_j = id_verletlist[p];
      auto p_j = _position.Get(pts_id_j) - offset_verletlist[pts_id_j];

      function(p_i, p_j, pts_id_j);
    }
  }


  template<typename Func>
  VTKM_EXEC void ExecuteOnNeighbor(const vtkm::Id& atoms_id, Func& function) const
  {
    auto p_i = _position.Get(atoms_id);
    vtkm::Id3 p_i_cell = InWhichCell(p_i);
    for (Id i = -_num_cycles; i <= _num_cycles; i++)
    {
      for (Id j = -_num_cycles; j <= _num_cycles; j++)
      {
        for (Id k = -_num_cycles; k <= _num_cycles; k++)
        {
          auto neighbor_ijk = p_i_cell + Id3{ i, j, k };
          auto ijk = PeriodicIndexOffset(neighbor_ijk);
          auto coord_offset = PeriodicCoordOffset(ijk - neighbor_ijk);
          auto num_pts = NumberPtsInCell(ijk);
          for (Id p = 0; p < num_pts; p++)
          {
            auto pts_id_j = PtsInCell(ijk, p);
            auto p_j = _position.Get(pts_id_j) - coord_offset;

            function(p_i, p_j, pts_id_j);
          }
        }
      }
    }
  }

  template<typename Func>
  VTKM_EXEC void ExecuteOnKNeighbor(const IdComponent& k_maxconst, const vtkm::Id& atoms_id,Func& function) const
  {
    vtkm::Id indexEwald = 0;
    for (vtkm::Id i = -k_maxconst; i <= k_maxconst; i++)
    {
      for (vtkm::Id j = -k_maxconst; j <= k_maxconst; j++)
      {
        for (vtkm::Id k = -k_maxconst; k <= k_maxconst; k++)
        {
          if (i != 0 || j != 0 || k != 0)
          {
            indexEwald++;
            vtkm::Vec3f M = { Real(i), Real(j), Real(k) };
            function(M, indexEwald);
          }
        }
      }
    }
  }
             
  // 周期性算最小距离
  VTKM_EXEC Vec3f MinDistance(const Vec3f& p1, const Vec3f& p2, const Real& _vlength) const
  {
    Vec3f vec = p1 - p2;
    // 处理原子间的距离
    for (vtkm::Id i = 0; i < 3; ++i)
    {
      while (vtkm::Abs(vec[i]) > _vlength / 2.0)
      {
        if (vec[i] < 0)
        {
          vec[i] += _vlength;
        }
        else
        {
          vec[i] -= _vlength;
        }
      }
    }
    return vec;
  }

  VTKM_EXEC Vec3f MinDistanceIf(const Vec3f& p1, const Vec3f& p2, const Real& _vlength) const
  {
    Vec3f vec = p1 - p2;
    // 处理原子间的距离
    for (vtkm::Id i = 0; i < 3; ++i)
    {
      if (vtkm::Abs(vec[i]) > _vlength / 2.0)
      {
        if (vec[i] < 0)
        {
          vec[i] += _vlength;
        }
        else
        {
          vec[i] -= _vlength;
        }
      }
    }
    return vec;
  }

  VTKM_EXEC Vec3f MinDistanceVec(const Vec3f& p1, const Vec3f& p2, const Vec3f& _box) const
  {
    Vec3f p12 = p1 - p2;
    //case1:  orthogonal box
    if (p12[0] < -_box[0] * 0.5)
    {
      p12[0] += _box[0];
    }
    else if (p12[0] > _box[0] * 0.5)
    {
      p12[0] -= _box[0];
    }

    if (p12[1] < -_box[1] * 0.5)
    {
      p12[1] += _box[1];
    }
    else if (p12[1] > _box[1] * 0.5)
    {
      p12[1] -= _box[1];
    }

    if (p12[2] < -_box[2] * 0.5)
    {
      p12[2] += _box[2];
    }
    else if (p12[2] > _box[2] * 0.5)
    {
      p12[2] -= _box[2];
    }

    //case2:  triclinic box
    return p12;
  }

  VTKM_EXEC vtkm::Vec3f GetPtsPosition(const Id& atoms_id) const
  {
    return _position.Get(atoms_id);
  }

  VTKM_EXEC vtkm::Id GetNumCycles() const { return this->_num_cycles; }

  VTKM_EXEC vtkm::Id3 InWhichCell(const vtkm::Vec3f& queryPoint) const
  {
    return (queryPoint - this->_min) / this->_dxdydz;
  }

  VTKM_EXEC vtkm::Id CellFlatId(const vtkm::Id3& ijk) const
  {
    vtkm::Id cell_id = ijk[0] + (ijk[1] * this->_dims[0]) + (ijk[2] * this->_dims[0] * this->_dims[1]);
    return cell_id;
  }

  VTKM_EXEC vtkm::Id NumberPtsInCell(const vtkm::Id3& ijk) const
  {
    auto cell_id = CellFlatId(ijk);
    vtkm::Id lower = this->_cell_lower.Get(cell_id);
    vtkm::Id upper = this->_cell_upper.Get(cell_id);

    return upper - lower;
  }

  VTKM_EXEC vtkm::Vec3f GetDxdydz() const { return this->_dxdydz; }

  VTKM_EXEC vtkm::Id PtsInCell(const vtkm::Id3& ijk, const vtkm::Id& i) const
  {
    auto cell_id = CellFlatId(ijk);
    vtkm::Id lower = this->_cell_lower.Get(cell_id);
    vtkm::Id upper = this->_cell_upper.Get(cell_id);

    vtkm::Id point_id = this->_point_ids.Get(lower + i);
    return point_id;
  }

  VTKM_EXEC vtkm::Id3 PeriodicIndexOffset(const vtkm::Id3& ijk) const
  {
    auto i = ijk[0];
    auto j = ijk[1];
    auto k = ijk[2];

    // 求余，使得i,j,k都在计算域
    i = (i + this->_dims[0]) % this->_dims[0];
    j = (j + this->_dims[1]) % this->_dims[1];
    k = (k + this->_dims[2]) % this->_dims[2];

    return { i, j, k };
  }
  VTKM_EXEC vtkm::Vec3f PeriodicCoordOffset(const vtkm::Id3& offset_ijk) const
  {
    return this->_dxdydz * offset_ijk;
  }

  VTKM_EXEC void UpdateOverRangePoint(vtkm::Vec3f& coord) const
  {
    for (int i = 0; i < 3; ++i)
    {
      if (coord[i] < _min[i])
      {
        //std::cout << "---------------------- " << std::endl;
        //std::cout << "coord: " << coord << std::endl;
        coord[i] += _max[i] - _min[i];
      }

      else if (coord[i] > _max[i])
      {
        //std::cout << "---------------------- " << std::endl;
        // std::cout << "coord: " << coord << std::endl;
        coord[i] -= _max[i] - _min[i];
      }
    }
  }

  VTKM_EXEC void UpdateFlagOverRangePoint(vtkm::Vec3f& coord, vtkm::Id3& coord_flag) const
  {
    for (int i = 0; i < 3; i++)
    {
      if (coord[i] < _min[i])
      {
        coord[i] += _max[i] - _min[i];
        coord_flag[i] -= 1;
      }

      else if (coord[i] > _max[i])
      {
        coord[i] -= _max[i] - _min[i];
        coord_flag[i] += 1;
      }
    }
  }

public:
  vtkm::FloatDefault _cut_off;
  vtkm::FloatDefault _rs;
  vtkm::Id _random_num;

private:
  vtkm::Vec3f _min;
  vtkm::Vec3f _max;
  //vtkm::FloatDefault _cut_off;
  //vtkm::FloatDefault _rs;
  vtkm::Id3 _dims;
  vtkm::Id _num_cycles;
  vtkm::Vec3f _dxdydz;

  IdPortalType _point_ids;
  IdPortalType _cell_lower;
  IdPortalType _cell_upper;
  Vec3fPortalType _position;

  //verletlist
  IdPortalType _num_verletlist;
  GroupVecPortalType _id_verletlist_group;
  CoordOffsetPortalType _offset_verletlist_group;
};
