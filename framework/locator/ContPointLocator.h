#pragma once
#include <vtkm/cont/ExecutionObjectBase.h>
#include <vtkm/cont/ArrayHandleIndex.h>
#include <vtkm/cont/ArrayCopy.h>
#include "locator/ExecPointLocator.h"

  using GroupVecPortalType =  vtkm::cont::ArrayHandleGroupVecVariable<
  vtkm::cont::ArrayHandle<vtkm::Id>,
  vtkm::cont::ArrayHandle<vtkm::Id>>;
using CoordOffsetPortalType =  vtkm::cont::ArrayHandleGroupVecVariable<
  vtkm::cont::ArrayHandle<vtkm::Vec3f>,
  vtkm::cont::ArrayHandle<vtkm::Id>>;

class ContPointLocator : public vtkm::cont::ExecutionObjectBase
{

public:
  void SetRange(const vtkm::Vec3f& left_bottom, const vtkm::Vec3f& right_top)
  {
    if (_left_bottom != left_bottom || _right_top != right_top)
    {
      this->_left_bottom = left_bottom;
      this->_right_top = right_top;
    }
    _num_cycles = std::ceil(_cut_off / _rs);
    this->_dims = ((_right_top - _left_bottom) / (_cut_off / _num_cycles));
  }

  void SetCutOff(const Real& cut_off)
  {
    if (_cut_off != cut_off)
    {
      this->_cut_off = cut_off;
    }
    _num_cycles = std::ceil(_cut_off / _rs);
    this->_dims = ((_right_top - _left_bottom) / (_cut_off / _num_cycles));
  }

  void SetRs(const Real& rs)
  {
    if (_rs != rs)
    {
      this->_rs = rs;
    }
    _num_cycles = std::ceil(_cut_off / _rs);
    this->_dims = ((_right_top - _left_bottom) / (_cut_off / _num_cycles));
  }

  void SetRandomNum(const vtkm::Id& random_num) { this->_random_num = random_num; }

  void SetPosition(
    const vtkm::cont::ArrayHandle<vtkm::Vec3f>& position,
    const vtkm::cont::ArrayHandle<vtkm::Id>& point_ids = vtkm::cont::ArrayHandle<vtkm::Id>{})
  {
    this->_position = position;

    if (0 == point_ids.GetNumberOfValues())
    {
      vtkm::cont::ArrayHandleIndex pointIndex(this->_position.GetNumberOfValues());
      vtkm::cont::ArrayCopy(pointIndex, this->_point_ids);
    }
    else
    {
      this->_point_ids = point_ids;
    }

    this->Build();
  }

  void SetVerletListInfo(const vtkm::cont::ArrayHandle<vtkm::Id>& num_verletlist,
                         const GroupVecPortalType& id_verletlist_group,
                         const CoordOffsetPortalType& offset_verletlist_group) 
  {
    _num_verletlist = num_verletlist;
    _id_verletlist_group = id_verletlist_group;
    _offset_verletlist_group = offset_verletlist_group;
  }

  VTKM_CONT
  ExecPointLocator PrepareForExecution(vtkm::cont::DeviceAdapterId device,
                                                   vtkm::cont::Token& token) const;
  
private:
  VTKM_CONT void Build();

  vtkm::Vec3f _left_bottom = { 0, 0, 0 };
  vtkm::Vec3f _right_top = { 10, 10, 10 };
  Real _cut_off=10;
  Real _rs = 5;
  Id _num_cycles;
  Id _random_num;
  vtkm::Id3 _dims; 
  vtkm::cont::ArrayHandle<vtkm::Id> _point_ids;
  vtkm::cont::ArrayHandle<vtkm::Id> _cell_lower;
  vtkm::cont::ArrayHandle<vtkm::Id> _cell_upper;
  vtkm::cont::ArrayHandle<vtkm::Vec3f> _position;
  //verletlist
  vtkm::cont::ArrayHandle<vtkm::Id> _num_verletlist;
  GroupVecPortalType _id_verletlist_group;
  CoordOffsetPortalType _offset_verletlist_group;
};
