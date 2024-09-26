//==================================================================================
//  RBMD 1.0 is developed for random batch molecular dynamics calculation.
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
//  Contact Email : [your - email@example.com]
//==================================================================================

#include "locator/ContPointLocator.h"
#include "Types.h"
#include <vtkm/cont/Algorithm.h>
#include <vtkm/cont/ArrayCopy.h>
#include <vtkm/cont/ArrayHandleIndex.h>
#include <vtkm/cont/Invoker.h>
#include <vtkm/worklet/WorkletMapField.h>

namespace internal
{

class BinPointsWorklet : public vtkm::worklet::WorkletMapField
{
public:
  using ControlSignature = void(FieldIn coord, FieldOut label);

  using ExecutionSignature = void(_1, _2);

  VTKM_CONT
  BinPointsWorklet(vtkm::Vec3f min, vtkm::Vec3f max, vtkm::Id3 dims)
    : _min(min)
    , _dims(dims)
    , _dxdydz((max - _min) / _dims)
  {
  }

  template <typename CoordVecType, typename IdType>
  VTKM_EXEC void operator()(const CoordVecType& coord, IdType& label) const
  {
    vtkm::Id3 ijk = (coord - _min) / _dxdydz;
    ijk = vtkm::Max(ijk, vtkm::Id3(0));
    ijk = vtkm::Min(ijk, this->_dims - vtkm::Id3(1));
    label = ijk[0] + ijk[1] * _dims[0] + ijk[2] * _dims[0] * _dims[1];
  }

private:
  vtkm::Vec3f _min;
  vtkm::Id3 _dims;
  vtkm::Vec3f _dxdydz;
};

} // vtkm::cont::internal

void ContPointLocator::Build()
{
  VTKM_LOG_SCOPE(vtkm::cont::LogLevel::Perf, "PointLocatorSparseGrid::Build");

  auto rmin = static_cast<vtkm::Vec3f>(this->_left_bottom);
  auto rmax = static_cast<vtkm::Vec3f>(this->_right_top);

  using internal::BinPointsWorklet;
  vtkm::cont::ArrayHandle<vtkm::Id> cellIds;
  BinPointsWorklet cellIdWorklet(rmin, rmax, this->_dims);
  vtkm::cont::Invoker invoke;
  invoke(cellIdWorklet, this->_position, cellIds);


  // Group points of the same cell together by sorting them according to the cell ids
  vtkm::cont::Algorithm::SortByKey(cellIds, this->_point_ids);

  // for each cell, find the lower and upper bound of indices to the sorted point ids.
  vtkm::cont::ArrayHandleCounting<vtkm::Id> cell_ids_counting(
    0, 1, this->_dims[0] * this->_dims[1] * this->_dims[2]);
  vtkm::cont::Algorithm::UpperBounds(cellIds, cell_ids_counting, this->_cell_upper);
  vtkm::cont::Algorithm::LowerBounds(cellIds, cell_ids_counting, this->_cell_lower);
}

  ExecPointLocator ContPointLocator::PrepareForExecution(
  vtkm::cont::DeviceAdapterId device,
  vtkm::cont::Token& token)const
{
  auto rmin = static_cast<vtkm::Vec3f>(this->_left_bottom);
  auto rmax = static_cast<vtkm::Vec3f>(this->_right_top);
  return ExecPointLocator(
    rmin,
    rmax,
    this->_cut_off,    
    this->_rs,
    this->_dims,
    this->_num_cycles,
    this->_random_num,
    this->_point_ids.PrepareForInput(device, token),
    this->_position.PrepareForInput(device, token),
    this->_cell_lower.PrepareForInput(device, token),
    this->_cell_upper.PrepareForInput(device, token),
    this->_num_verletlist.PrepareForInput(device, token),
    this->_id_verletlist_group.PrepareForInput(device, token),
    this->_offset_verletlist_group.PrepareForInput(device, token));
}

