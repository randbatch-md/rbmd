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
//  The post-processing data produced by VASPKIT may be used in any publications 
//  provided that its use is explicitly acknowledged. A suitable reference for VASPKIT is:
//  [1] Gao W, Zhao T, Guo Y, et al.RBMD: A molecular dynamics package enabling to simulate 
//  10 million all - atom particles in a single graphics processing unit[J].arXiv preprint arXiv : 2407.09315, 2024.
// 
//  Contact Email : [support_wz@sciai.com.cn]
//==================================================================================

#include "MadelungInitCondition.h"
#include "FieldName.h"
#include <vtkm/cont/DataSetBuilderUniform.h>
#include <vtkm/worklet/WorkletPointNeighborhood.h>
#include <vtkm/cont/Invoker.h>
#include <vtkm/cont/ArrayCopy.h>

struct MadelungInitCondition::SetChargeWorklet : vtkm::worklet::WorkletPointNeighborhood
{
  using ControlSignature = void(CellSetIn, FieldOut charge);
  using ExecutionSignature = void(Boundary, _2);

  template<typename T>
  VTKM_EXEC void operator()(const vtkm::exec::BoundaryState& boundary, T& charge) const
  {
    if ((boundary.IJK[0] + boundary.IJK[1] + boundary.IJK[2]) % 2 == 0)
      charge = -1;
    else
    {
      charge = 1;
    }
  }
};

MadelungInitCondition::MadelungInitCondition(const Configuration& cfg)
  : InitCondition(cfg)
  , _dims(GetVectorValue<int>("dims"))
  , _x_range(GetVectorValue<Real>("x_range"))
  , _y_range(GetVectorValue<Real>("y_range"))
  , _z_range(GetVectorValue<Real>("z_range"))
{

}

void MadelungInitCondition::Execute() 
{
  UpdateField();
}

void MadelungInitCondition::UpdateField() 
{
  //update position
  auto position = _para.GetFieldAsArrayHandle<Vec3f>(field::position);
  Vec3f start{ _x_range[0], _y_range[0], _z_range[0] };
  Vec3f end{ _x_range[1], _y_range[1], _z_range[1] };
  Vec3f space = end - start;
  Id3 point_dim = { _dims[0], _dims[1], _dims[2]};
  for (size_t i = 0; i < 3; i++)
  {
    space[i] = (end[i] - start[i]) / _dims[i];
  }
  auto data_set = vtkm::cont::DataSetBuilderUniform::Create(point_dim, start, space);
  vtkm::cont::ArrayCopy(data_set.GetCoordinateSystem().GetData(), position);
  
  //update charge
  auto charge = _para.GetFieldAsArrayHandle<Real>(field::charge);
  vtkm::cont::Invoker{}(SetChargeWorklet{}, data_set.GetCellSet(),charge);
}
