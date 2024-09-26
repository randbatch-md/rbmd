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

#include "LJInitCondition.h"
#include "FieldName.h"
#include "UnitFactor.h"

LJInitCondition::LJInitCondition(const Configuration& cfg)
  : MeshFreeCondition(cfg)
{
  InitField();
  SetParameters();
}

void LJInitCondition::Execute() 
{

  AddMoleculeInfo();
  UpdateField();
}

void LJInitCondition::UpdateField()
{
  MeshFreeCondition::UpdateField();
}

void LJInitCondition::AddMoleculeInfo() 
{
  auto N = _position.GetNumberOfValues();
  std::vector<Real> special_weights(N, 1.0);
  std::vector<Id> special_offsets(N, 1);

  vtkm::Id offsetSize;
  vtkm::cont::ArrayHandle<vtkm::Id> offsetsArray = vtkm::cont::ConvertNumComponentsToOffsets(
    vtkm::cont::make_ArrayHandle(special_offsets), offsetSize);
  auto special_offsets_array = _para.GetFieldAsArrayHandle<Id>(field::special_offsets);
  vtkm::cont::ArrayCopy(offsetsArray, special_offsets_array);

  auto special_ids_array = _para.GetFieldAsArrayHandle<Id>(field::special_ids);
  vtkm::cont::ArrayHandleIndex atomsIdIndex(_position.GetNumberOfValues());
  vtkm::cont::ArrayCopy(atomsIdIndex, special_ids_array);

  auto special_weights_array = _para.GetFieldAsArrayHandle<Real>(field::special_weights);
  vtkm::cont::ArrayCopy(vtkm::cont::make_ArrayHandle(special_weights), special_weights_array);
}

void LJInitCondition::InitField()
{
  MeshFreeCondition::InitField();
  _para.AddField(field::pts_type, ArrayHandle<Id>{});
  _para.AddField(field::position_flag, ArrayHandle<Id3>{});
  _para.AddField(field::center_position, ArrayHandle<Vec3f>{});
  _para.AddField(field::target_position, ArrayHandle<Vec3f>{});
  _para.AddField(field::epsilon, ArrayHandle<Real>{});
  _para.AddField(field::sigma, ArrayHandle<Real>{});
}

void LJInitCondition::SetParameters()
{
  auto xLength = _x_range[1] - _x_range[0];
  auto yLength = _y_range[1] - _y_range[0];
  auto zLength = _z_range[1] - _z_range[0];
  Vec3f box = { xLength, yLength, zLength };
  auto num_pos = _dims[0] * _dims[1] * _dims[2];
  auto range = vtkm::Vec<vtkm::Range, 3>{ { _x_range[0], _x_range[1] },
                                          { _y_range[0], _y_range[1] },
                                          { _z_range[0], _z_range[1] } };
  _para.SetParameter(PARA_VLENGTH, xLength);
  _para.SetParameter(PARA_BOX, box);
  _para.SetParameter(PARA_VOLUME, box[0] * box[1] * box[2]);
  _para.SetParameter(PARA_RHO, num_pos / (box[0] * box[1] * box[2]));
  _para.SetParameter(PARA_RANGE, range);
  _para.SetParameter(gtest::velocity_type, Get<std::string>("velocity_type"));
  _para.SetParameter(PARA_UNIT, Get<std::string>("unit"));
  auto cut_off = _para.GetParameter<Real>(PARA_CUTOFF);
  auto bin_number = Id3{
    static_cast<int>(xLength / cut_off),
    static_cast<int>(yLength / cut_off),
    static_cast<int>(zLength / cut_off),
  };
  _para.SetParameter(PARA_BIN_NUMBER, bin_number);
}
