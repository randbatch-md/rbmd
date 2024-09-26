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

#pragma once
#include <vtkm/Pair.h>
#include <vtkm/Hash.h>
#include <vtkm/Types.h>
#include <vtkm/cont/ArrayHandleGroupVec.h>
#include <vtkm/cont/ArrayHandleGroupVecVariable.h>
//namespace galois{
  using vtkm::FloatDefault;
  using vtkm::HashType;
  using vtkm::Id;
  using vtkm::Id2;
  using vtkm::Id3;
  using vtkm::IdComponent;
  using vtkm::Pair;
  using vtkm::UInt16;
  using vtkm::UInt32;
  using vtkm::UInt8;
  using vtkm::Vec3f;
  using vtkm::Vec2f;

  using vtkm::Vec;
  using Real = vtkm::FloatDefault;
  using GroupVecType = vtkm::cont::ArrayHandleGroupVecVariable<vtkm::cont::ArrayHandle<vtkm::Id>,
                                                               vtkm::cont::ArrayHandle<vtkm::Id>>;

  using GroupNumType = vtkm::cont::ArrayHandleGroupVec<vtkm::cont::ArrayHandle<vtkm::Id>, 2>;

  using CoordOffsetType =
    vtkm::cont::ArrayHandleGroupVecVariable<vtkm::cont::ArrayHandle<vtkm::Vec3f>,
                                            vtkm::cont::ArrayHandle<vtkm::Id>>;
  using Vec7f = vtkm::Vec<vtkm::FloatDefault, 7>;

  using GroupRealType =
    vtkm::cont::ArrayHandleGroupVecVariable<vtkm::cont::ArrayHandle<Real>,
                                            vtkm::cont::ArrayHandle<vtkm::Id>>;
  //using RealGradient = vtkm::Vec3f;
//}
