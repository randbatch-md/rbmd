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
#include <vtkm/cont/ArrayHandle.h>
#include <vtkm/cont/CoordinateSystem.h>
#include <vtkm/cont/CellSetExplicit.h>
#include <vtkm/cont/DynamicCellSet.h>
#include <vtkm/cont/DataSet.h>
#include <vtkm/cont/Invoker.h>
#include "Types.h"

//namespace galois
//{
using vtkm::cont::ArrayHandle;
using vtkm::cont::CoordinateSystem;
using vtkm::cont::CellSetExplicit;
using vtkm::cont::CellSetStructured;
using vtkm::cont::Invoker;
using vtkm::cont::DynamicCellSet;
using vtkm::cont::DataSet;
using vtkm::cont::UnknownCellSet;
//}
