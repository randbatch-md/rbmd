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

#pragma once
#include <vtkm/cont/ExecutionObjectBase.h>
#include "staticTable/ExecStaticTable.h"
#include "Types.h"


class ContStaticTable : vtkm::cont::ExecutionObjectBase
{
public:
  void SetTableIndex1(const vtkm::cont::ArrayHandle<IdComponent>& table_index1)
  {
    this->_table_index_1 = table_index1;
  }

  void SetTableRij1(const vtkm::cont::ArrayHandle<Real>& table_rij1)
  {
    this->_table_rij_1 = table_rij1;
  }

  void SetTabledRij1(const vtkm::cont::ArrayHandle<Real>& table_drij1)
  {
    this->_table_drij_1 = table_drij1;
  }

  void SetTableFunctionRij1(const vtkm::cont::ArrayHandle<Real>& table_function_rij1)
  {
    this->_table_function_rij_1 = table_function_rij1;
  }

  void SetTabledFunctionRij1(const vtkm::cont::ArrayHandle<Real>& table_dfunction_rij1)
  {
    this->_table_dfunction_rij_1 = table_dfunction_rij1;
  }

  void SetTableIndex2(const vtkm::cont::ArrayHandle<IdComponent>& table_index2)
  {
    this->_table_index_2 = table_index2;
  }

  void SetTableRij2(const vtkm::cont::ArrayHandle<Real>& table_rij2)
  {
    this->_table_rij_2 = table_rij2;
  }

  void SetTabledRij2(const vtkm::cont::ArrayHandle<Real>& table_drij2)
  {
    this->_table_drij_2 = table_drij2;
  }

  void SetTableFunctionRij2(const vtkm::cont::ArrayHandle<Real>& table_function_rij2)
  {
    this->_table_function_rij_2 = table_function_rij2;
  }

  void SetTabledFunctionRij2(const vtkm::cont::ArrayHandle<Real>& table_dfunction_rij2)
  {
    this->_table_dfunction_rij_2 = table_dfunction_rij2;
  }

  VTKM_CONT
  ExecStaticTable PrepareForExecution(vtkm::cont::DeviceAdapterId device,
                                       vtkm::cont::Token& token) const;

  private:
  vtkm::cont::ArrayHandle<IdComponent> _table_index_1;
  vtkm::cont::ArrayHandle<Real> _table_rij_1;
  vtkm::cont::ArrayHandle<Real> _table_drij_1;
  vtkm::cont::ArrayHandle<Real> _table_function_rij_1;
  vtkm::cont::ArrayHandle<Real> _table_dfunction_rij_1;
  vtkm::cont::ArrayHandle<IdComponent> _table_index_2;
  vtkm::cont::ArrayHandle<Real> _table_rij_2;
  vtkm::cont::ArrayHandle<Real> _table_drij_2;
  vtkm::cont::ArrayHandle<Real> _table_function_rij_2;
  vtkm::cont::ArrayHandle<Real> _table_dfunction_rij_2;


};