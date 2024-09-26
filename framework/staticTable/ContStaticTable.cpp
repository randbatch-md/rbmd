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

#include "staticTable/ContStaticTable.h"
#include <vtkm/cont/ArrayHandle.h>

VTKM_CONT ExecStaticTable ContStaticTable::PrepareForExecution(vtkm::cont::DeviceAdapterId device,
                                       vtkm::cont::Token& token) const
{
  return ExecStaticTable(this->_table_index_1.PrepareForInput(device, token),
                           this->_table_rij_1.PrepareForInput(device, token),
                           this->_table_drij_1.PrepareForInput(device, token),
                           this->_table_function_rij_1.PrepareForInput(device, token),
                           this->_table_dfunction_rij_1.PrepareForInput(device, token),
                           this->_table_index_2.PrepareForInput(device, token),
                           this->_table_rij_2.PrepareForInput(device, token),
                           this->_table_drij_2.PrepareForInput(device, token),
                           this->_table_function_rij_2.PrepareForInput(device, token),
                           this->_table_dfunction_rij_2.PrepareForInput(device, token)
                           );
}