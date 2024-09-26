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

#include "topology/ContTopology.h"
#include <vtkm/cont/ArrayHandle.h>
#include  "ContTopology.h"

VTKM_CONT ExecTopology ContTopology::PrepareForExecution(vtkm::cont::DeviceAdapterId device,
                                                         vtkm::cont::Token& token) const
{
  return ExecTopology(this->_group_vec_array.PrepareForInput(device, token), 
                      this->_molecular_id.PrepareForInput(device, token),
                      this->_atoms_type.PrepareForInput(device, token),
                      this->_charge.PrepareForInput(device, token),
                      this->_epsilon.PrepareForInput(device, token),
                      this->_sigma.PrepareForInput(device, token) );
                      // this->_group_vec_neighbour.PrepareForInput(device, token) 
                      
}