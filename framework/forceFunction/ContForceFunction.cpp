﻿//==================================================================================
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

#include "forceFunction/ContForceFunction.h"
#include <vtkm/cont/ArrayHandle.h>

VTKM_CONT ExecForceFunction ContForceFunction::PrepareForExecution(vtkm::cont::DeviceAdapterId device,
                                       vtkm::cont::Token& token) const
{
  return ExecForceFunction(this->_cut_Off,
                           this->_alpha,
                           this->_volume,
                           this->_Vlength,
                           this->_box,
                           this->Kmax,
                           this->RBEP,
                           this->_rhomax,
                           this->_nrho,
                           this->_drho,
                           this->_nr,
                           this->_dr);
}
