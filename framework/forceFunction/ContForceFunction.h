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
#include "forceFunction/ExecForceFunction.h"
#include "Types.h"


class ContForceFunction : vtkm::cont::ExecutionObjectBase
{
public:

  void SetParameters(const Real& cut_off,
                     const Real& alpha,
                     const Real& volume,
                     const Real& vlength,
                     const Vec3f& box,
                     const IdComponent& kmax)
  {
    this->_cut_Off = cut_off;
    this->_alpha = alpha;
    this->_volume = volume;
    this->_Vlength = vlength;
    this->_box = box;
    this->Kmax = kmax;
  }

  void SetEAMParameters(const Real& rhomax,
                        const Id& nrho,
                        const Real& drho,
                        const Id& nr,
                        const Real& dr)
  {
    this->_rhomax = rhomax;
    this->_nrho = nrho;
    this->_drho = drho;
    this->_nr = nr;
    this->_dr = dr;
  }

  VTKM_CONT
  ExecForceFunction PrepareForExecution(vtkm::cont::DeviceAdapterId device,
                                       vtkm::cont::Token& token) const;

  private:
  Real _cut_Off;
  Real _alpha;
  Real _volume;
  Real _Vlength;
  Vec3f _box;
  IdComponent Kmax;
  IdComponent RBEP;

  Real _rhomax;
  Id _nrho;
  Id _nr;
  Real _drho;
  Real _dr;
};