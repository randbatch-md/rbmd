//==================================================================================
//  RBMD 2.1.0 is developed for random batch molecular dynamics calculation.
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
#ifndef vtk_m_Random_h
#define vtk_m_Random_h

#include <vtkm/Types.h>
#include <vtkm/cont/testing/MakeTestDataSet.h>
#include <random>
#include "Types.h"

namespace vtkm
{

namespace internal
{

template<typename DistributionType>
class Random
{
public:
  VTKM_SUPPRESS_EXEC_WARNINGS
  VTKM_EXEC_CONT
  Random(vtkm::Id seed = 0)
    : Distribution()
    , Generator(static_cast<typename DistributionType::result_type>(seed))
  {
  }

  VTKM_SUPPRESS_EXEC_WARNINGS
  VTKM_EXEC_CONT
  typename DistributionType::result_type operator()()
  {
    return this->Distribution(this->Generator);
  }


public:
  VTKM_EXEC_CONT
  vtkm::FloatDefault RandUniform()
  {
    static internal::Random<DistributionType> randomGen;
    return randomGen();
  }

  VTKM_EXEC_CONT
  vtkm::FloatDefault RandFloat(Real lower, Real upper)
  {
    static internal::Random<DistributionType> randomGen;
    DistributionType::param_type param(lower, upper);
    randomGen.Distribution.param(param);
    return randomGen();
  }

private:
  DistributionType Distribution;
  std::mt19937 Generator;
};

} // namespace internal
} // namespace vtkm

#endif // vtk_m_Random_h