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