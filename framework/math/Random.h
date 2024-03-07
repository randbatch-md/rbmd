//#pragma once
//#ifndef vtk_m_Random_h
//#define vtk_m_Random_h
//
//#include <vtkm/Types.h>
//
//#include <random>
//
//namespace vtkm
//{
//
//namespace internal
//{
//
//template<typename DistributionType>
//class Random
//{
//public:
//  VTKM_SUPPRESS_EXEC_WARNINGS
//  VTKM_EXEC_CONT
//  Random(vtkm::Id seed = 0)
//    : Distribution()
//    , Generator(static_cast<typename DistributionType::result_type>(seed))
//  {
//  }
//
//  VTKM_SUPPRESS_EXEC_WARNINGS
//  VTKM_EXEC_CONT
//  typename DistributionType::result_type operator()()
//  {
//    return this->Distribution(this->Generator);
//  }
//
//
//public:
//  VTKM_EXEC_CONT
//  vtkm::FloatDefault RandUniform()
//  {
//    static internal::Random<std::uniform_real_distribution<vtkm::FloatDefault>> randomGen;
//    return randomGen();
//  }
//
//  VTKM_EXEC_CONT
//  vtkm::FloatDefault RandNormal()
//  {
//    static internal::Random<std::normal_distribution<vtkm::FloatDefault>> randomGen;
//    return randomGen();
//  }
//
//  VTKM_EXEC_CONT
//  vtkm::FloatDefault RandFloat(Real lower, Real upper)
//  {
//    static internal::Random<std::uniform_real_distribution<vtkm::FloatDefault>> randomGen;
//    std::uniform_real_distribution<Real>::param_type param(lower, upper);
//    randomGen.Distribution.param(param);
//    return randomGen();
//  }
//
//  VTKM_EXEC_CONT
//  vtkm::Int32 RandInt32(vtkm::Int32 lower, vtkm::Int32 upper)
//  {
//    static internal::Random<std::uniform_int_distribution<vtkm::Int32>> randomGen;
//    std::uniform_int_distribution<vtkm::Int32>::param_type param(lower, upper);
//    randomGen.Distribution.param(param);
//    return randomGen();
//  }
//
//  VTKM_EXEC_CONT
//  vtkm::UInt32 RandUInt32(vtkm::UInt32 lower, vtkm::UInt32 upper)
//  {
//    static internal::Random<std::uniform_int_distribution<vtkm::UInt32>> randomGen;
//    std::uniform_int_distribution<vtkm::UInt32>::param_type param(lower, upper);
//    randomGen.Distribution.param(param);
//    return randomGen();
//  }
//
//  VTKM_EXEC_CONT
//  vtkm::Int64 RandInt64(vtkm::Int64 lower, vtkm::Int64 upper)
//  {
//    static internal::Random<std::uniform_int_distribution<vtkm::Int64>> randomGen;
//    std::uniform_int_distribution<vtkm::Int64>::param_type param(lower, upper);
//    randomGen.Distribution.param(param);
//    return randomGen();
//  }
//
//  VTKM_EXEC_CONT
//  vtkm::UInt64 RandUInt64(vtkm::UInt64 lower, vtkm::UInt64 upper)
//  {
//    static internal::Random<std::uniform_int_distribution<vtkm::UInt64>> randomGen;
//    std::uniform_int_distribution<vtkm::UInt64>::param_type param(lower, upper);
//    randomGen.Distribution.param(param);
//    return randomGen();
//  }
//
//private:
//  DistributionType Distribution;
//  std::mt19937 Generator;
//};
//
//} // namespace internal
//
//
//
//} // namespace vtkm
//
//#endif // vtk_m_Random_h

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