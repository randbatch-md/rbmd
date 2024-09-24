#pragma once
#include <vtkm/TypeTraits.h>
#include <vtkm/Types.h>
#include <vtkm/VecTraits.h>

#include <limits> // must be found with or without CUDA.
#ifndef VTKM_CUDA
#include <cmath>
#include <cstring>
#include <limits.h>
#include <math.h>
#include <stdlib.h>
#endif // !VTKM_CUDA

#if !defined(VTKM_CUDA_DEVICE_PASS)
#define VTKM_USE_STL
#include <algorithm>
#endif

#ifdef VTKM_MSVC
#include <intrin.h>                // For bitwise intrinsics (__popcnt, etc)
#include <vtkm/internal/Windows.h> // for types used by MSVC intrinsics.
#ifndef VTKM_CUDA
#include <math.h>
#endif // VTKM_CUDA
#endif // VTKM_MSVC

#define VTKM_CUDA_MATH_FUNCTION_32(func) func##f
#define VTKM_CUDA_MATH_FUNCTION_64(func) func

namespace vtkm
{

//-----------------------------------------------------------------------------
namespace detail
{
template<typename T>
struct FloatingPointReturnType1
{
  using ctype = typename vtkm::VecTraits<T>::ComponentType;
  using representable_as_float_type =
    std::integral_constant<bool,
                           ((sizeof(ctype) < sizeof(float)) ||
                            std::is_same<ctype, vtkm::Float32>::value)>;
  using Type = typename std::
    conditional<representable_as_float_type::value, vtkm::Float32, vtkm::Float64>::type;
};
} // namespace detail
inline VTKM_EXEC_CONT vtkm::Float32 ERF(vtkm::Float32 x)
{
#ifdef VTKM_CUDA
  return VTKM_CUDA_MATH_FUNCTION_32(erf)(x);
#else
  return std::erf(x);
#endif
}

inline VTKM_EXEC_CONT vtkm::Float64 ERF(vtkm::Float64 x)
{
#ifdef VTKM_CUDA
  return VTKM_CUDA_MATH_FUNCTION_64(erf)(x);
#else
  return std::erf(x);
#endif
}
template<typename T>
static inline VTKM_EXEC_CONT typename detail::FloatingPointReturnType1<T>::Type ERF(const T& x)
{
  using RT = typename detail::FloatingPointReturnType1<T>::Type;
  return vtkm::ERF(static_cast<RT>(x));
}
template<typename T, vtkm::IdComponent N>
static inline VTKM_EXEC_CONT vtkm::Vec<typename detail::FloatingPointReturnType1<T>::Type, N> ERF(
  const vtkm::Vec<T, N>& x)
{
  vtkm::Vec<typename detail::FloatingPointReturnType1<T>::Type, N> result;
  for (vtkm::IdComponent index = 0; index < N; index++)
  {
    result[index] = vtkm::ERF(x[index]);
  }
  return result;
}
template<typename T>
static inline VTKM_EXEC_CONT vtkm::Vec<typename detail::FloatingPointReturnType1<T>::Type, 4> ERF(
  const vtkm::Vec<T, 4>& x)
{
  return vtkm::Vec<typename detail::FloatingPointReturnType1<T>::Type, 4>(
    vtkm::ERF(x[0]), vtkm::ERF(x[1]), vtkm::ERF(x[2]), vtkm::ERF(x[3]));
}
template<typename T>
static inline VTKM_EXEC_CONT vtkm::Vec<typename detail::FloatingPointReturnType1<T>::Type, 3> ERF(
  const vtkm::Vec<T, 3>& x)
{
  return vtkm::Vec<typename detail::FloatingPointReturnType1<T>::Type, 3>(
    vtkm::ERF(x[0]), vtkm::ERF(x[1]), vtkm::ERF(x[2]));
}
template<typename T>
static inline VTKM_EXEC_CONT vtkm::Vec<typename detail::FloatingPointReturnType1<T>::Type, 2> ERF(
  const vtkm::Vec<T, 2>& x)
{
  return vtkm::Vec<typename detail::FloatingPointReturnType1<T>::Type, 2>(vtkm::ERF(x[0]),
                                                                          vtkm::ERF(x[1]));
}
}
