#include "InitGlobal.h"

#ifdef VTKM_HIP
#include <Kokkos_Core.hpp>
#endif // VTKM_HIP

InitGlobal::InitGlobal(int argc, char** argv)
{
  KokkosInit(argc, argv);
}

void InitGlobal::KokkosInit(int argc, char** argv)
{
#ifdef VTKM_HIP
  Kokkos::initialize(argc, argv);
#endif // VTKM_HIP
}
