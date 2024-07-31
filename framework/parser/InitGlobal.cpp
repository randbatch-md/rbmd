#include "InitGlobal.h"
//#include <mpi/mpi.h> // 添加mpi头文件 （这里需要再cmake中配置 mpi环境）
#ifdef VTKM_HIP
#include <Kokkos_Core.hpp>
#endif // VTKM_HIP

InitGlobal::InitGlobal(std::string& device,int argc, char** argv)
{
  if (device == "duc")
  {
    KokkosInit(argc, argv);
  }
  else if (device == "mpi")
  {
    MPIInit(argc, argv);
  }
}

void InitGlobal::KokkosInit(int argc, char** argv)
{
#ifdef VTKM_HIP
  Kokkos::initialize(argc, argv);
#endif // VTKM_HIP
}

void InitGlobal::MPIInit(int argc, char** argv)
{
  //MPI_Init(&argc, &argv);
}
