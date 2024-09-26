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

#include "InitGlobal.h"
//#include <mpi/mpi.h> // 添加mpi头文件 （这里需要再cmake中配置 mpi环境）
#ifdef VTKM_HIP
#include <Kokkos_Core.hpp>
#endif // VTKM_HIP

InitGlobal::InitGlobal(const std::string& device,int argc, char** argv)
{
  if (device == "dcu")
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
