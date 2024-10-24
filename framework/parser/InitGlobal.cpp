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
