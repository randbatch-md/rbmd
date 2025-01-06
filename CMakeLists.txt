##==================================================================================
##  RBMD 2.2.0 is developed for random batch molecular dynamics calculation.
##
##  Copyright(C) 2024 SHANGHAI JIAOTONG UNIVERSITY CHONGQING RESEARCH INSTITUTE
##
##  This program is free software : you can redistribute it and /or modify
##  it under the terms of the GNU General Public License as published by
##  the Free Software Foundation, either version 3 of the License, or
##  (at your option) any later version.
##
##  This program is distributed in the hope that it will be useful,
##  but WITHOUT ANY WARRANTY; without even the implied warranty of
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
##  GNU General Public License for more details.
##
##  You should have received a copy of the GNU General Public License
##  along with this program.If not, see < https://www.gnu.org/licenses/>.
##
##  The post-processing data produced by VASPKIT may be used in any publications 
##  provided that its use is explicitly acknowledged. A suitable reference for VASPKIT is:
##  [1] Gao W, Zhao T, Guo Y, et al.RBMD: A molecular dynamics package enabling to simulate 
##  10 million all - atom particles in a single graphics processing unit[J].arXiv preprint arXiv : 2407.09315, 2024.
## 
##  Contact Email : [support_wz@sciai.com.cn]
##==================================================================================

cmake_minimum_required(VERSION 3.12...3.15 FATAL_ERROR)

set(CMAKE_CXX_STANDARD 17)

set(project_name rbmd)
project(${project_name} CXX)

set(_this ${project_name})

# 第三方库
add_subdirectory(contrib) 
# 框架核心
add_subdirectory(framework)

# 单元测试 (配置文件发生变化单元测试已注释)
#enable_testing()
#add_subdirectory(tests)
#option(ENABLE_TESTS "Enable unit tests" OFF)
#if(ENABLE_TESTS)
#  add_subdirectory(tests/googletest) 
#  enable_testing()
#endif() 


# 解算器
#add_subdirectory(module)