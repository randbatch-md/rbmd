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

# 依赖vtkm
cmake_minimum_required(VERSION 3.8)
project(framework)

#安装cxxopts
include(./cmake/cxxopts.cmake)
setup_cxxopts()

#安装json
include(./cmake/jsoncpp.cmake)
setup_jsoncpp()

if(NOT ${CMAKE_SYSTEM_NAME} STREQUAL "Windows")
set(VTKm_DIR "${CMAKE_CURRENT_LIST_DIR}/../../vtkm_out/lib/cmake/vtkm-1.9" CACHE PATH "VTK directory override" FORCE)
endif()

find_package(VTKm REQUIRED QUIET)

set(target ${project_name})

set(host_src
parser/Configuration.cpp
application/Application.cpp
application/MDApplication.cpp
executioner/Executioner.cpp
output/Output.cpp
output/FileOutput.cpp
output/ConsoleOutput.cpp
output/Progress.cpp
model/FormatTable.cpp
InitCondition/MeshFreeCondition.cpp
InitCondition/ModelFileInitCondition.cpp
InitCondition/MeshFreeFileInitCondition.cpp
hyperparameters/Coulomb.cpp
hyperparameters/Extend.cpp
hyperparameters/ForceField.cpp
hyperparameters/Neighbor.cpp
)

set(device_src
output/ThermoOutput.cpp
output/TempOutput.cpp
output/RDFOutput.cpp
output/MSDOutput.cpp
output/VACFOutput.cpp
output/worklet/OutPutWorklet.cpp
output/TrajectoryOutput.cpp
InitCondition/LJInitCondition.cpp
InitCondition/MadelungInitCondition.cpp
locator/ContPointLocator.cxx
forceFunction/ContForceFunction.cpp
topology/ContTopology.cpp
staticTable/ContStaticTable.cpp
run/ExecutionNPT.cpp
run/ExecutionNVT.cpp
run/ExecutionNVE.cpp
run/ExecutionMD.cpp
run/worklet/RunWorklet.cpp
parser/InitGlobal.cpp
)

# 头文件引入
add_library(_compiler_flags INTERFACE)
target_include_directories(_compiler_flags INTERFACE
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/application
  ${CMAKE_CURRENT_SOURCE_DIR}/executioner
  ${CMAKE_CURRENT_SOURCE_DIR}/mesh
  ${CMAKE_CURRENT_SOURCE_DIR}/meshmodifier
  ${CMAKE_CURRENT_SOURCE_DIR}/output
  ${CMAKE_CURRENT_SOURCE_DIR}/parser
  ${CMAKE_CURRENT_SOURCE_DIR}/userobject
  ${CMAKE_CURRENT_SOURCE_DIR}/model
  ${CMAKE_CURRENT_SOURCE_DIR}/InitCondition
  ${CMAKE_CURRENT_SOURCE_DIR}/locator
  ${CMAKE_CURRENT_SOURCE_DIR}/forceFunction
  ${CMAKE_CURRENT_SOURCE_DIR}/topology
  ${CMAKE_CURRENT_SOURCE_DIR}/hyperparameters
  ${CMAKE_CURRENT_SOURCE_DIR}/run
  ${JSON_INSTALL_PATH}/include
  ${CXXOPTS_INSTALL_PATH}/include
  )
  # 定义lib
add_library(librbmd STATIC ${host_src} ${device_src})
target_link_libraries(librbmd PUBLIC _compiler_flags)

#链接json库
target_link_directories(librbmd PUBLIC ${JSON_INSTALL_PATH}/lib)
target_link_libraries(librbmd PRIVATE jsoncpp)

#加载 VTK-m 
target_link_libraries(librbmd PUBLIC vtkm_cont vtkm_worklet vtkm_io)
target_link_libraries(librbmd PRIVATE ${VTK_LIBRARIES})
target_link_libraries(librbmd PUBLIC fmt::fmt)
vtkm_add_target_information(librbmd DROP_UNUSED_SYMBOLS MODIFY_CUDA_FLAGS DEVICE_SOURCES ${device_src})

  # 定义target
add_executable(${target} main.cpp)



target_link_libraries(${target} PUBLIC librbmd )
vtkm_add_target_information(${target} DROP_UNUSED_SYMBOLS MODIFY_CUDA_FLAGS DEVICE_SOURCES main.cpp)

install(TARGETS ${target} 
  DESTINATION "${CMAKE_INSTALL_PREFIX}/bin"
  )
  
  # 在可执行文件的目标目录下创建 lib 文件夹
set_target_properties(librbmd PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/lib"
)

# 确保 lib 文件夹在构建时存在
file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/lib")

   # 设置 VTKM 许可证文件路径
 set(VTKM_LICENSE_PATH "${CMAKE_CURRENT_LIST_DIR}/../licenses/LICENSE.TXT")
 
 # 在生成 librbmd 静态库后复制许可证文件
 add_custom_command(TARGET librbmd POST_BUILD
     COMMAND ${CMAKE_COMMAND} -E copy ${VTKM_LICENSE_PATH} $<TARGET_FILE_DIR:librbmd>
 )
 
 # 在生成可执行文件后复制许可证文件
 add_custom_command(TARGET ${target} POST_BUILD
     COMMAND ${CMAKE_COMMAND} -E copy ${VTKM_LICENSE_PATH} $<TARGET_FILE_DIR:${target}>
 )
