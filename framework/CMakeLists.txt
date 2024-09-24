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


#set(CUDAToolkit_DIR /usr/local/cuda-11.8)
#set(CUDAToolkit_DIR /public/software/compiler/nvidia/cuda/12.1)
find_package(CUDAToolkit REQUIRED)
set(CMAKE_CUDA_ARCHITECTURES  3.5 5.0 8.0 8.6)
set(CUDNN_INCLUDE_DIR /work/home/zhaoteng_sciai/zkh/cudnn/include)
set(CUDNN_LIBRARY /work/home/zhaoteng_sciai/zkh/cudnn/lib64/libcudnn.so)
set(CUDNN_DIR /work/home/zhaoteng_sciai/zkh/cudnn)

include_directories(${CUDNN_INCLUDE_DIR})
link_directories(${CUDNN_LIBRARY})

#find_package(CUDNN REQUIRED)
#find_package(CUSPARSELT REQUIRED)
if(CUDNN_FOUND)
    add_definitions(-DUSE_CUDNN=1)
endif()

if(CUSPARSELT_FOUND)
    add_definitions(-DUSE_CUSPARSELT=1)
endif()

#set(CMAKE_PREFIX_PATH "/public/home/zkh/zkh_mace/libtorch/share/cmake/Torch")
#set(Torch_DIR "/public/software/life/libtorch-gpu/share/cmake/Torch")
set(TORCH_CUDA_ARCH_LIST "6.1;7.0;8.0;8.6;8.9;9.0")
set(Torch_DIR "/work/home/zhaoteng_sciai/zkh/libtorch/share/cmake/Torch")
find_package(Torch REQUIRED)
set(LIBTORCH_PATH "/work/home/zhaoteng_sciai/zkh/libtorch")
#message(STATUS "Torch CXX Flags: ${TORCH_CXX_FLAGS}")
#set(LIBTORCH_PATH "/public/software/life/libtorch-gpu")
#set(LIBTORCH_PATH "/public/home/zhaoteng_sciai/zkh/libtorch")
include_directories(${LIBTORCH_PATH}/include)
include_directories(${LIBTORCH_PATH}/include/torch/csrc/api/include)
include_directories(${LIBTORCH_PATH}/include/torch/csrc/autograd)
include_directories(${LIBTORCH_PATH}/include/torch/csrc/cuda)
link_directories(${LIBTORCH_PATH}/lib)

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
mace/maceload.cpp
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
  ${CMAKE_CURRENT_SOURCE_DIR}/mace
  ${JSON_INSTALL_PATH}/include
  ${CXXOPTS_INSTALL_PATH}/include
  )
  # 看来要在res 文件下建立cmake文件
  # 定义lib
add_library(libsemd STATIC ${host_src} ${device_src})
target_link_libraries(libsemd PUBLIC _compiler_flags)

#链接json库
#set(CMAKE_PREFIX_PATH ${JSON_INSTALL_PATH}/lib/cmake)
#find_package(jsoncpp REQUIRED)
target_link_directories(libsemd PUBLIC ${JSON_INSTALL_PATH}/lib)
target_link_libraries(libsemd PRIVATE jsoncpp)
#find_package(fmt REQUIRED)
#加载 VTK-m 
target_link_libraries(libsemd PUBLIC vtkm_cont vtkm_worklet vtkm_io)
target_link_libraries(libsemd PRIVATE ${VTK_LIBRARIES})
message(STATUS fmt::fmt)
target_link_libraries(libsemd PUBLIC fmt::fmt)
vtkm_add_target_information(libsemd DROP_UNUSED_SYMBOLS MODIFY_CUDA_FLAGS DEVICE_SOURCES ${device_src})

#add_library(precxxtorch INTERFACE)
#target_compile_definitions(precxxtorch INTERFACE GLIBCXX_USE_CXX11_ABI=0)
#定义target
#target_include_directories(precxxtorch INTERFACE
target_include_directories(libsemd PUBLIC
    ${LIBTORCH_PATH}/include
    ${LIBTORCH_PATH}/include/torch/csrc/api/include
    ${LIBTORCH_PATH}/include/torch/csrc/autograd
    ${LIBTORCH_PATH}/include/torch/csrc/cuda
)
#target_link_directories(libsemd PUBLIC ${JSON_INSTALL_PATH}/lib)
#target_link_libraries(libsemd PRIVATE jsoncpp)

# 设置链接目录
#target_link_directories(precxxtorch INTERFACE ${LIBTORCH_PATH}/lib)
#target_link_directories(libsemd PUBLIC ${LIBTORCH_PATH}/lib)
# 链接 libtorch 库
message(STATUS "Torch Libraries: ${TORCH_LIBRARIES}")
#target_link_libraries(precxxtorch INTERFACE
target_link_libraries(libsemd PRIVATE
    ${TORCH_LIBRARIES}
    #${LIBTORCH_PATH}/lib/libtorch.so
    #torch
    #torch_library
    #/data/zkh_mace/libtorch/lib/libc10.so
    #/data/zkh_mace/libtorch/lib/libkineto.a
    #/usr/lib/x86_64-linux-gnu/libcuda.so
    #/usr/local/cuda-11.8/lib64/libnvrtc.so
    #/usr/local/cuda-11.8/lib64/libnvToolsExt.so
    #/usr/local/cuda-11.8/lib64/libcudart.so
    #/data/zkh_mace/libtorch/lib/libc10_cuda.so
    #torch# torch_cuda
    #torch_cpu
    #c10
    #c10_cuda
)
#torch_add_target_information(libsemd  MODIFY_CUDA_FLAGS DEVICE_SOURCES ${device_src})

message(STATUS "CUDA Libraries: ${CUDA_LIBRARIES}")
message(STATUS "cuda include dirs: ${CUDA_INCLUDE_DIRS}")
target_link_libraries(libsemd PRIVATE CUDA::cudart)
target_link_libraries(libsemd PRIVATE ${CUDA_LIBRARIES} ${CUDNN_LIBRARY})
target_include_directories(libsemd PRIVATE ${CUDA_INCLUDE_DIRS} ${CUDNN_INCLUDE_DIR})
add_executable(${target} main.cpp)
#file(COPY ${CMAKE_SOURCE_DIR}/framework/res/Gnear_table_geq2_2.txt DESTINATION ${CMAKE_BINARY_DIR}/framework/res)
#file(COPY ${CMAKE_SOURCE_DIR}/framework/res/Gnear_table_leq2_2.txt DESTINATION ${CMAKE_BINARY_DIR}/framework/res)
#file(COPY "../res" DESTINATION ${CMAKE_BINARY_DIR}/framework)
#file(COPY "../res" DESTINATION ${CMAKE_BINARY_DIR}/framework)


#target_link_libraries(libsemd PRIVATE precxxtorch)

target_link_libraries(${target} PUBLIC libsemd )
#target_link_libraries(${target} PUBLIC libsemd precxxtorch)
vtkm_add_target_information(${target} DROP_UNUSED_SYMBOLS MODIFY_CUDA_FLAGS DEVICE_SOURCES main.cpp)

install(TARGETS ${target} 
  DESTINATION "${CMAKE_INSTALL_PREFIX}/bin"
  )

#install(DIRECTORY "../res"
#        DESTINATION "${CMAKE_INSTALL_PREFIX}/bin")
        #FILES_MATCHING
        #PATTERN "Gnear_table_geq2_2.txt"
        #PATTERN "Gnear_table_leq2_2.txt")

#install(DIRECTORY "../res/"
#        DESTINATION "${CMAKE_INSTALL_PREFIX}/bin"
#        FILES_MATCHING
#        PATTERN "LJFluid.json"
#        PATTERN "rbmd.data"
#        PATTERN "rbmd.json")
