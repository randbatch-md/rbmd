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

set(CUDAToolkit_DIR /usr/local/cuda-11.8)
find_package(CUDAToolkit REQUIRED)
set(CMAKE_CUDA_ARCHITECTURES "86")
#set(Torch_DIR "/data/zkh_mace/libtorch/share/cmake/Torch")
set(Torch_DIR "/data/zkh_mace/libtorch/share/cmake/Torch")
find_package(Torch REQUIRED)
set(LIBTORCH_PATH "/data/zkh_mace/libtorch")
include_directories(${LIBTORCH_PATH}/include)
include_directories(${LIBTORCH_PATH}/include/torch/csrc/api/include)
include_directories(${LIBTORCH_PATH}/include/torch/csrc/autograd)
include_directories(${LIBTORCH_PATH}/include/torch/csrc/cuda)

#find_package(CUDAToolkit REQUIRED)
#set(CMAKE_CUDA_ARCHITECTURES "86")
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

#加载 VTK-m 
target_link_libraries(libsemd PUBLIC vtkm_cont vtkm_worklet vtkm_io)
target_link_libraries(libsemd PRIVATE ${VTK_LIBRARIES})
target_link_libraries(libsemd PUBLIC fmt::fmt)
vtkm_add_target_information(libsemd DROP_UNUSED_SYMBOLS MODIFY_CUDA_FLAGS DEVICE_SOURCES ${device_src})


#定义target
target_include_directories(libsemd PRIVATE
    ${LIBTORCH_PATH}/include
    ${LIBTORCH_PATH}/include/torch/csrc/api/include
    ${LIBTORCH_PATH}/include/torch/csrc/autograd
    ${LIBTORCH_PATH}/include/torch/csrc/cuda
)
# 设置链接目录
target_link_directories(libsemd PRIVATE ${LIBTORCH_PATH}/lib)

# 链接 libtorch 库
#message(STATUS "Torch Libraries: ${TORCH_LIBRARIES}")
target_link_libraries(libsemd PRIVATE
    ${TORCH_LIBRARIES}#torch
    #torch_library
    #/data/zkh_mace/libtorch/lib/libc10.so
    #/data/zkh_mace/libtorch/lib/libkineto.a
    #/usr/lib/x86_64-linux-gnu/libcuda.so
    #/usr/local/cuda-11.8/lib64/libnvrtc.so
    #/usr/local/cuda-11.8/lib64/libnvToolsExt.so
    #/usr/local/cuda-11.8/lib64/libcudart.so
    #/data/zkh_mace/libtorch/lib/libc10_cuda.so
    #torch torch_cuda
    #torch_cpu
    #c10
    #c10_cuda
    #${LIBTORCH_PATH}/lib/libtorch.so
    #${LIBTORCH_PATH}/lib/libtorch_cuda.so
    #${LIBTORCH_PATH}/lib/libtorch_cpu.so
    #${LIBTORCH_PATH}/lib/libc10.so
    #${LIBTORCH_PATH}/lib/libc10_cuda.so
    #${LIBTORCH_PATH}/lib/libc10d_cuda_test.so
    #${LIBTORCH_PATH}/lib/libbackend_with_compiler.so
    #${LIBTORCH_PATH}/lib/libcaffe2_nvrtc.so
    #${LIBTORCH_PATH}/lib/libtorch_cuda_linalg.so
    #${LIBTORCH_PATH}/lib/libtorch_global_deps.so
    #${LIBTORCH_PATH}/lib/libtorchbind_test.so
    #${LIBTORCH_PATH}/lib/libtorch_python.so
    #${LIBTORCH_PATH}/lib/libpytorch_jni.so
    #${LIBTORCH_PATH}/lib/libshm.so
    #${LIBTORCH_PATH}/lib/libnnapi_backend.so
    #${LIBTORCH_PATH}/lib/libjitbackend_test.so
    #${LIBTORCH_PATH}/lib/libfbjni.so

)
#target_link_libraries(libsemd PRIVATE   torch torch_cuda torch_cpu c10 c10_cuda )
#target_link_libraries(libsemd PRIVATE    torch    torch_cpu    torch_cuda  )
target_link_libraries(libsemd PRIVATE CUDA::cudart)
target_link_libraries(libsemd PRIVATE ${CUDA_LIBRARIES})
target_include_directories(libsemd PRIVATE ${CUDA_INCLUDE_DIRS})
add_executable(${target} main.cpp)
#file(COPY ${CMAKE_SOURCE_DIR}/framework/res/Gnear_table_geq2_2.txt DESTINATION ${CMAKE_BINARY_DIR}/framework/res)
#file(COPY ${CMAKE_SOURCE_DIR}/framework/res/Gnear_table_leq2_2.txt DESTINATION ${CMAKE_BINARY_DIR}/framework/res)
#file(COPY "../res" DESTINATION ${CMAKE_BINARY_DIR}/framework)
#file(COPY "../res" DESTINATION ${CMAKE_BINARY_DIR}/framework)



target_link_libraries(${target} PUBLIC libsemd )
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
