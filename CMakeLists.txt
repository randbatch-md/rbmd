
cmake_minimum_required(VERSION 3.12...3.15 FATAL_ERROR)

set(CMAKE_CXX_STANDARD 17)

set(project_name rbmd)
project(${project_name} CXX)

set(_this ${project_name})

# 第三方库
add_subdirectory(contrib) 
# 框架核心
add_subdirectory(framework)

# 单元测试
enable_testing()
add_subdirectory(tests)
#option(ENABLE_TESTS "Enable unit tests" OFF)
#if(ENABLE_TESTS)
#  add_subdirectory(tests/googletest) 
#  enable_testing()
#endif() 


# 解算器
#add_subdirectory(module)