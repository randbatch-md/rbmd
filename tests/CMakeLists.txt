set(_this test_semd)



set(_sources
   TestMadelung.cpp
   TestLennardJonesSystem.cpp
   TestSalineSolutionSystem.cpp
   TestH2OSystem.cpp
   TestNaClSystem.cpp
)

find_package(VTKm REQUIRED QUIET)

add_executable(${_this} 
  ${_sources}

)

link_libraries(framework)
target_link_libraries(${_this}
PUBLIC libsemd gtest_main
)

add_test(
    NAME ${_this}
    COMMAND ${_this}
)