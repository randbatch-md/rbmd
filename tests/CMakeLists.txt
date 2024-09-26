##==================================================================================
##  RBMD 1.0 is developed for random batch molecular dynamics calculation.
##  
##  Copyright (C) 2024 SHANGHAI JIAOTONG UNIVERSITY CHONGQING RESEARCH INSTITUTE
##  
##  This program is free software: you can redistribute it and/or modify
##  it under the terms of the GNU General Public License as published by
##  the Free Software Foundation, either version 3 of the License, or
##  (at your option) any later version.
##  
##  This program is distributed in the hope that it will be useful,
##  but WITHOUT ANY WARRANTY; without even the implied warranty of
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##  GNU General Public License for more details.
##  
##  You should have received a copy of the GNU General Public License
##  along with this program.  If not, see <https://www.gnu.org/licenses/>.
##  
##  Contact Email: [your-email@example.com]
##==================================================================================

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