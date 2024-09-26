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

#pragma once
#include <gtest/gtest.h>
class System;
class TestLennardJonesSystem : public ::testing::Test
{
public:
  TestLennardJonesSystem();

protected:
  void LennardJonesSystemFormater(const std::string& type = "serial",
                                  const int dims = 10,
                                  const double dt = 0.002,
                                  const int num_steps = 10);
  void ComputeLennardJonesSystem();

protected:
  float _kbt;
};