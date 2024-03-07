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