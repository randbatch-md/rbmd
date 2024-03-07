#pragma once
#include <gtest/gtest.h>
class System;
class TestSalineSolutionSystem : public ::testing::Test
{
public:
  TestSalineSolutionSystem();

protected:
  void SalineSolutionSystemFormater(const std::string& type = "serial",
                                    const int dims = 10,
                                    const int force_type = 1,
                                    const double dt = 0.002,
                                    const int num_steps = 10);
  void ComputeSalineSolutionSystem();

protected:
  float _kbt;
};