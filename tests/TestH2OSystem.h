#pragma once
#include <gtest/gtest.h>
class System;
class TestH2OSystem : public ::testing::Test
{
public:
  TestH2OSystem();

protected:
  void H2OSystemFormater(const std::string& type = "serial",
                         const int force_type = 1,
                         const double dt = 0.2,
                         const int num_steps = 10);
  void ComputeH2OSystem();

protected:
  float _kbt;
  float _ave_kbt;
  float _ave_pe;
  float _ave_ke;
};