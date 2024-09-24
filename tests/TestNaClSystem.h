#pragma once
#include <gtest/gtest.h>
class System;
class TestNaClSystem : public ::testing::Test
{
public:
  TestNaClSystem();

protected:
  void NaClSystemFormater(const std::string& type = "serial",
                          const int force_type=1,
                         const double dt = 0.2,
                         const int num_steps = 10);
  void ComputeNaClSystem();

protected:
  float _kbt;
  float _ave_kbt;
  float _ave_pe;
  float _ave_ke;
};