#pragma once
#include <gtest/gtest.h>
class System;
class TestMadelung : public ::testing::Test
{
public:
  TestMadelung();

protected:
  void MadelungFormater(const std::string& type = "serial", const int dims = 101);
  void ComputeMadelung();

protected:
  float _madelung_result;
};