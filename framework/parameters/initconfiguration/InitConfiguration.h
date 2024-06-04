#pragma once
#include "Parameters.h"

class InitConfiguration : public Parameters
{
public:
  InitConfiguration(const Configuration& cfg)
    : Parameters(cfg){};

  ~InitConfiguration() = default;

  void Execute() override {};
  void SetPara() override {};

protected:

  //virtual void UpdateField() = 0;
};