#pragma once
#include "Parameters.h"

class HyperParameters : public Parameters
{
public:
  HyperParameters(const Configuration& cfg)
    : Parameters(cfg){};

  ~HyperParameters() = default;

   void Execute() override {};
   void SetPara() override {};

protected:
  //virtual void UpdateField() = 0;
};