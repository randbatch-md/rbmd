#pragma once
#include "HyperParameters.h"

class Extend : public HyperParameters
{
public:
  Extend(const Configuration& cfg);
  ~Extend(){};
  void Execute() override;
  void SetPara() override;

protected:
};