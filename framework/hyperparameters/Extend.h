#pragma once
#include "HyperParameters.h"

class Extend : public HyperParameters
{
public:
  Extend(const Configuration& cfg);
  ~Extend(){};
  void Execute() override;
  void SetParameters() override;

protected:
};