#pragma once
#include "HyperParameters.h"

class ForceField :public HyperParameters
{
public:
  ForceField(const Configuration& cfg);
  ~ForceField(){};
  void Execute() override;

protected:

};