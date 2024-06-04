#pragma once
#include "HyperParameters.h"

class ForceField :public HyperParameters
{
public:
  ForceField(const Configuration& cfg);
  ~ForceField(){};
  void Execute() override;
  void SetParameters() override;

protected:

};