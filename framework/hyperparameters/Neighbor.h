#pragma once
#include "HyperParameters.h"

class Neighbor : public HyperParameters
{
public:
  Neighbor(const Configuration& cfg);
  ~Neighbor(){};
  void Execute() override;
  void SetParameters() override;

protected:
};