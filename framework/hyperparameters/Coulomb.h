#pragma once
#include "HyperParameters.h"

class Coulomb : public HyperParameters
{
public:
  Coulomb(const Configuration& cfg);
  ~Coulomb(){};
  void Execute() override;
  void SetParameters() override;

protected:
};