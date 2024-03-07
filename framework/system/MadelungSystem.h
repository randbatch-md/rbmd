#pragma once
#include "System.h"

class MadelungSystem : public System
{
public:
  struct InitialCondtionWorklet;
  struct ComputeCoeffient;

public:
  MadelungSystem(const Configuration& cfg);
  virtual ~MadelungSystem(){};
  void Init() override;

protected:
  void Evolve() override;
  void InitField() override;

private:
  void TimeIntegration(){};
  void InitialCondition(){};
  void virtual UpdateResidual(){};

private:
  Real _madelung_result;
};