#pragma once
#include "Parameters.h"

class ExecutionParameters : public Parameters
{
public:
  ExecutionParameters(const Configuration& cfg)
    : Parameters(cfg){};

  ~ExecutionParameters() = default;

  void Execute();
  void SetPara();

protected:
  //virtual void UpdateField() = 0;
};