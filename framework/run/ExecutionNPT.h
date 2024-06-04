#pragma once
#include "Execution.h"

class ExecutionNPT : public Execution
{
public:
  ExecutionNPT(const Configuration& cfg);

  ~ExecutionNPT() = default;

  void Execute() override;
  void SetParameters();

protected:
};