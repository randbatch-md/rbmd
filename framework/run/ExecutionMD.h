#pragma once
#include "Execution.h"

class ExecutionMD : public Execution
{
public:
  ExecutionMD(const Configuration& cfg);

  ~ExecutionMD() = default;

  void Init() override;
  void Execute() override;
  void SetParameters() override;

protected:
};