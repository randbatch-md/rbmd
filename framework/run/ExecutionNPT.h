#pragma once
#include "ExecutionMD.h"

class ExecutionNPT : public ExecutionMD
{
public:
  ExecutionNPT(const Configuration& cfg);

  ~ExecutionNPT() = default;

  void Init() override;
  void Execute() override;
  void SetParameters() override;

protected:
};