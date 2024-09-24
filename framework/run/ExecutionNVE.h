#pragma once
#include "ExecutionMD.h"

class ExecutionNVE : public ExecutionMD
{
public:
  ExecutionNVE(const Configuration& cfg);

  ~ExecutionNVE() = default;

  void Init() override;
  void Execute() override;
  void SetParameters();

protected:
};