#pragma once
#include "ExecutionMD.h"

class ExecutionNVT : public ExecutionMD
{
public:
  ExecutionNVT(const Configuration& cfg);

  ~ExecutionNVT() = default;

  void Init() override;
  void Execute() override;
  void SetParameters();

protected:
};