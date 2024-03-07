#pragma once
#include "InitCondition.h"
#include "MeshFreeSystem.h"

class MeshFreeFileInitCondition : public InitCondition
{
public:
  MeshFreeFileInitCondition(const Configuration& cfg);
  virtual ~MeshFreeFileInitCondition() = default;

  void Execute() override {};
  void UpdateField() override {};

protected:
  std::ifstream _file;
};
