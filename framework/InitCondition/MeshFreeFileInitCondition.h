#pragma once
#include "InitCondition.h"

class MeshFreeFileInitCondition : public InitCondition
{
public:
  MeshFreeFileInitCondition(const Configuration& cfg);
  virtual ~MeshFreeFileInitCondition() = default;

  void Execute() override {};
  void UpdateField() override {};
  void InitField() override;
  void SetParameters() override;

protected:
  std::ifstream _file;
};
