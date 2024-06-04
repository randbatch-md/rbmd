#pragma once
#include "InitCondition.h"
#include "MeshFreeSystem.h"

class MeshFreeFileInitCondition : public InitCondition
{
public:
  MeshFreeFileInitCondition(const Configuration& cfg);
  virtual ~MeshFreeFileInitCondition() = default;
  //void InitParameter() override;
  void Execute() override {};
  void UpdateField() override {};
  void InitField() override;
  void SetParameters() override;

protected:
  std::ifstream _file;
  std::string _file_name;
};
