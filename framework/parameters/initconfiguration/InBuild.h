#pragma once
#include "InitConfiguration.h"

class InBuild : public InitConfiguration
{
public:
  InBuild(const Configuration& cfg);
  ~InBuild() = default;
  void Execute() override;
  void SetPara() override;

  void InitField();

protected:
};