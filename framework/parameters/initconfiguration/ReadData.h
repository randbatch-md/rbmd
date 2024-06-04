#pragma once
#include "InitConfiguration.h"

class ReadData : public InitConfiguration
{
public:
  ReadData(const Configuration& cfg);
  ~ReadData() = default;
  void Execute() override;
  void SetPara() override;

  void InitField();

protected:
};