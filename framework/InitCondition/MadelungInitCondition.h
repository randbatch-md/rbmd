#pragma once

#include "InitCondition.h"

class MadelungInitCondition : public InitCondition
{
public:
  struct SetChargeWorklet;

public:
  MadelungInitCondition(const Configuration& cfg);
  ~MadelungInitCondition() = default;

protected:
  void Execute() override;
  void UpdateField() override;

private:
  std::vector<int> _dims;
  std::vector<Real> _x_range;
  std::vector<Real> _y_range;
  std::vector<Real> _z_range;
};