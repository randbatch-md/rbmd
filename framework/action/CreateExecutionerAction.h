#pragma once
#include "Action.h"

class CreateExecutionerAction : public Action
{
public:
  CreateExecutionerAction(Application& app);
  ~CreateExecutionerAction() = default;

  virtual void Execute() override;
};
