#pragma once

#include "Action.h"

class CreateSystemAction : public Action
{
public:
  CreateSystemAction(Application& app);
  ~CreateSystemAction() = default;

  virtual void Execute() override;
};
