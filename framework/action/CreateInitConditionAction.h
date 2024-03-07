#pragma once

#include "Action.h"

class CreateInitConditionAction : public Action
{
public:
  CreateInitConditionAction(Application& app);
  ~CreateInitConditionAction() = default;

  void Execute() override;
};