#pragma once

#include "Action.h"

class AddOutputAction : public Action
{
public:
  AddOutputAction(Application& app);
  ~AddOutputAction() = default;

  void Execute() override;
};
