#pragma once

#include "Action.h"

class SetupDeviceAction : public Action
{
public:
  SetupDeviceAction(Application& app);
  ~SetupDeviceAction() = default;

  void Execute() override;
};
