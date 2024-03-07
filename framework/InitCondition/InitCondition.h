#pragma once
#include "Object.h"
#include "DataObject.h"
#include "System.h"

class InitCondition : public Object
{
public:
  InitCondition(const Configuration& cfg)
    : Object(cfg)
    , _app(*Get<Application*>("_app"))
    , _system(*(_app.GetSystem())){};
  ~InitCondition() = default;

  virtual void Execute() = 0;

protected:
  virtual void UpdateField()=0;

protected:
  Application& _app;
  System& _system;
};