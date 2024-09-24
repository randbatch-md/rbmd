#pragma once
#include "Object.h"
#include "DataObject.h"
#include "Para.h"

class InitCondition : public Object
{
public:
  InitCondition(const Configuration& cfg)
    : Object(cfg)
    , _app(*Get<Application*>("_app"))
    , _para(*(_app.GetParameter()))
  {
  };

  ~InitCondition() = default;

  virtual void Execute() = 0;
  virtual void InitField() = 0;
  virtual void SetParameters() = 0;

protected:
  virtual void UpdateField()=0;

protected:
  Application& _app;
  Para& _para;
};