#pragma once
#include "DataObject.h"
#include "Object.h"
#include "Para.h"
class HyperParameters : public Object
{
public:
  HyperParameters(const Configuration& cfg)
    : Object(cfg)
    , _app(*Get<Application*>("_app"))
    , _para(*(_app.GetParameter())){};

  ~HyperParameters() = default;

  virtual void Execute() = 0;

protected:

protected:
  Application& _app;
  Para& _para;
};