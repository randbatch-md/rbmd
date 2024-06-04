#pragma once
#include "DataObject.h"
#include "Object.h"
#include "Para.h"
class Parameters : public Object
{
public:
  Parameters(const Configuration& cfg)
    : Object(cfg)
    , _app(*Get<Application*>("app"))
    , _para(*(_app.GetParameter())){};

  ~Parameters() = default;

  virtual void Execute(){};
  virtual void SetPara(){};

protected:
  //virtual void UpdateField() = 0;

protected:
  Application& _app;
  Para& _para;
};