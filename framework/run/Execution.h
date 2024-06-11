#pragma once
#include "DataObject.h"
#include "Object.h"
#include "Para.h"
#include "Timer.h"
#include "FieldName.h"
class Execution : public Object
{
public:
  Execution(const Configuration& cfg)
    : Object(cfg)
    , _app(*Get<Application*>("_app"))
    , _para(*(_app.GetParameter())){};

  ~Execution() = default;

  virtual void Init(){};
  virtual void Execute() = 0;
  virtual void SetParameters(){};
  auto& GetTimer() { return this->_timer; }

protected:
  Application& _app;
  Para& _para;
  Timer _timer;
};