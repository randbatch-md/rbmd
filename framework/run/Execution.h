#pragma once
#include "DataObject.h"
#include "Object.h"
#include "Para.h"
class Execution : public Object
{
public:
  Execution(const Configuration& cfg)
    : Object(cfg)
    , _app(*Get<Application*>("app"))
    , _para(*(_app.GetParameter())){};

  ~Execution() = default;

  virtual void Execute() = 0;
  virtual void SetParameters()
  {
    _para.SetParameter(PARA_ENSEMBLE, Get<std::string>("ensemble"));
    _para.SetParameter(PARA_TEMP_CTRL_TYPE, Get<std::string>("temp_ctrl_type"));
    _para.SetParameter(PARA_PRESS_CTRL_TYPE, Get<std::string>("press_ctrl_type"));
    _para.SetParameter(PARA_TIMESTEP, Get<Real>("timestep"));
    _para.SetParameter(PARA_NUM_STEPS, Get<Real>("num_steps"));
    _para.SetParameter(PARA_TEMPERATURE, GetVectorValue<Real>("temperature"));
    _para.SetParameter(PARA_PRESSURE, GetVectorValue<Real>("pressure"));
  };

protected:
protected:
  Application& _app;
  Para& _para;
};