#include "ExecutionMD.h"
#include "FieldName.h"

ExecutionMD::ExecutionMD(const Configuration& cfg)
  : Execution(cfg)
{
  Execution::SetParameters();
}

void ExecutionMD::Init() {}

void ExecutionMD::Execute()
{
  //_timer.Start();
}

void ExecutionMD::SetParameters() 
{
  _para.SetParameter(PARA_ENSEMBLE, Get<std::string>("ensemble"));
  _para.SetParameter(PARA_TEMP_CTRL_TYPE, Get<std::string>("temp_ctrl_type"));
  _para.SetParameter(PARA_PRESS_CTRL_TYPE, Get<std::string>("press_ctrl_type"));
  _para.SetParameter(PARA_TIMESTEP, Get<Real>("timestep"));
  _para.SetParameter(PARA_NUM_STEPS, Get<Real>("num_steps"));
  _para.SetParameter(PARA_TEMPERATURE, GetVectorValue<Real>("temperature"));
  _para.SetParameter(PARA_PRESSURE, GetVectorValue<Real>("pressure"));
}