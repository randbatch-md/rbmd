#include "Executioner.h"
#include "Application.h"
#include "Logging.h"
#include "Output.h"
#include <iostream>    
#include <iomanip>  

Executioner::Executioner(const Configuration& cfg)
  : Object(cfg)
  , _app(*(Get<Application*>("_app")))
  , _para(*(_app.GetParameter()))
  , _run(*(_app.GetRun()))
  , _timer(_run.GetTimer())
  , _dt(_para.GetParameter<Real>(PARA_TIMESTEP))
  , _start_time(0)
  , _current_time(0)
  , _num_steps(_para.GetParameter<Real>(PARA_NUM_STEPS))
  , _end_time(std::numeric_limits<Real>::max())
{
  _current_step = 0;
  _end_time = _dt * _num_steps;
  _timer.Start();
  _time_count.Start();
}

void Executioner::Init() 
{
  _run.Init();

  for (auto& output : _app.GetOutputWarehouse())
  {
    output->Init();
  }
}

void Executioner::Execute()
{
  PreExecute();

  while (KeepGoing())
  {
    Timer time;
    time.Start();
    PreStep();
    TakeStep();
    PostStep();
  }

  PostExecute();
}

void Executioner::PostExecute()
{
  console::Success("Running successfully.  Total time: ", _time_count.GetElapsedTime(), "s");
  //console::Info("当前步: ", _current_step-1, "总步数: ", _num_steps);
  console::Info("Time step: ", _dt, "fs      ", "Total simulation time: ", _end_time, "fs");
  std::ofstream log_file("rbmd.log",std::ios::app);
  try
  {
    log_file << fmt::format(u8"运行完成.\n运行时间:   {0}\n当前步:   {1}\n总步数:   {2}\n当前时间: "
                            u8"  {3}\n结束时间:   {4}\n",
                            _timer.GetElapsedTime(),
                            _current_step - 1,
                            _num_steps,
                            _current_time - _dt,
                            _end_time);
    log_file.close();
  }
  catch (const std::exception& e)
  {
    log_file.close();
    console::Error(e.what());
  }
}

void Executioner::PreStep() 
{
  ComputeDT();

  for (auto& output : _app.GetOutputWarehouse())
  {
    output->Execute();
  }
}

void Executioner::TakeStep() 
{
  _current_step++;
  _current_time += _dt;
  _run.Execute();
}

void Executioner::PostStep() {}

bool Executioner::KeepGoing()
{
  bool keep_going = false;

  if (_current_step <= _num_steps)
    keep_going = true;

  return keep_going;
}
