#pragma once
#include "Object.h"
#include "Types.h"
#include "Timer.h"
#include "Execution.h"
#include "Para.h"

class Application;

class Executioner : public Object
{
public:
  Executioner(const Configuration& cfg);

  ~Executioner(){};

  void Init();
  void PreExecute(){};
  void Execute();
  void PostExecute();

  Real& Dt() { return _dt; }
  int& CurrentStep() { return _current_step; }
  int& NumStep() { return _num_steps; }
  Real& CurrentTime() { return _current_time; }
  Real& EndTime() { return _end_time; }
  
private:
  void PreStep();
  void ComputeDT(){};
  void TakeStep();
  void PostStep();
  bool KeepGoing();

protected:
  Application& _app;
  Para& _para;
  Execution& _run;
  Timer& _timer;
  Timer _time_count;
  int _num_steps;
  int _current_step;
  Real _dt, _current_time;
  Real _start_time, _end_time;
};