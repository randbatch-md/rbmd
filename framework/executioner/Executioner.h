//==================================================================================
//  RBMD 1.0 is developed for random batch molecular dynamics calculation.
//
//  Copyright(C) 2024 SHANGHAI JIAOTONG UNIVERSITY CHONGQING RESEARCH INSTITUTE
//
//  This program is free software : you can redistribute it and /or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with this program.If not, see < https://www.gnu.org/licenses/>.
//
//  Contact Email : [your - email@example.com]
//==================================================================================

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
  int _num_steps;
  int _current_step;
  Real _dt, _current_time;
  Real _start_time, _end_time;
};