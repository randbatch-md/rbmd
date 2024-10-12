﻿//==================================================================================
//  RBMD 2.1.0 is developed for random batch molecular dynamics calculation.
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
//  The post-processing data produced by VASPKIT may be used in any publications 
//  provided that its use is explicitly acknowledged. A suitable reference for VASPKIT is:
//  [1] Gao W, Zhao T, Guo Y, et al.RBMD: A molecular dynamics package enabling to simulate 
//  10 million all - atom particles in a single graphics processing unit[J].arXiv preprint arXiv : 2407.09315, 2024.
// 
//  Contact Email : [support_wz@sciai.com.cn]
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
  Timer _time_count;
  int _num_steps;
  int _current_step;
  Real _dt, _current_time;
  Real _start_time, _end_time;
};