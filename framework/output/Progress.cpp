﻿//==================================================================================
//  RBMD 2.2.0 is developed for random batch molecular dynamics calculation.
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

#include "Progress.h"

Progress::Progress(const Configuration& cfg)
  : FileOutput(cfg)
  , _progress("Progress.csv")
  , _judgement_step(_executioner.NumStep() / 10)
{
}

Progress::~Progress() 
{
  _progress.close();
}

void Progress::Execute() 
{
  auto current_step = _executioner.CurrentStep();
  if (current_step==1)
  {
    _progress << 0.05 << std::endl;
  }
  else if (_judgement_step != 0 && current_step % _judgement_step == 0 &&
           current_step != _executioner.NumStep() && current_step > 20)
  {
    try
    {
      _progress << (Real)current_step / _executioner.NumStep() << std::endl;
    }
    catch (const std::exception& e)
    {
      _progress.close();
      console::Error(e.what());
    }
  }
  if (current_step == _executioner.NumStep())
  {
    _progress << 1 << std::endl;
  }
}
