//==================================================================================
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

#include "ConsoleOutput.h"
#include "Executioner.h"

const std::string HEADER_STEP_NAME = "Step";
const std::string HEADER_TIME_NAME = "Time";
const std::string HEADER_CUMULATIVE_TIME_NAME = "Cumulative_Time";
ConsoleOutput::ConsoleOutput(const Configuration& cfg) : 
	Output(cfg)
  , _timer(_app.GetRun()->GetTimer())
  , _cumulative_time(0)
  , _parser(*(_app.GetParser()))
{
  _format_table = std::make_unique<FormatTable>();
  AddHeader(HEADER_STEP_NAME);
  AddHeader(HEADER_TIME_NAME);
  AddHeader(HEADER_CUMULATIVE_TIME_NAME);
}

void ConsoleOutput::Init() 
{
  Output::Init();
}

void ConsoleOutput::Execute() 
{

    if (_format_table->IsEmpty())
      return;
    EssientialOutput();
    _format_table->Print();
    _format_table->LogFile();
}

bool ConsoleOutput::AddHeader(const std::string& name)
{
  return _format_table->AddHeader(name);
}

void ConsoleOutput::EssientialOutput() 
{
  auto& executioner = *(_app.GetExecutioner());
  _time = _timer.GetElapsedTime();
  int current_step = executioner.CurrentStep();
  AddData(HEADER_STEP_NAME, current_step);
  AddData(HEADER_TIME_NAME, _time);
  if (current_step)
  {
    _cumulative_time += _time;
    AddData(HEADER_CUMULATIVE_TIME_NAME, _cumulative_time);
  }
  else
    AddData(HEADER_CUMULATIVE_TIME_NAME, 0);
}
