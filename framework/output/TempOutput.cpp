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
//  The post-processing data produced by VASPKIT may be used in any publications 
//  provided that its use is explicitly acknowledged. A suitable reference for VASPKIT is:
//  [1] Gao W, Zhao T, Guo Y, et al.RBMD: A molecular dynamics package enabling to simulate 
//  10 million all - atom particles in a single graphics processing unit[J].arXiv preprint arXiv : 2407.09315, 2024.
// 
//  Contact Email : [support_wz@sciai.com.cn]
//==================================================================================

#include "Application.h"
#include "ConsoleOutput.h"
#include "Executioner.h"
#include "FieldName.h"
#include "TempOutput.h"
#include "UnitFactor.h"
#include "forceFunction/ContForceFunction.h"
#include "math/Math.h"
#include "output/worklet/OutPutWorklet.h"
#include <ctime>
#include <vtkm/cont/Algorithm.h>

const std::string HEADER_KBT_NAME = "KBT";

TempOutput::TempOutput(const Configuration& cfg)
  : ConsoleOutput(cfg)
  , _file("temperature.rbmd")
  , _interval(Get<int>("interval", 1))
  , _temperature(0.0)
{
}

TempOutput::~TempOutput()
{
  _file.close();
}

void TempOutput::Init()
{
  ConsoleOutput::Init();
  AddHeader(HEADER_KBT_NAME);
}

void TempOutput::Execute()
{
  _temperature = _para.GetParameter<Real>(PARA_TEMPT);
  WriteToFile();

}

bool TempOutput::ShouldOutput()
{
  auto current_step = _executioner->CurrentStep();
  if (current_step == _executioner->NumStep())
    return true;
  else
    return current_step % _interval == 0;
}


void TempOutput::WriteToFile()
{
  if (_executioner->CurrentStep() == 0)
  {
    _file << "Step"
          << " "
          << "Temperature" << std::endl;
  }
  if (ShouldOutput())
  {
    try
    {
      _file << _executioner->CurrentStep() << " "
            << _temperature << std::endl;
    }
    catch (const std::exception& e)
    {
      _file.close();
      console::Error(e.what());
    }
  }
}

