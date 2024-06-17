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

//RegisterObject(TempOutput);

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
          << ", "
          << "Temperature" << std::endl;
  }
  if (ShouldOutput())
  {
    try
    {
      _file << _executioner->CurrentStep() << ", "
            << _temperature << std::endl;
    }
    catch (const std::exception& e)
    {
      _file.close();
      console::Error(e.what());
    }
  }
}

