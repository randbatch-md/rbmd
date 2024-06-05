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
  , _binary(Get<bool>("binary", false))
  , _file(_name + ".csv")
  , _interval(Get<int>("interval", 1))
  , _out_initial(Get<bool>("out_initial", false))
  , _output_file(Get<bool>("output_file"))
  , _compute(Get<bool>("compute"))
  , _temperature(0.0)
{
}

TempOutput::~TempOutput()
{
  _file.close();
}

void TempOutput::Init()
{
  if (_compute) //_output_screen &&
  {
    ConsoleOutput::Init();
    AddHeader(HEADER_KBT_NAME);

  }
}

void TempOutput::Execute()
{
  if (_compute)
  {
    _temperature = _system.GetParameter<Real>(PARA_TEMPT);

    WriteToFile();
  }

}

bool TempOutput::ShouldOutput()
{
  auto current_step = _executioner->CurrentStep();
  if (current_step == 0)
    return _out_initial;
  else if (current_step == _executioner->NumStep())
    return true;
  else
    return current_step % _interval == 0;
}


void TempOutput::WriteToFile()
{
  if (ShouldOutput() && _output_file)
  {
    if (_executioner->CurrentStep() == 1)
    {
      _file << "Step"
            << ", "
            << "temperature" 
            << std::endl;
    }
    try
    {
      _file << _executioner->CurrentStep() << ", " << _temperature << std::endl;
    }
    catch (const std::exception& e)
    {
      _file.close();
      console::Error(e.what());
    }
  }
}

