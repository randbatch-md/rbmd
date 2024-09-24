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
