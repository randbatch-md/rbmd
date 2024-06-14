#include "FileOutput.h"
#include "Application.h"
#include "Executioner.h"

FileOutput::FileOutput(const Configuration& cfg)
  : Output(cfg)
  , _system(*(_app.GetSystem()))
  , _executioner(*(_app.GetExecutioner()))
  , _interval(Get<int>("interval", 1))
{
}

void FileOutput::Init() 
{
  Output::Init();
}

void FileOutput::Execute() {}

bool FileOutput::ShouldOutput()
{
  auto current_step = _executioner.CurrentStep();
  if (current_step == _executioner.NumStep())
    return true;
  else
    return current_step % _interval == 0;
}
