#include "ConsoleOutput.h"
#include "Executioner.h"

//RegisterObject(ConsoleOutput);
const std::string HEADER_STEP_NAME = "STEP";
const std::string HEADER_TIME_NAME = "TIME";
const std::string HEADER_CUMULATIVE_TIME_NAME = "CUMULATIVE_TIME";
ConsoleOutput::ConsoleOutput(const Configuration& cfg) : 
	Output(cfg)
  , _name(Get<std::string>("type"))
  , _timer(_app.GetSystem()->GetTimer())
  , _cumulative_time(0)
  , _output_screen(Get<bool>("output_screen"))
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
  if (_output_screen)
  {
    if (_format_table->IsEmpty())
      return;
    EssientialOutput();
    _format_table->Print();
    _format_table->LogFile();
  }
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
