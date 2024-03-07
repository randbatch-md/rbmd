#include "MDApplication.h"
#include "AddOutputAction.h"
#include "CreateExecutionerAction.h"
#include "CreateSystemAction.h"
#include "SetupDeviceAction.h"
#include "CreateInitConditionAction.h"
#include "Executioner.h"
#include "CommandLine.h"

MDApplication::MDApplication(int argc, char** argv)
  : Application(argc, argv)
{
  _ifile = _command_line->GetValue<std::string>("j");
  _parser = std::make_shared<JsonParser>(_ifile);
}

void MDApplication::PrintLogo()
{
  std::string logo = R"(
SEMD
)";
  //std::cout << logo << std::endl;
}
void MDApplication::Run()
{
  PrintLogo();
  ParseCLI();
  CreateActions();
  SetupDevice();

  RunExecutioner();
}

void MDApplication::RunExecutioner()
{
  _executioner->Init();
  _executioner->Execute();
}

void MDApplication::CreateActions()
{
  _awh.push_back(std::make_shared<SetupDeviceAction>(*this));
  _awh.push_back(std::make_shared<CreateSystemAction>(*this));
  _awh.push_back(std::make_shared<CreateInitConditionAction>(*this));
  _awh.push_back(std::make_shared<CreateExecutionerAction>(*this));
  _awh.push_back(std::make_shared<AddOutputAction>(*this));

  for (auto& action : _awh)
  {
    action->Execute();
  }
}
