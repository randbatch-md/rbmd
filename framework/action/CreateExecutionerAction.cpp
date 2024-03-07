#include "Application.h"
#include "CreateExecutionerAction.h"
#include "Executioner.h"
#include "ConsoleOutput.h"

  CreateExecutionerAction::CreateExecutionerAction(Application& app)
  : Action(app)
  {

  };

void CreateExecutionerAction::Execute()
{
  auto& execut_node = _parser.GetJsonNode("executioner");
  Configuration cfg;
  cfg.Add<Application*>("_app", &_app);
  cfg.Add<Json::Value*>("_json_node", &execut_node);

  auto type = "Executioner";
  _app.GetExecutioner() = std::make_shared<Executioner>(cfg);
  //console::Success("CreateExecutioner");
}
