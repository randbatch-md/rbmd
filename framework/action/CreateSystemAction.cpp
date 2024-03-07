#include "Application.h"
#include "CreateSystemAction.h"
#include "System.h"
#include "JsonParser.h"

  CreateSystemAction::CreateSystemAction(Application& app)
  : Action(app)
  {

  };


void CreateSystemAction::Execute()
{
  auto& system_node = _parser.GetJsonNode("system");
  Configuration cfg;
  cfg.Add<Application*>("_app", &_app);
  cfg.Add<Json::Value*>("_json_node", &system_node);

  auto type = system_node["type"].asString();
  _app.GetSystem() = ObjectFactory::Make<System>(type, cfg);
  //console::Success("CreateSystemAction");
}
