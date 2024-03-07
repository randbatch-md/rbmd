#include "CreateInitConditionAction.h"
#include "Application.h"
#include "Configuration.h"
#include "InitCondition.h"
#include "JsonParser.h"
#include "System.h"

CreateInitConditionAction::CreateInitConditionAction(Application& app)
  : Action(app)
{
}

void CreateInitConditionAction::Execute()
{
  auto& init_node = _parser.GetJsonNode("init_condition");
  std::vector<std::string> initPuts;
  if (init_node.isObject())
  {
    initPuts = init_node.getMemberNames();
  }

  for (const auto& initPut : initPuts)
  {
    auto& initput_child_node = init_node[initPut.c_str()];
    Configuration cfg;
    cfg.Add<Application*>("_app", &_app);
    cfg.Add<Json::Value*>("_json_node", &initput_child_node);

    auto type = initput_child_node["type"].asString();
    _app.AddInitCondition(ObjectFactory::Make<InitCondition>(type, cfg));
  }

  //console::Success("CreateInitConditionAction");
}
