#include "AddOutputAction.h"
#include "Application.h"
#include "Output.h"

 AddOutputAction::AddOutputAction(Application& app)
  : Action(app)
{

}


void AddOutputAction::Execute()
{
  if (!_parser.HasNode("outputs"))
  {
    return;
  }

  auto& outputs_node = _parser.GetJsonNode("outputs");
  std::vector<std::string> outPuts;
  if (outputs_node.isObject())
  {
    outPuts = outputs_node.getMemberNames();
  }

  auto sys = _app.GetSystem();
  if (!sys)
    console::Error("System未创建.");

  for (const auto& outPut : outPuts)
  {
    auto& output_child_node = outputs_node[outPut.c_str()];
    Configuration cfg;
    cfg.Add<Application*>("app", &_app);
    cfg.Add<Json::Value*>("_json_node", &output_child_node);

    auto type = output_child_node["type"].asString();
    _app.AddOutput(ObjectFactory::Make<Output>(type, cfg));
  }

  //console::Success("AddOutputAction");
}