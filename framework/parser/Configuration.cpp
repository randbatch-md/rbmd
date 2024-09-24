#include "Configuration.h"

bool Configuration::Have(const std::string& name) const
{
  return HaveInSelf(name) || HaveInJsonNode(name);
}

bool Configuration::HaveInSelf(const std::string& name) const
{
  return _parameters.find(name) != _parameters.end();
}

bool Configuration::HaveInJsonNode(const std::string& name) const
{
  if (HaveInSelf("_json_node"))
  {
    Json::Value* json_node = linb::any_cast<Json::Value*>(_parameters.at("_json_node"));
    return json_node->isMember(name);
  }

  return false;
}
