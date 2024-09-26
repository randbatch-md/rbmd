//==================================================================================
//  RBMD 1.0 is developed for random batch molecular dynamics calculation.
//
//  Copyright(C) 2024 SHANGHAI JIAOTONG UNIVERSITY CHONGQING RESEARCH INSTITUTE
//
//  This program is free software : you can redistribute it and /or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with this program.If not, see < https://www.gnu.org/licenses/>.
//
//  Contact Email : [your - email@example.com]
//==================================================================================

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
