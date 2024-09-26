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

#pragma once

#include "Any.h"
#include "JsonParser.h"
#include <map>
#include <type_traits>
#include "JsonParser.h"

using linb::any_cast;
struct PointType
{
};

struct ValueType
{
};

template<typename T>
struct Trait
{
  using type = T;
};

class Configuration
{
  friend class Object;

public:
  Configuration(){};
  ~Configuration(){};


  template<typename T>
  void Add(const std::string& name, const T& value)
  {
    if (Have(name))
    {
      console::Error(name, "已经存在");
    }
    _parameters[name] = value;
  }


  template<typename T>
  std::enable_if_t<std::is_pointer<T>::value, T> Get_(const std::string& name) const
  {
    if (HaveInSelf(name))
    {
      T pointer = linb::any_cast<T>(_parameters.at(name));
      if (pointer == nullptr)
      {
        console::Error("空指针");
      }
      return pointer;
    }
    console::Error(name, "不存在");
  }

  template<typename T>
  std::enable_if_t<!std::is_pointer<T>::value, T> Get_(const std::string& name) const
  {
    if (HaveInJsonNode(name))
    {
      Json::Value* json_node = linb::any_cast<Json::Value*>(_parameters.at("_json_node"));
      return (*json_node)[name].as<T>();
    }

    console::Error(name, "不存在");
  }

  // 添加用来数据判断
  template<typename T>
  std::enable_if_t<!std::is_pointer<T>::value, T> Get_Judge(const std::string& name) const
  {
    if (HaveInJsonNode(name))
    {
      Json::Value* json_node = linb::any_cast<Json::Value*>(_parameters.at("_json_node"));
      return (*json_node)[name].as<T>();
    }
    else
    {
      return "null";
    }
  }

  template<typename T>
  T Get(const std::string& name) const
  {
    return Get_<T>(name);
  }

  // 添加用来数据判断
  template<typename T>
  T GetJudge(const std::string& name) const
  {
    return Get_Judge<T>(name);
  }

  template<typename T>
  T Get(const std::string& name, const T& value) const
  {
    if (HaveInJsonNode(name))
    {
      Json::Value* json_node = linb::any_cast<Json::Value*>(_parameters.at("_json_node"));
      return (*json_node)[name].as<T>();
    }
    else
    {
      return value;
    }

    console::Error(name, "不存在");
  }

  template<typename T>
  std::vector<T> GetVectorValue(const std::string& name) const
  {
    if (HaveInJsonNode(name))
    {
      Json::Value* json_node = linb::any_cast<Json::Value*>(_parameters.at("_json_node"));
      auto& node = (*json_node)[name];
      if (node.isArray())
      {
        std::vector<T> value;
        for (auto i = 0; i < node.size(); ++i)
        {
          value.push_back(node[i].as<T>());
        }
        return value;
      }
      else
      {
        console::Error(name, "不是数组");
      }
    }

    console::Error(name, "不存在");
  }

  template<typename T>
  std::vector<std::vector<T>> GetVectorOfVectorsValue(const std::string& name) const
  {
    if (HaveInJsonNode(name))
    {
      Json::Value* json_node = linb::any_cast<Json::Value*>(_parameters.at("_json_node"));
      auto& node = (*json_node)[name];
      if (node.isArray())
      {
        std::vector<std::vector<T>> value;
        for (auto i = 0; i < node.size(); ++i)
        {
          if (node[i].isArray())
          {
            std::vector<T> innerValue;
            for (auto j = 0; j < node[i].size(); ++j)
            {
              innerValue.push_back(node[i][j].as<T>());
            }
            value.push_back(innerValue);
          }
          else
          {
            console::Error(name, "内部元素不是数组");
            return {}; // 返回空的 vector
          }
        }
        return value;
      }
      else
      {
        console::Error(name, "不是数组");
      }
    }

    console::Error(name, "不存在");
    return {}; // 返回空的 vector
  }

  // 添加用来数据判断
  template<typename T>
  std::vector<T> GetVectorValueJudge(const std::string& name) const
  {
    if (HaveInJsonNode(name))
    {
      Json::Value* json_node = linb::any_cast<Json::Value*>(_parameters.at("_json_node"));
      auto& node = (*json_node)[name];
      if (node.isArray())
      {
        std::vector<T> value;
        for (auto i = 0; i < node.size(); ++i)
        {
          value.push_back(node[i].as<T>());
        }
        return value;
      }
    }
    else
    {
      return std::vector<T>();
    }
  }

  bool Have(const std::string& name) const;
  bool HaveInSelf(const std::string& name) const;
  bool HaveInJsonNode(const std::string& name) const;

private:
  std::map<std::string, linb::any> _parameters;
};
