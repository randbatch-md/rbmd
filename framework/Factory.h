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

#include "Logging.h"
#include <functional>
#include <map>
#include <memory>
#include <type_traits>

template<typename Base, typename... Args>
struct Factory
{
  template<typename Derived>
  static bool Register(const std::string& name)
  {
    static_assert(std::is_base_of<Base, Derived>::value, "");
    if (data().count(name) != 0)
    {
      console::Warning(name, "已经注册");
      return false;
    }
    data()[name] = [](Args&&... args) -> std::shared_ptr<Base> {
      return std::make_shared<Derived>(std::forward<Args>(args)...);
    };

    return true;
  }

  static std::shared_ptr<Base> Make(const std::string& name, Args&&... args)
  {
    if (data().count(name) == 0)
    {
      console::Error(name, "未注册");
    }

    return data().at(name)(std::forward<Args>(args)...);
  }

  template<typename Derived>
  static std::shared_ptr<Derived> Make(const std::string& name, Args&&... args)
  {
    auto base = Make(name, std::forward<Args>(args)...);
    auto derived = std::dynamic_pointer_cast<Derived>(base);
    if (!derived)
    {
      console::Error(name, "转换失败");
    }
    return derived;
  }
  //using ConstructionFunc = Base* (*)(Args&&...);
  //using ConstructionFunc = std::shared_ptr<Base>(Args&&...);
  using ConstructionFunc = std::function<std::shared_ptr<Base>(Args&&...)>;
  using MapType = std::map<std::string, ConstructionFunc>;

  template<typename Derived>
  struct Registar
  {
    static volatile bool registered;
  };

//private:
  static MapType& data()
  {
    static MapType s;
    return s;
  }
};


