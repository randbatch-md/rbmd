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

#include "Configuration.h"
#include "Factory.h"
#include "Application.h"
#include <vtkm/Types.h>

class Object// : public Factory<Object, Configuration&>
{
public:
  Object(const Configuration& cfg)
    : _cfg(cfg)
  {
  }

  virtual ~Object() = default;

  template<typename T>
  T Get(const std::string& name) const
  {
    return _cfg.Get<T>(name);
  }

  //添加用来做判断
  template<typename T>
  T GetJudge(const std::string& name) const
  {
    return _cfg.GetJudge<T>(name);
  }

  template<typename T>
  T Get(const std::string& name, const T& value) const
  {
    return _cfg.Get<T>(name, value);
  }

  template<typename T>
  std::vector<T> GetVectorValue(const std::string& name)
  {
    return _cfg.GetVectorValue<T>(name);
  }

  template<typename T>
  std::vector<std::vector<T>> GetVectorOfVectorsValue(const std::string& name)
  {
    return _cfg.GetVectorOfVectorsValue<T>(name);
  }

  //添加用来做判断
  template<typename T>
  std::vector<T> GetVectorValueJudge(const std::string& name)
  {
    return _cfg.GetVectorValueJudge<T>(name);
  }

protected:
  const Configuration& _cfg;
};

