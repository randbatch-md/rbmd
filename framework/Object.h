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

