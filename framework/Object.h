﻿#pragma once

#include "Configuration.h"
#include "Factory.h"
#include "Application.h"

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

protected:
  const Configuration& _cfg;
};

//#define RegisterObject(ObjectName)                                                                 \
//  template<>                                                                                       \
//  template<>                                                                                       \
//  volatile bool Object::Registar<ObjectName>::registered =                                         \
//    Object::Register<ObjectName>(#ObjectName)
//using ObjectFactory = Factory<Object, Configuration&>;
//
//#define RegisterObject(ObjectName)                                                                 \
//  volatile bool registered_##ObjectName = ObjectFactory::Register<ObjectName>(#ObjectName)

