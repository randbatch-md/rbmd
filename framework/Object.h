#pragma once
#include <any>
#include <vtkm/cont/ArrayHandle.h>
#include <vtkm/cont/Field.h>
#include "Configuration.h"
#include "Factory.h"
#include "Application.h"
//#include "Types.h"
//#include "DataObject.h"
//#include "Logging.h"

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
    
  /*
  template<typename T>
  void AddField(const std::string& name, const vtkm::cont::ArrayHandle<T>& arrayHandle)
  {
    _field.insert(std::make_pair(
      name, vtkm::cont::Field(name, vtkm::cont::Field::Association::Any, arrayHandle)));
  }

  template<typename T>
  vtkm::cont::ArrayHandle<T> GetFieldAsArrayHandle(const std::string& name)
  {
    try
    {
      return _field[name].GetData().AsArrayHandle<ArrayHandle<T>>();
    }
    catch (const std::exception& e)
    {
      console::Error(e.what());
    }
  }

  template<typename T>
  void SetParameter(const std::string& name, const T& field)
  {
    _parameters[name] = field;
  }

  template<typename T>
  T GetParameter(const std::string& name)
  {
    try
    {
      return std::any_cast<T>(_parameters.at(name));
    }
    catch (const std::exception& e)
    {
      console::Error(e.what());
    }
  }
  */
protected:
  const Configuration& _cfg;
  ////std::shared_ptr<Configuration> _cfg;
  //
  //std::map<std::string, vtkm::cont::Field> _field;
  //std::map<std::string, std::any> _parameters;
};

////#define RegisterObject(ObjectName)                                                                 \
////  template<>                                                                                       \
////  template<>                                                                                       \
////  volatile bool Object::Registar<ObjectName>::registered =                                         \
////    Object::Register<ObjectName>(#ObjectName)
//using ObjectFactory = Factory<Object, Configuration&>;
//
//#define RegisterObject(ObjectName)                                                                 \
//  volatile bool registered_##ObjectName = ObjectFactory::Register<ObjectName>(#ObjectName)

