#pragma once
#include "Application.h"
#include <vtkm/Types.h>
#include "vtkm/cont/ArrayHandle.h"
#include <any>
#include <vtkm/cont/Field.h>
#include "Types.h"

class Input : public Application
{
public:
  Input(int argc, char** argv);
  virtual ~Input()=default;
  //void ExecuteCommand(const Object& obj);
  void InitConfigurationInBuild();
  void InitConfigurationReadData();
  void ExecuteCommandom() override;


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

protected:
  //std::shared_ptr<Configuration> _cfg;
  std::map<std::string, vtkm::cont::Field> _field;
  std::map<std::string, std::any> _parameters;
  //Application& _app;

};