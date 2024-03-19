#pragma once

#include <any>
#include "DataObject.h"
#include "Logging.h"
#include "Object.h"
#include "Timer.h"

class Application;
class InitCondition;
class System : public Object
{
public:
  System(const Configuration& cfg);
  virtual ~System(){};

  virtual void Init();
  virtual void Evolve() = 0;

  template<typename T>
  void AddField(const std::string& name, const ArrayHandle<T>& arrayHandle)
  {
    _field.insert(std::make_pair(
      name, vtkm::cont::Field(name, vtkm::cont::Field::Association::Any, arrayHandle)));
  }

  template<typename T>
  ArrayHandle<T> GetFieldAsArrayHandle(const std::string& name)
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
  bool HaveParameter(const std::string& name) const
  {
    return _parameters.find(name) != _parameters.end();
  }

    // 定义一个函数用于读取一行数据
  template<typename T>
  void read_line(std::ifstream& file, T& data)
  {
    file >> data;
  }

  // 定义一个函数用于读取一维数组
  template<typename T>
  void read_array(std::ifstream& file, std::vector<T>& array, int size)
  {
    for (int i = 0; i < size; ++i)
    {
      read_line(file, array[i]);
    }
  }

  auto& App() { return this->_app; }
  auto& GetTimer() { return this->_timer; }

protected:
  virtual void InitField(){};

protected:
  std::map<std::string, vtkm::cont::Field> _field;
  std::map<std::string, std::any> _parameters;
  Application& _app;
  Timer _timer;
};
