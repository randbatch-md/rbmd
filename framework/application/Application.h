#pragma once
#include "Object.h"
#include <memory>
#include <vector>

class CommandLine;
class JsonParser;
class System;
class Executioner;
class Action;
class Output;
class ConsoleOutput;
class InitCondition;
class MyTestFixture;
class InitGlobal;

namespace vtkm
{
namespace cont
{
class DeviceAdapterId;
}
}

class Application
{
  friend class MyTestFixture;
public:
  Application(int argc, char** argv);
  virtual ~Application();

  virtual void PrintLogo();

  // 解析命令行
  virtual void ParseCLI();

  virtual void CreateActions();

  virtual void Run();
  void SetupDevice();

  void AddOutput(std::shared_ptr<Output> output);
  void AddInitCondition(std::shared_ptr<InitCondition> initCondition);

  auto &DeviceTag() { return this->_device; }
  auto& GetParser() { return this->_parser; }
  auto& GetSystem() { return this->_system; }
  auto& GetInitCondition() { return this->_init_condition; }
  auto& GetExecutioner() { return this->_executioner; }
  auto& GetOutputWarehouse() { return this->_owh; }
  auto& GetInitConditionWarehouse() { return this->_init_condition_wh; }

  template<typename T>
  T& GetSystem()
  {
    static_assert(std::is_base_of<System, T>::value, "convert error");
    return static_cast<T&>(*_system);
  }

protected:
  std::vector<std::shared_ptr<Action>> _awh;     // action warehouse
  std::vector<std::shared_ptr<Output>> _owh;     // output warehouse
  std::unique_ptr<CommandLine> _command_line;
  std::unique_ptr<InitGlobal> _init_global;
  std::shared_ptr<JsonParser> _parser;
  std::shared_ptr<System> _system;
  std::shared_ptr<Executioner> _executioner;
  std::shared_ptr<InitCondition> _init_condition;
  std::vector<std::shared_ptr<InitCondition>> _init_condition_wh; //init condition ware house

  // 不使用unique_ptr，因为需要#include DeviceAdapterId
  std::shared_ptr<vtkm::cont::DeviceAdapterId> _device;
};
