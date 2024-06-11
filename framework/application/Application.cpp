#include "Application.h"
#include "CommandLine.h"
//#include "Register.h"
#include <vtkm/cont/RuntimeDeviceTracker.h>
#include "InitGlobal.h"
#include "SetupDeviceAction.h"

Application::Application(int argc, char** argv)
{
  int iargc = 1;
  if (argc == 2)
  {
    if (std::string(argv[iargc]) == "--version")
    {
      OutputVersion();
    }
    else if (std::string(argv[iargc]) == "-help")
    {
      HelpMessages();
    }
    else
    {
      ErrerMessages();
    }
  }
  else if (argc > 2)
  {
    if (std::string(argv[iargc]) == "-j" && argc == 3)
    {
      _command_line = std::make_unique<CommandLine>(argc, argv);
      _init_global = std::make_unique<InitGlobal>(argc, argv);
      //RegisterObjectGlobal();
    }
    else if (std::string(argv[iargc]) == "-j" && argc > 3)
    {
      std::cout << "--- The formate is: '-j'+ '*.json'---" << std::endl;
      exit(0);
    }
    else
    {
      ErrerMessages();
    }
  }
  else
  {
    ErrerMessages();
  }
  //_command_line = std::make_unique<CommandLine>(argc, argv);
  //_init_global = std::make_unique<InitGlobal>(argc, argv);
  //RegisterObjectGlobal();
}

Application::~Application() {}

void Application::SetupDevice()
{
  auto& tracker = vtkm::cont::GetRuntimeDeviceTracker();

  try
  {
    tracker.ForceDevice(*_device);

    if (!tracker.CanRunOn(*_device))
    {
      console::Error(
        "不能在 Device Tag: ", _device->GetName(), "上运行，", "选项：serial|cuda|tbb|openmp|hip");
    }
  }
  catch (const std::exception&)
  {
    console::Error(
      "不能在 Device Tag: ", _device->GetName(), "上运行，", "选项：serial|cuda|tbb|openmp|hip");
  }

  console::Info("Device Tag: ", _device->GetName());
}

  void Application::OutputVersion()
{
  const std::string RBMD_VERSION = "1.0.0";
  std::cout << RBMD_VERSION << std::endl;
  exit(0);
}

void Application::HelpMessages()
{
  std::cout << "---Reference website: https://www.randbatch.com/guide/CASE_STUDIES.html---"
            << std::endl;
  exit(0);
}

void Application::ErrerMessages()
{
  std::cout << "---Please enter the configuration command---" << std::endl
            << "---You can enter '-help' to query---" << std::endl;
  exit(0);
}
void Application::Run()
{

  //vtkm::cont::ScopedRuntimeDeviceTracker track(tbb);

  PrintLogo(); // 打印logo

  ParseCLI(); //解析命令行

  CreateActions(); // 解析配置文件，创建Actions
}
void Application::AddOutput(std::shared_ptr<Output> output)
{
  _owh.push_back(output);
}
void Application::AddInitCondition(std::shared_ptr<InitCondition> initCondition)
{
  _init_condition_wh.push_back(initCondition);
}
void Application::PrintLogo()
{
  auto logo = R"(
                                                                                       
   @@@@@@@@@     .@@.     =@@@@@@]        @@`      @\       @ 
         =@`     @/=@     =@     \@`     //=@      @\@.     @ 
        /@.     =@  @\    =@     =@^    =@. \\     @ ,@`    @ 
      ,@/      =@`  ,@^   =@   ,/@/    ,@`   @^    @  .@\   @ 
     /@`      ,@\]]]]/@`  =@[[\@`      @\]]]]/@`   @    \@. @ 
   ,@/        @/      \@  =@   ,@^    /@      =@.  @     =@`@ 
  =@\]]]]]]] /@        @\ =@    .@\  =@.       @\  @      ,@@ 
  ---------------------------------------------------------------

)";
  std::cout << logo;
}

void Application::ParseCLI() {}

void Application::CreateActions() {}
