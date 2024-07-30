#include "Application.h"
#include "parser/InitGlobal.h"
#include "CommandLine.h"
#include <vtkm/cont/RuntimeDeviceTracker.h>

#define VTKM_NO_ERROR_ON_MIXED_CUDA_CXX_TAG
#include <vtkm/cont/DeviceAdapter.h>

#undef VTKM_NO_ERROR_ON_MIXED_CUDA_CXX_TAG


Application::Application(int argc, char** argv)
{
  int iarg = 1;  
  if (std::string(argv[iarg]) == "-help" || std::string(argv[iarg]) == "-h")
    {
      if (2 == argc)
        HelpMessages();      
      else
      {
        ErrerMessages();
      }
    } 
  else if (std::string(argv[iarg]) == "--version" || std::string(argv[iarg]) == "-v")
    {
      if (2 == argc)
        VersionMessages();
      else
      {
        ErrerMessages();
      }
    }   
  else if (argc - iarg >= 2 && std::string(argv[iarg])=="cuda")
    {
      if (4 == argc)
      {
        _device.reset(new vtkm::cont::DeviceAdapterTagCuda);
      }
      else
      {
        ErrerMessages();
      }
    }
  else if (argc - iarg >= 2 && std::string(argv[iarg]) == "dcu")
  {
    if (4 == argc)
    {
      _device.reset(new vtkm::cont::DeviceAdapterTagKokkos);
      _init_global = std::make_unique<InitGlobal>(std::string(argv[iarg]),argc, argv);
    }
    else
    {
      ErrerMessages();
    }
  }
  else if (argc - iarg >= 2 && std::string(argv[iarg]) == "mpi")
  {
    _device.reset(new vtkm::cont::DeviceAdapterTagOpenMP); 
  }
  else if (argc - iarg >= 2 && std::string(argv[iarg]) == "-j")
  {
    if (3 == argc)
    {
      _device.reset(new vtkm::cont::DeviceAdapterTagSerial);
    }
    else
    {
      ErrerMessages();
    }
  }
  else
    {
      std::cout << (argv[iarg]) << std::endl;

      ErrerMessages();
    }
  
  _command_line = std::make_unique<CommandLine>(argc, argv);
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
        "不能在 Device Tag: ", _device->GetName(), "上运行，", "选项：serial|cuda|dcu|mpi");
    }
  }
  catch (const std::exception&)
  {
    console::Error(
      "不能在 Device Tag: ", _device->GetName(), "上运行，", "选项：serial|cuda|dcu|mpi");
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
  exit(1);
}
void Application::VersionMessages()
{
  std::cout << "RBMD_VERSION = 2.0.0 " << std::endl;
  exit(1);
}

void Application::ErrerMessages()
{
  std::cout << "---Please enter the configuration command---" << std::endl
            << "---You can enter '-help' to query---" << std::endl;
  exit(0);
}
void Application::Run()
{
  PrintLogo(); 

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
