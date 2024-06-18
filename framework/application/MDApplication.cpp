
#include "MDApplication.h"
#include "Executioner.h"
#include "CommandLine.h"
#include "ModelFileInitCondition.h"
#include "LJInitCondition.h"
#include "ForceField.h"
#include "Neighbor.h"
#include "Coulomb.h"
#include "Extend.h"
#include "ExecutionNPT.h"
#include "ExecutionNVE.h"
#include "ExecutionNVT.h"

#include "H2OSystem.h"
#include "ThermoOutput.h"
#include "TempOutput.h"
#include "RDFOutput.h"
#include "MSDOutput.h"
#include "VACFOutput.h"
#include "TrajectoryOutput.h"
MDApplication::MDApplication(int argc, char** argv)
  : Application(argc, argv)
{
  _ifile = _command_line->GetValue<std::string>("j");
  _parser = std::make_shared<JsonParser>(_ifile);
  _cfg = std::make_shared<Configuration>();
  _parameter = std::make_shared<Para>();   // 这里要提前将app中的Para生成好，后面才会存在一个para中；
  _parameter->SetParameter(PARA_FAR_FORCE, false);
  _parameter->SetParameter(PARA_DIHEDRALS_FORCE, false);
}

void MDApplication::PrintLogo()
{
  std::string logo = R"(
SEMD
)";
  //std::cout << logo << std::endl;
}
void MDApplication::Run()
{
  PrintLogo();
  ParseCLI();
  CreateCommandom();
  //CreateActions();
  //SetupDevice();

  RunExecutioner();
}

void MDApplication::RunExecutioner()
{
  _init_condition->Execute();
  _executioner->Init();
  _executioner->Execute();
}

void MDApplication::CreateActions()
{
 // _awh.push_back(std::make_shared<SetupDeviceAction>(*this));
 // _awh.push_back(std::make_shared<CreateSystemAction>(*this));
 // _awh.push_back(std::make_shared<CreateInitConditionAction>(*this));
 // _awh.push_back(std::make_shared<CreateExecutionerAction>(*this));
 // _awh.push_back(std::make_shared<AddOutputAction>(*this));
 //
 // for (auto& action : _awh)
 // {
 //   action->Execute();
 // }
}

void MDApplication::CreateCommandom()
{
  // 初始化配置文件节点
  HyperParametersCommandom();
  InitConfigurationCommandom();
  //HyperParametersCommandom();
  ExecutionCommandom();
  OutputsCommandom();
}

void MDApplication::InitConfigurationCommandom()
{
  auto& init_node = _parser->GetJsonNode("init_configuration");
  std::vector<std::string> inits;
  if (init_node.isObject())
  {
    inits = init_node.getMemberNames();
    auto init_string = inits[0].c_str();
    auto& init_child_node = init_node[init_string];
    Configuration cfg;
    cfg.Add<Application*>("_app", this);
    cfg.Add<Json::Value*>("_json_node", &init_child_node);
    if (inits[0] == "read_data")
    {
      //_cfg->Add<Json::Value*>("read_data", &init_child_node);
      //cfg.Add<Application*>("_app", this);
      //cfg.Add<Json::Value*>("_json_node", &init_child_node);
      _init_condition = std::make_shared<ModelFileInitCondition>(cfg);
      //_init_condition->Execute();
      _parameter->SetParameter(PARA_INIT_WAY, (std::string)"read_data");   


      //_parameters = std::make_shared<ReadData>(cfg);
      //_parameters->Execute();
    }
    else
    {
      //_cfg->Add<Json::Value*>("inbuild", &init_child_node);
      //cfg.Add<Json::Value*>("_json_node", &init_child_node);
      //cfg.Add<Application*>("_app", this);
      _init_condition = std::make_shared<LJInitCondition>(cfg); // 调用一次这个就可以初始化
      //_init_condition->Execute(); // 这里是初始化了各个参数；
      _parameter->SetParameter(PARA_INIT_WAY, (std::string) "inbuild");      
      auto _init_way = _parameter->GetParameter<std::string>(PARA_INIT_WAY);


      //_init_condition->InitParameter();
      //_init_condition->Execute();
      //_input->InitConfigurationInBuild();

      //_parameters = std::make_shared<InBuild>(cfg);
      //_parameters->Execute();
    }
  }
}

void MDApplication::HyperParametersCommandom()
{
  auto& parameter_node = _parser->GetJsonNode("hyper_parameters");
  std::vector<std::string> parameters;
  parameters = parameter_node.getMemberNames();
  for (const auto& parameter : parameters)
  {
    Configuration cfg;
    auto parameter_string = parameter.c_str();
    auto& parameter_child_node = parameter_node[parameter_string];
    cfg.Add<Application*>("_app", this);
    cfg.Add<Json::Value*>("_json_node", &parameter_child_node);
    if (parameter == "force_field")
    {
      _hyper_parameters = std::make_shared<ForceField>(cfg);
      _hyper_parameters->Execute();
    }
    else if (parameter == "neighbor")
    {
      _hyper_parameters = std::make_shared<Neighbor>(cfg);
      _hyper_parameters->Execute();
    }
    else if (parameter == "coulomb")
    {
      _hyper_parameters = std::make_shared<Coulomb>(cfg);
      _hyper_parameters->Execute();
      _parameter->SetParameter(PARA_FAR_FORCE, true);
    }
    else if (parameter == "extend")
    {
      _hyper_parameters = std::make_shared<Extend>(cfg);
      _hyper_parameters->Execute();

      if (!_parameter->GetParameter<std::vector<int>>(PARA_SPECIAL_BONDS).empty())
        _parameter->SetParameter(PARA_DIHEDRALS_FORCE, true);   
    }
    else
      std::cout << "hyper_parameters is wrong" << std::endl;
  }
}

void MDApplication::ExecutionCommandom()
{
  Configuration cfg;
  auto& execution_node = _parser->GetJsonNode("execution");
  //_cfg->Add<Json::Value*>("execution", &execution_node);
  cfg.Add<Application*>("_app", this);
  cfg.Add<Json::Value*>("_json_node", &execution_node); // 为了提取一级标题的内容；

  auto ensemble = cfg.Get<std::string>("ensemble");
  if (ensemble == "NPT")
  {
    _run = std::make_shared<ExecutionNPT>(cfg);
  }
  else if (ensemble == "NVE")
  {
    _run = std::make_shared<ExecutionNVE>(cfg);
  }
  else if (ensemble == "NVT")
  {
    //_run = std::make_shared<ExecutionTest>(cfg);
    _run = std::make_shared<ExecutionNVT>(cfg);
  }
  else if (ensemble == "TEST")
  {
    //_run = std::make_shared<ExecutionTest>(cfg);
  }
  else
  {
    std::cout << "Json File Error:the ensemble of execution is unknown" << std::endl;
    exit(0);
  }
  _executioner = std::make_shared<Executioner>(cfg);

  // // Execution中还有2个二级标题（temperature 和 pressure）暂时改为了数组
  // std::vector<std::string> executions;
  // executions = execution_node.getMemberNames();
  // for (auto& execution : executions)
  // {
  //   Configuration cfg;
  //   auto execution_string = execution.c_str();
  //   auto& execution_child_node = execution_node[execution_string];
  //   cfg.Add<Application*>("_app", this); // 这个app 可以注释掉？！！！
  //   cfg.Add<Json::Value*>("_json_node", &execution_child_node);
  // }
}

void MDApplication::OutputsCommandom()
{
  auto& output_node = _parser->GetJsonNode("outputs");
  std::vector<std::string> outputs;
  outputs = output_node.getMemberNames();
  for (const auto& output : outputs)
  {
    Configuration cfg;
    auto output_string = output.c_str();
    auto& output_child_node = output_node[output_string];
    cfg.Add<Application*>("_app", this);
    cfg.Add<Json::Value*>("_json_node", &output_child_node);
    std::cout << output_string << std::endl;
    if (output == "thermo_out")
    {
        _Output = std::make_shared<ThermoOutput>(cfg);
        _owh.push_back(_Output);
        _Output = std::make_shared<TempOutput>(cfg);
        _owh.push_back(_Output);

    }
    else if (output == "rdf_out")
    {
        _Output = std::make_shared<RDFOutput>(cfg);
       _owh.push_back(_Output);
    }
    else if (output == "msd_out")
    {
        _Output = std::make_shared<MSDOutput>(cfg);
        _owh.push_back(_Output);
    }
    else if (output == "vacf_out")
    {
        _Output = std::make_shared<VACFOutput>(cfg);
        _owh.push_back(_Output);
    }
    else if (output == "trajectory_out")
    {
        _Output = std::make_shared<TrajectoryOutput>(cfg);
        _owh.push_back(_Output);
    }
    else
        std::cout << "outputs is wrong" << std::endl;
  }
}