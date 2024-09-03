
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
  _parameter = std::make_shared<Para>(); 
  _parameter->SetParameter(PARA_FAR_FORCE, false);
  _parameter->SetParameter(PARA_DIHEDRALS_FORCE, false);
}

void MDApplication::PrintLogo()
{
  std::string logo = R"( SEMD )";

  std::ofstream log_file("rbmd.log");
  try
  {
    log_file << "RBMD (V2.0)" << std::endl
             << "========================================================" << std::endl
             << _parser->GetFileStr() << std::endl;
    log_file.close();
  }
  catch (const std::exception& e)
  {
    log_file.close();
    console::Error(e.what());
  }
}
void MDApplication::Run()
{
  PrintLogo();
  ParseCLI();
  CreateCommand();
  SetupDevice();

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
}

void MDApplication::CreateCommand()
{
  HyperParametersCommand();
  InitConfigurationCommand();
  ExecutionCommand();
  OutputsCommand();
}

void MDApplication::InitConfigurationCommand()
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
      _init_condition = std::make_shared<ModelFileInitCondition>(cfg);
      _parameter->SetParameter(PARA_INIT_WAY, (std::string)"read_data");   
    }
    else
    {
      _init_condition = std::make_shared<LJInitCondition>(cfg); 
      _parameter->SetParameter(PARA_INIT_WAY, (std::string) "inbuild");     
    }
  }
}

void MDApplication::HyperParametersCommand()
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

void MDApplication::ExecutionCommand()
{
  Configuration cfg;
  auto& execution_node = _parser->GetJsonNode("execution");
  cfg.Add<Application*>("_app", this);
  cfg.Add<Json::Value*>("_json_node", &execution_node); 

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
    _run = std::make_shared<ExecutionNVT>(cfg);
  }
  else if (ensemble == "TEST")
  {
  }
  else
  {
    std::cout << "Json File Error:the ensemble of execution is unknown" << std::endl;
    exit(0);
  }
  _executioner = std::make_shared<Executioner>(cfg);
}

void MDApplication::OutputsCommand()
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