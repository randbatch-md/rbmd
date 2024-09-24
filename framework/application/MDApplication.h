#pragma once
#include <string>
#include "Application.h"

class MDApplication : public Application
{
public:
  MDApplication(int argc, char** argv);

  void PrintLogo() override;

  void Run() override;
  void RunExecutioner();
  void CreateActions() override;
  void CreateCommand();
  void InitConfigurationCommand();
  void HyperParametersCommand();
  void ExecutionCommand();
  void OutputsCommand();

private:
  std::string _ifile;
};