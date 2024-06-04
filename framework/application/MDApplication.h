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
  void CreateCommandom();
  void InitConfigurationCommandom(); // 可以直接在这里初始化一套流程运算
  void HyperParametersCommandom();
  void ExecutionCommandom();
  void OutputsCommandom();
  void ExecuteCommandom() override;

private:
  std::string _ifile;
};