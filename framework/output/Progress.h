#pragma once
#include <iostream>
#include <fstream>
#include "FileOutput.h"
#include "Executioner.h"
class Progress : public FileOutput
{
public:
  Progress(const Configuration& cfg);
  ~Progress();

  void Init() override{};
  void Execute() override;

protected:
  std::ofstream _progress;
  int _judgement_step;
};