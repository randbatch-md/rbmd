#pragma once
#include "ConsoleOutput.h"
#include <vtkm/cont/ArrayHandleGroupVecVariable.h>


class TempOutput : public ConsoleOutput
{
public:
  TempOutput(const Configuration& cfg);
  virtual ~TempOutput();

  void Init() override;
  void Execute() override;

private:
  void AddDataToTable();
  void WriteToFile();
  bool ShouldOutput();

private:
  std::ofstream _file;

  int _interval;
  bool _output_file;

  Real _temperature;
  Real _pressure;

};
