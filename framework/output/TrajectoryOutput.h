#pragma once

#include "FileOutput.h"
#include <vtkm/cont/ArrayHandleGroupVecVariable.h>

class TrajectoryOutput : public FileOutput
{
public:
  TrajectoryOutput(const Configuration& cfg);
  virtual ~TrajectoryOutput() { _trajectory_file.close(); };

  void Init() override;
  void Execute() override;

protected:
  bool ShouldOutput() override;

private:
  int _interval;

  vtkm::IdComponent _atom_num;
  std::ofstream _trajectory_file;
};