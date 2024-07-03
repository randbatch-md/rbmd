#pragma once
#include "FileOutput.h"
#include "ConsoleOutput.h"
#include <vtkm/cont/ArrayHandleGroupVecVariable.h>
class MSDOutput : public FileOutput
{
public:
  MSDOutput(const Configuration& cfg);
  virtual ~MSDOutput(){};

  void Init() override;
  void Execute() override;
  bool ShouldOutput() override;

  void ExecuteMSD();

protected:
  int _interval;
  Executioner& _executioner;
  std::ofstream _MSD_file;

  vtkm::cont::ArrayHandle<Id3> temp_position_flag;
  vtkm::cont::ArrayHandle<vtkm::Vec3f> _original_position;
  vtkm::Vec4f _MSD_value_ave;

  vtkm::cont::ArrayHandle<Vec3f> _MSD_position;
  vtkm::cont::ArrayHandle<Vec3f> _temp_MSD_position;

private:
  Real _Vlength;
  vtkm::cont::ArrayHandle<Vec3f> _position;
  IdComponent _start_step;
  IdComponent _end_step;
};
