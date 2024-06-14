#pragma once
#include "FileOutput.h"
#include "ConsoleOutput.h"
#include <vtkm/cont/ArrayHandleGroupVecVariable.h>

class VACFOutput : public FileOutput
{
public:
  VACFOutput(const Configuration& cfg);
  virtual ~VACFOutput(){};

  void Init() override;
  void Execute() override;
  bool ShouldOutput() override;

  void ExecuteMSD();

protected:
  int _interval;
  Executioner& _executioner;
  std::ofstream _VACF_file;

  ArrayHandle<Id3> temp_position_flag;
  vtkm::cont::ArrayHandle<vtkm::Vec3f> _original_position;
  vtkm::Vec4f _MSD_value_ave;

  vtkm::cont::ArrayHandle<vtkm::Vec3f> _original_velocity;
  vtkm::Vec4f _VACF_value_ave;

private:
  Real _Vlength;
  ArrayHandle<Vec3f> _position;
  ArrayHandle<vtkm::Vec3f> _velocity;
  IdComponent _start_step;
  IdComponent _end_step;
};
