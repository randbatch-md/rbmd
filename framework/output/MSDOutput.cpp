#include "Application.h"
#include "MSDOutput.h"
#include "Executioner.h"
#include <vtkm/cont/Algorithm.h>
#include "output/worklet/OutPutWorklet.h"
#include "MeshFreeSystem.h"
#include "ConsoleOutput.h"
#include "FieldName.h"
//
//RegisterObject(MSDOutput);

MSDOutput::MSDOutput(const Configuration& cfg)
  : FileOutput(cfg)
  , _interval(Get<int>("interval", 1))
  , _executioner(*(_app.GetExecutioner()))
  , _MSD_file("msd.rbmd")
  , _start_step(Get<IdComponent>("start_step"))
  , _end_step(Get<IdComponent>("end_step"))
{
}

void MSDOutput::Init() 
{  
  _position = _para.GetFieldAsArrayHandle<vtkm::Vec3f>(field::position);
  _original_position.AllocateAndFill(_position.GetNumberOfValues(), 0);
  temp_position_flag.AllocateAndFill(_position.GetNumberOfValues(), 0);
  _temp_MSD_position.Allocate(_position.GetNumberOfValues());
 
  _Vlength = _para.GetParameter<Real>(PARA_VLENGTH);
  _box = _para.GetParameter<Vec3f>(PARA_BOX);
}

void MSDOutput::Execute() 
{
  if (_executioner.CurrentStep() == _start_step)
  {
    _original_position.DeepCopyFrom(_position);
    _MSD_position.DeepCopyFrom(_position);
    _temp_MSD_position.DeepCopyFrom(_position);
    auto&& position_flag = _para.GetFieldAsArrayHandle<Id3>(field::position_flag);

    _MSD_value_ave = { 0, 0, 0, 0 };
  }

  if (_executioner.CurrentStep() >= _start_step && _executioner.CurrentStep() <= _end_step)
  {
    ExecuteMSD();
  }

  if (_executioner.CurrentStep() == 0)
  {
    _MSD_file << "Step"
              << " "
              << "Time"
              << " "
              << "MSD_x"
              << " "
              << "MSD_y"
              << " "
              << "MSD_z"
              << " "
              << "MSD_total"
              << " "
              << "sqrt(MSD_total) " << std::endl;
  }

  if (ShouldOutput())
  {
    try
    {
      _MSD_file << _executioner.CurrentStep() << " "
                << _executioner.CurrentStep() * _para.GetParameter<Real>(PARA_TIMESTEP) << " "
                << _MSD_value_ave[0] << " " << _MSD_value_ave[1] << " " << _MSD_value_ave[2] << " "
                << _MSD_value_ave[3] << " " << vtkm::Sqrt(_MSD_value_ave[3]) << std::endl;
    }
    catch (const std::exception& e)
    {
      _MSD_file.close();
      console::Error(e.what());
    }
  }
}

bool MSDOutput::ShouldOutput()
{
  if (_executioner.CurrentStep() < 1)
  {
    return false;
  }
  return _executioner.CurrentStep() % _interval == 0;
}

void MSDOutput::ExecuteMSD() 
{
  auto n = _position.GetNumberOfValues();  
  vtkm::cont::ArrayHandle<vtkm::Vec4f> MSD_vector;
  
  auto&& position_flag = _para.GetFieldAsArrayHandle<Id3>(field::position_flag);

  if (_executioner.CurrentStep() == _start_step)
  {
    position_flag.AllocateAndFill(_position.GetNumberOfValues(),0);
  }
  MSD_vector.AllocateAndFill(_position.GetNumberOfValues(),0);
  OutPut::ComputeMSD(_box, _original_position, _position, position_flag, MSD_vector);

  vtkm::Vec4f MSD_value_total = vtkm::cont::Algorithm::Reduce(MSD_vector, vtkm::TypeTraits<vtkm::Vec4f>::ZeroInitialization());  
  _MSD_value_ave = MSD_value_total / n;

  _temp_MSD_position.DeepCopyFrom(_MSD_position);  
}
