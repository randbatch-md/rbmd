#include "Application.h"
#include "VACFOutput.h"
#include "Executioner.h"
#include <vtkm/cont/Algorithm.h>
#include "output/worklet/OutPutWorklet.h"
#include "ConsoleOutput.h"
#include "FieldName.h"

VACFOutput::VACFOutput(const Configuration& cfg)
  : FileOutput(cfg)
  , _interval(Get<int>("interval", 1))
  , _executioner(*(_app.GetExecutioner()))
  , _VACF_file("vacf.rbmd")
  , _start_step(Get<IdComponent>("start_step"))
  , _end_step(Get<IdComponent>("end_step"))
{
}

void VACFOutput::Init() 
{  
  _position = _para.GetFieldAsArrayHandle<vtkm::Vec3f>(field::position);
  _velocity = _para.GetFieldAsArrayHandle<vtkm::Vec3f>(field::velocity);
  _original_velocity.AllocateAndFill(_position.GetNumberOfValues(), 0);
}

void VACFOutput::Execute() 
{
  if (_executioner.CurrentStep() == _start_step)
  {
    _original_velocity.DeepCopyFrom(_velocity);
  }

  if (_executioner.CurrentStep() >= _start_step && _executioner.CurrentStep() <= _end_step)
  {
    ExecuteMSD();
  }

   if (_executioner.CurrentStep() == 0)
   {
       _VACF_file << "Step"
                  << " "
                  << "Time"
                  << " "
                  << "VACFx"
                  << " "
                  << "VACFy"
                  << " "
                  << "VACFz"
                  << " "
                  << "VACF_total" << std::endl;
    }

  if (ShouldOutput())
  {
    try
    {
      _VACF_file << _executioner.CurrentStep() << " "
                 << _executioner.CurrentStep() * _para.GetParameter<Real>(PARA_TIMESTEP) << " "
                 << _VACF_value_ave[0] << " " << _VACF_value_ave[1] << " " << _VACF_value_ave[2]
                 << " " << _VACF_value_ave[3] << std::endl;
    }
    catch (const std::exception& e)
    {
      _VACF_file.close();
      console::Error(e.what());
    }
  }
}

bool VACFOutput::ShouldOutput()
{
  if (_executioner.CurrentStep() < 1)
  {
    return false;
  }
  return _executioner.CurrentStep() % _interval == 0;
}

void VACFOutput::ExecuteMSD() 
{
  auto n = _position.GetNumberOfValues(); 

  vtkm::cont::ArrayHandle<vtkm::Vec4f> VACF_vector;
  VACF_vector.AllocateAndFill(n,0);
  
  OutPut::ComputeVACF(_original_velocity, _velocity, VACF_vector);

  vtkm::Vec4f VACF_value_total = vtkm::cont::Algorithm::Reduce(VACF_vector, vtkm::TypeTraits<vtkm::Vec4f>::ZeroInitialization());
  _VACF_value_ave = VACF_value_total / n;
}
