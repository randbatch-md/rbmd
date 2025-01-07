//==================================================================================
//  RBMD 2.2.0 is developed for random batch molecular dynamics calculation.
//
//  Copyright(C) 2024 SHANGHAI JIAOTONG UNIVERSITY CHONGQING RESEARCH INSTITUTE
//
//  This program is free software : you can redistribute it and /or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with this program.If not, see < https://www.gnu.org/licenses/>.
//
//  The post-processing data produced by VASPKIT may be used in any publications 
//  provided that its use is explicitly acknowledged. A suitable reference for VASPKIT is:
//  [1] Gao W, Zhao T, Guo Y, et al.RBMD: A molecular dynamics package enabling to simulate 
//  10 million all - atom particles in a single graphics processing unit[J].arXiv preprint arXiv : 2407.09315, 2024.
// 
//  Contact Email : [support_wz@sciai.com.cn]
//==================================================================================

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
