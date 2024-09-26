//==================================================================================
//  RBMD 1.0 is developed for random batch molecular dynamics calculation.
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
//  Contact Email : [your - email@example.com]
//==================================================================================

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
  Vec3f _box;
  vtkm::cont::ArrayHandle<Vec3f> _position;
  IdComponent _start_step;
  IdComponent _end_step;
};
