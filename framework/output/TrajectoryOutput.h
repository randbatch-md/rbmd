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