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
//  The post-processing data produced by VASPKIT may be used in any publications 
//  provided that its use is explicitly acknowledged. A suitable reference for VASPKIT is:
//  [1] Gao W, Zhao T, Guo Y, et al.RBMD: A molecular dynamics package enabling to simulate 
//  10 million all - atom particles in a single graphics processing unit[J].arXiv preprint arXiv : 2407.09315, 2024.
// 
//  Contact Email : [support_wz@sciai.com.cn]
//==================================================================================

#include "Executioner.h"
#include "TrajectoryOutput.h"
#include "Application.h"
#include <vtkm/cont/Algorithm.h>
#include "FieldName.h"

TrajectoryOutput::TrajectoryOutput(const Configuration& cfg)
  : FileOutput(cfg)
  , _trajectory_file("rbmd.trj", std::ios::ate)
  , _interval(Get<int>("interval", 1))
{
}

void TrajectoryOutput::Init() 
{
  FileOutput::Init();
}

void TrajectoryOutput::Execute()
{
  if (ShouldOutput())
  {
    try
    {
      auto position = _para.GetFieldAsArrayHandle<Vec3f>(field::position);
      auto portal_pos = position.ReadPortal();
      vtkm::Vec<vtkm::Range, 3> range = _para.GetParameter<vtkm::Vec<vtkm::Range, 3>>(PARA_RANGE);
      auto atom_type = _para.GetFieldAsArrayHandle<Id>(field::pts_type);
      auto portal_atomtype = atom_type.ReadPortal();
      _trajectory_file << "ITEM: TIMESTEP" << std::endl
                       << _executioner.CurrentStep() << std::endl
                       << "ITEM: NUMBER OF ATOMS" << std::endl
                       << position.GetNumberOfValues() << std::endl
                       << "ITEM: BOX BOUNDS pp pp pp" << std::endl
                       << range[0].Min << " " << range[0].Max << std::endl
                       << range[1].Min << " " << range[1].Max << std::endl
                       << range[2].Min << " " << range[2].Max << std::endl
                       << "ITEM: ATOMS id type x y z" << std::endl;
      for (auto i = 0; i < position.GetNumberOfValues(); ++i)
      {
        _trajectory_file << i + 1 << " " << portal_atomtype.Get(i) + 1 << " "
                         << portal_pos.Get(i)[0] << " " << portal_pos.Get(i)[1] << " "
                         << portal_pos.Get(i)[2] << std::endl;
      }
    }
    catch (const std::exception& e)
    {
      _trajectory_file.close();
      console::Error(e.what());
    }
  }
}

bool TrajectoryOutput::ShouldOutput()
{
  if (_executioner.CurrentStep() < 1)
  {
    return false;
  }
  return _executioner.CurrentStep() % _interval == 0;
}