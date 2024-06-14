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
      int a = _executioner.CurrentStep();
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