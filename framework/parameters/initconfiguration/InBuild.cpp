#include "FieldName.h"
#include "InBuild.h"
InBuild::InBuild(const Configuration& cfg)
  : InitConfiguration(cfg) 
{
}

void InBuild::Execute()
{
  SetPara();
}

void InBuild::SetPara()
{
  auto dims = GetVectorValue<int>("dims"); 
  auto num_pos = dims[0] * dims[1] * dims[2];
  auto x_range = GetVectorValue<Real>("x_range");
  auto y_range = GetVectorValue<Real>("y_range");
  auto z_range = GetVectorValue<Real>("z_range");
  auto x_length = x_range[1] - x_range[0];
  auto y_length = y_range[1] - y_range[0];
  auto z_length = z_range[1] - z_range[0];
  auto range = vtkm::Vec<vtkm::Range, 3>{ { x_range[0], x_range[1] },{ y_range[0], y_range[1] },{ z_range[0], z_range[1] } };

  _para.SetParameter(PARA_VLENGTH, x_length);
  _para.SetParameter(PARA_VOLUME, x_length * y_length * z_length);
  _para.SetParameter(PARA_RHO, num_pos / (x_length * y_length * z_length));
  _para.SetParameter(PARA_RANGE, range);
  _para.SetParameter(gtest::velocity_type, Get<std::string>("velocity_type"));
}

void InBuild::InitField() 
{
  _para.AddField(field::charge, ArrayHandle<Real>{});
  _para.AddField(field::velocity, ArrayHandle<Vec3f>{});
  _para.AddField(field::mass, ArrayHandle<Real>{});
  _para.AddField(field::molecule_id, ArrayHandle<Id>{});
  _para.AddField(field::atom_id, ArrayHandle<Id>{});
  _para.AddField(field::position, ArrayHandle<Vec3f>{});

  _para.AddField(field::pts_type, ArrayHandle<Id>{});
  _para.AddField(field::position_flag, ArrayHandle<Id3>{});
  _para.AddField(field::center_position, ArrayHandle<Vec3f>{});
  _para.AddField(field::target_position, ArrayHandle<Vec3f>{});
  _para.AddField(field::epsilon, ArrayHandle<Real>{});
  _para.AddField(field::sigma, ArrayHandle<Real>{});
}
