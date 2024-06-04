#include "MeshFreeCondition.h"
#include "worklet/SystemWorklet.h"
#include "math/Utiles.h"
#include "FieldName.h"
#include <vtkm/cont/ArrayCopy.h>
#include <vtkm/cont/ArrayHandleIndex.h>
#include <vtkm/cont/DataSetBuilderUniform.h>

MeshFreeCondition::MeshFreeCondition(const Configuration& cfg)
  : InitCondition(cfg)
  , _dims(GetVectorValue<int>("dims"))
  , _x_range(GetVectorValue<Real>("x_range"))
  , _y_range(GetVectorValue<Real>("y_range"))
  , _z_range(GetVectorValue<Real>("z_range"))
{
}

MeshFreeCondition::~MeshFreeCondition() {}

void MeshFreeCondition::Execute()
{}

void MeshFreeCondition::UpdateField()
{
  InitPosition();
  InitParameters();
  InitMassAndVelocity();
  InitPositionFlag();
  InitId();
}

//void MeshFreeCondition::InitParameter()
//{
//  auto& app = _app.GetSystem();
//  auto xLength = _x_range[1] - _x_range[0];
//  auto yLength = _y_range[1] - _y_range[0];
//  auto zLength = _z_range[1] - _z_range[0];
//  auto num_pos = _dims[0] * _dims[1] * _dims[2];
//   SetParameter(PARA_VLENGTH, xLength);
//   SetParameter(PARA_VOLUME, xLength * yLength * zLength);
//   SetParameter(PARA_RHO, num_pos / (xLength * yLength * zLength));
//  auto cut_off =  GetParameter<Real>(PARA_CUTOFF);
//  auto bin_number = Id3{
//    static_cast<int>(xLength / cut_off),
//    static_cast<int>(yLength / cut_off),
//    static_cast<int>(zLength / cut_off),
//  };
//  auto range = vtkm::Vec<vtkm::Range, 3>{ { _x_range[0], _x_range[1] },
//                                          { _y_range[0], _y_range[1] },
//                                          { _z_range[0], _z_range[1] } };
//   SetParameter(PARA_BIN_NUMBER, bin_number);
//   SetParameter(PARA_RANGE, range);
//   SetParameter(gtest::velocity_type, Get<std::string>("velocity_type"));
//}

void MeshFreeCondition::InitPosition()
{
  _position = _para.GetFieldAsArrayHandle<Vec3f>(field::position);
  Vec3f start{ _x_range[0], _y_range[0], _z_range[0] };
  Vec3f end{ _x_range[1], _y_range[1], _z_range[1] };
  Vec3f space = end - start;
  Id3 point_dim = { _dims[0] + 1, _dims[1] + 1, _dims[2] + 1 };
  for (size_t i = 0; i < 3; i++)
  {
    space[i] = (end[i] - start[i]) / _dims[i];
  }
  auto data_set = vtkm::cont::DataSetBuilderUniform::Create(point_dim, start, space);

  SystemWorklet::InitPosition(data_set.GetCellSet(), data_set.GetCoordinateSystem(), _position);
}

void MeshFreeCondition::InitMassAndVelocity()
{
  ArrayHandle<Vec3f> velocity_temp;
  auto num = _position.GetNumberOfValues();

  velocity_temp.Allocate(num);
  auto write_portal = velocity_temp.WritePortal();
  auto velocity_type = _para.GetParameter<std::string>(gtest::velocity_type);
  if (velocity_type == "GAUSS")
  {
    for (auto i = 0; i < num; ++i)
    {
       Real u10 = RandomValue<Real>(0.0, 1.0);
       Real u20 = RandomValue<Real>(0.0, 1.0);
       Real z00 = vtkm::Sqrt(-2.0 * vtkm::Log(u10)) * vtkm::Cos(2.0 * vtkm::Pif() * u20);
       Real u11 = RandomValue<Real>(0.0, 1.0);
       Real u21 = RandomValue<Real>(0.0, 1.0);
       Real z01 = vtkm::Sqrt(-2.0 * vtkm::Log(u11)) * vtkm::Cos(2.0 * vtkm::Pif() * u21);
       Real u12 = RandomValue<Real>(0.0, 1.0);
       Real u22 = RandomValue<Real>(0.0, 1.0);
       Real z02 = vtkm::Sqrt(-2.0 * vtkm::Log(u12)) * vtkm::Cos(2.0 * vtkm::Pif() * u22);
       Real sigmaz = 1.0;
       Vec3f gaussianr = Vec3f(z00 * sigmaz, z01 * sigmaz, z02 * sigmaz);
       write_portal.Set(i, gaussianr);
    }
  }
  else if (velocity_type == "ZERO" || "TEST")
  {
    for (auto i = 0; i < num; ++i)
    {
      write_portal.Set(i, { 0, 0, 0 });
    }
  }
  auto velocity = _para.GetFieldAsArrayHandle<Vec3f>(field::velocity);
  auto mass = _para.GetFieldAsArrayHandle<Real>(field::mass);
  SystemWorklet::InitCondtion(_position, velocity_temp, mass, velocity);
}

void MeshFreeCondition::InitPositionFlag()
{
  auto num = _position.GetNumberOfValues();

  auto position_flag = _para.GetFieldAsArrayHandle<Id3>(field::position_flag);
  position_flag.AllocateAndFill(num, { 0, 0, 0 });
}

void MeshFreeCondition::InitId()
{
  auto atom_id = _para.GetFieldAsArrayHandle<Id>(field::atom_id);
  auto molecule_id = _para.GetFieldAsArrayHandle<Id>(field::molecule_id);
  vtkm::cont::ArrayHandleIndex atomsIdIndex(_position.GetNumberOfValues());
  vtkm::cont::ArrayCopy(atomsIdIndex, atom_id);
  vtkm::cont::ArrayHandleIndex moleculeIdIndex(_position.GetNumberOfValues());
  vtkm::cont::ArrayCopy(moleculeIdIndex, molecule_id);
}

void MeshFreeCondition::InitParameters()
{
  auto num_pos = _position.GetNumberOfValues();

  //pts type
  auto pts_type = _para.GetFieldAsArrayHandle<Id>(field::pts_type);
  pts_type.AllocateAndFill(num_pos, 0);

  //eps
  auto epsilon = _para.GetFieldAsArrayHandle<Real>(field::epsilon);
  epsilon.AllocateAndFill(1, 1.0);

  //sigma
  auto sigma = _para.GetFieldAsArrayHandle<Real>(field::sigma);
  sigma.AllocateAndFill(1, 1.0);
}

void MeshFreeCondition::InitField()
{
  _para.AddField(field::charge, ArrayHandle<Real>{});
  _para.AddField(field::velocity, ArrayHandle<Vec3f>{});
  _para.AddField(field::mass, ArrayHandle<Real>{});
  _para.AddField(field::molecule_id, ArrayHandle<Id>{});
  _para.AddField(field::atom_id, ArrayHandle<Id>{});
  _para.AddField(field::position, ArrayHandle<Vec3f>{});
}

void MeshFreeCondition::SetParameters() {}
