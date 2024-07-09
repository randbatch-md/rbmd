#include "LJInitCondition.h"
#include "system/worklet/SystemWorklet.h"
#include "FieldName.h"
#include "UnitFactor.h"

//RegisterObject(LJInitCondition);
LJInitCondition::LJInitCondition(const Configuration& cfg)
  : MeshFreeCondition(cfg)
{
  InitField();
  SetParameters();
}

void LJInitCondition::Execute() 
{
  //InitField();
  //SetParameters();
  AddMoleculeInfo();
  UpdateField();
 // auto step = 0;
 // while (step < _max_steps)
 // {
 //   DoInit();
 //   step++;
 // }
 // //console::Info("Init System End!");
}

void LJInitCondition::UpdateField()
{
  MeshFreeCondition::UpdateField();
}

void LJInitCondition::AddMoleculeInfo() 
{
  auto N = _position.GetNumberOfValues();
  std::vector<Real> special_weights(N, 1.0);
  std::vector<Id> special_offsets(N, 1);

  vtkm::Id offsetSize;
  vtkm::cont::ArrayHandle<vtkm::Id> offsetsArray = vtkm::cont::ConvertNumComponentsToOffsets(
    vtkm::cont::make_ArrayHandle(special_offsets), offsetSize);
  auto special_offsets_array = _para.GetFieldAsArrayHandle<Id>(field::special_offsets);
  vtkm::cont::ArrayCopy(offsetsArray, special_offsets_array);

  auto special_ids_array = _para.GetFieldAsArrayHandle<Id>(field::special_ids);
  vtkm::cont::ArrayHandleIndex atomsIdIndex(_position.GetNumberOfValues());
  vtkm::cont::ArrayCopy(atomsIdIndex, special_ids_array);

  auto special_weights_array = _para.GetFieldAsArrayHandle<Real>(field::special_weights);
  vtkm::cont::ArrayCopy(vtkm::cont::make_ArrayHandle(special_weights), special_weights_array);
}

void LJInitCondition::InitField()
{
  MeshFreeCondition::InitField();
  _para.AddField(field::pts_type, ArrayHandle<Id>{});
  _para.AddField(field::position_flag, ArrayHandle<Id3>{});
  _para.AddField(field::center_position, ArrayHandle<Vec3f>{});
  _para.AddField(field::target_position, ArrayHandle<Vec3f>{});
  _para.AddField(field::epsilon, ArrayHandle<Real>{});
  _para.AddField(field::sigma, ArrayHandle<Real>{});
}

void LJInitCondition::SetParameters()
{
  auto xLength = _x_range[1] - _x_range[0];
  auto yLength = _y_range[1] - _y_range[0];
  auto zLength = _z_range[1] - _z_range[0];
  Vec3f box = { xLength, yLength, zLength };
  auto num_pos = _dims[0] * _dims[1] * _dims[2];
  auto range = vtkm::Vec<vtkm::Range, 3>{ { _x_range[0], _x_range[1] },
                                          { _y_range[0], _y_range[1] },
                                          { _z_range[0], _z_range[1] } };
  _para.SetParameter(PARA_VLENGTH, xLength);
  _para.SetParameter(PARA_BOX, box);
  _para.SetParameter(PARA_VOLUME, box[0] * box[1] * box[2]);
  _para.SetParameter(PARA_RHO, num_pos / (box[0] * box[1] * box[2]));
  _para.SetParameter(PARA_RANGE, range);
  _para.SetParameter(gtest::velocity_type, Get<std::string>("velocity_type"));
  _para.SetParameter(PARA_UNIT, Get<std::string>("unit"));
  _para.SetParameter(ATOM_STYLE, Get<std::string>("atom_style"));
  auto cut_off = _para.GetParameter<Real>(PARA_CUTOFF);
  auto bin_number = Id3{
    static_cast<int>(xLength / cut_off),
    static_cast<int>(yLength / cut_off),
    static_cast<int>(zLength / cut_off),
  };
  _para.SetParameter(PARA_BIN_NUMBER, bin_number);

  // 注意：这里要放到后面填充！！！！！因为cut_off暂时不知道；
  // auto cut_off = GetParameter<Real>(PARA_CUTOFF);
  // auto bin_number = Id3{
  //   static_cast<int>(xLength / cut_off),
  //   static_cast<int>(yLength / cut_off),
  //   static_cast<int>(zLength / cut_off),
  // };
  // SetParameter(PARA_BIN_NUMBER, bin_number);
}

/*
void LJInitCondition::DoInit() 
{

  // stage1:
  ComputeForce();
  UpdateVelocity();

  // stage2:
  UpdatePosition();

  // stage3:
  ComputeForce();
  UpdateVelocity();
}

void LJInitCondition::ComputeForce()
{
  auto pts_type = _para.GetFieldAsArrayHandle<Id>(field::pts_type);
  auto atom_id = _para.GetFieldAsArrayHandle<Id>(field::atom_id);
  try
  {
    ContForceFunction force_function;

    ContTopology topology;
    SetTopology(topology);

    ContPointLocator locator;
    SetLocator(locator);
    auto cut_off = _para.GetParameter<Real>(PARA_CUTOFF);
    SystemWorklet::LJForceWithPeriodicBC(cut_off, atom_id, locator, topology, force_function, _LJforce);
  }
  catch (const std::exception& e)
  {
    std::cout << e.what() << std::endl;
  }
}

void LJInitCondition::UpdateVelocity() 
{
  try
  {
    auto velocity = _para.GetFieldAsArrayHandle<Vec3f>(field::velocity);
    auto mass = _para.GetFieldAsArrayHandle<Real>(field::mass);

    auto unit_factor = _para.GetParameter<UnitFactor>(PARA_UNIT_FACTOR);
    SystemWorklet::UpdateVelocity(_dt, unit_factor._fmt2v, _LJforce, mass, velocity);
  }
  catch (const std::exception& e)
  {
    std::cout << e.what() << std::endl;
  }
}

void LJInitCondition::UpdatePosition() 
{
  auto velocity = _para.GetFieldAsArrayHandle<Vec3f>(field::velocity);

  ContPointLocator locator;
  SetLocator(locator);
  SystemWorklet::UpdatePosition(_dt, velocity, locator, _position);
}
*/