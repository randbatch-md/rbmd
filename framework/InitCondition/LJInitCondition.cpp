#include "LJInitCondition.h"
#include "system/worklet/SystemWorklet.h"
#include "FieldName.h"
#include "UnitFactor.h"

//RegisterObject(LJInitCondition);
LJInitCondition::LJInitCondition(const Configuration& cfg)
  : MeshFreeCondition(cfg)
  , _max_steps(Get<vtkm::IdComponent>("max_steps"))
  , _dt(Get<Real>("dt"))
{
}

void LJInitCondition::Execute() 
{
  UpdateField();
  MeshFreeCondition::Execute();
  auto step = 0;
  while (step < _max_steps)
  {
    DoInit();
    step++;
  }
  //console::Info("Init System End!");
}

void LJInitCondition::UpdateField()
{
  MeshFreeCondition::UpdateField();
}

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
  auto pts_type = _system.GetFieldAsArrayHandle<Id>(field::pts_type);
  auto atom_id = _system.GetFieldAsArrayHandle<Id>(field::atom_id);
  try
  {
    ContForceFunction force_function;

    ContTopology topology;
    SetTopology(topology);

    ContPointLocator locator;
    SetLocator(locator);
    auto cut_off = _system.GetParameter<Real>(PARA_CUTOFF);
    SystemWorklet::Class2LJForceWithPeriodicBC(cut_off, atom_id, locator, topology, force_function, _LJforce);
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
    auto velocity = _system.GetFieldAsArrayHandle<Vec3f>(field::velocity);
    auto mass = _system.GetFieldAsArrayHandle<Real>(field::mass);

    auto unit_factor = _system.GetParameter<UnitFactor>(PARA_UNIT_FACTOR);
    SystemWorklet::UpdateVelocity(_dt, unit_factor._fmt2v, _LJforce, mass, velocity);
  }
  catch (const std::exception& e)
  {
    std::cout << e.what() << std::endl;
  }
}

void LJInitCondition::UpdatePosition() 
{
  auto velocity = _system.GetFieldAsArrayHandle<Vec3f>(field::velocity);

  ContPointLocator locator;
  SetLocator(locator);
  SystemWorklet::UpdatePosition(_dt, velocity, locator, _position);
}
