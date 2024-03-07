#pragma once
#include "System.h"
#include "DataObject.h"
#include "locator/ContPointLocator.h"
#include "forceFunction/ContForceFunction.h"
#include "topology/ContTopology.h"
#include "UnitFactor.h"
#include <vtkm/cont/ArrayCopy.h>
#include "staticTable/ContStaticTable.h"

class MeshFreeSystem : public System
{
public:
  MeshFreeSystem(const Configuration& cfg);
  virtual ~MeshFreeSystem() = default;

  void Init() override;
  void Evolve() override;

protected:
  virtual void PreSolve(){};
  virtual void Solve(){};
  virtual void PostSolve(){};
  void InitField() override;

private:
  void InitPara();
  void InitPointLocator();

protected:
  std::string _unit;
  UnitFactor _unit_factor;
  ArrayHandle<Vec3f> _position;
  ContPointLocator _locator;

  ContStaticTable _static_table;
};