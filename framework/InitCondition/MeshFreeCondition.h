#pragma once
#include "InitCondition.h"
#include "forceFunction/ContForceFunction.h"
#include "locator/ContPointLocator.h"
#include "topology/ContTopology.h"
#include "DataObject.h"
#include "FieldName.h"

class MeshFreeCondition: public InitCondition
{
public:
  MeshFreeCondition(const Configuration& cfg);
  ~MeshFreeCondition();
  void Execute() override;
  void UpdateField() override;
  void InitField() override;
  void SetParameters() override;

private:
  void InitPosition();
  void InitMassAndVelocity();
  void InitPositionFlag();
  void InitId();
  void InitParameters();

protected:
  ArrayHandle<Vec3f> _position;
  std::vector<int> _dims;
  std::vector<Real> _x_range;
  std::vector<Real> _y_range;
  std::vector<Real> _z_range;
};