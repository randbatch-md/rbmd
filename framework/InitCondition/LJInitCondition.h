#pragma once
#include "MeshFreeCondition.h"

class LJInitCondition : public MeshFreeCondition
{
public:
  LJInitCondition(const Configuration& cfg);
  virtual ~LJInitCondition()=default;
  
  void Execute() override;
  void UpdateField() override;
  void AddMoleculeInfo();
  void InitField() override;
  void SetParameters() override;

private:
  ArrayHandle<Vec3f> _LJforce;
  IdComponent _max_steps;
  Real _dt;
};
