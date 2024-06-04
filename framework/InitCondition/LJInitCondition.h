#pragma once
#include "MeshFreeCondition.h"

class LJInitCondition : public MeshFreeCondition
{
public:
  LJInitCondition(const Configuration& cfg);
  virtual ~LJInitCondition()=default;
  
  void Execute() override;
  void UpdateField() override;
  //void InitParameter() override;
  void InitField() override;
  void SetParameters() override;

private:
  void DoInit();
  void ComputeForce();
  void UpdateVelocity();
  void UpdatePosition();

private:
  ArrayHandle<Vec3f> _LJforce;
  IdComponent _max_steps;
  Real _dt;
};
