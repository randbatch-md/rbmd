#include "FieldName.h"
#include "Neighbor.h"
Neighbor::Neighbor(const Configuration& cfg)
  : HyperParameters(cfg)
{
}

void Neighbor::Execute()
{
  _para.SetParameter(PARA_NEIGHBOR_TYPE, Get<std::string>("type"));
  _para.SetParameter(PARA_CUTOFF, Get<Real>("cut_off"));

  // 注意：这里如果是 verletlist方法 locator中 locator.SetRs(GetParameter<Real>(PARA_R_CORE));要修改合理
  if (Get<std::string>("type") == "RBL")
  {
    _para.SetParameter(PARA_R_CORE, Get<Real>("r_core"));
    _para.SetParameter(PARA_NEIGHBOR_SAMPLE_NUM, Get<Real>("neighbor_sample_num"));
  }
  else
  {
    _para.SetParameter(PARA_R_CORE, Get<Real>("cut_off")); //Get<Real>("cut_off"));
  }
}
