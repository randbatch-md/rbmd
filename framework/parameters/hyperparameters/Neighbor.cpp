#include "FieldName.h"
#include "Neighbor.h"
Neighbor::Neighbor(const Configuration& cfg)
  : HyperParameters(cfg)
{
}

void Neighbor::Execute()
{
  SetPara();
}

void Neighbor::SetPara()
{
  _para.SetParameter(PARA_NEIGHBOR_TYPE, Get<std::string>("type"));
  _para.SetParameter(PARA_CUTOFF, Get<Real>("cut_off"));
  auto cut_off = _para.GetParameter<Real>(PARA_CUTOFF);
  auto xLength = _para.GetParameter<Real>(PARA_VLENGTH);
  auto yLength = xLength;
  auto zLength = xLength;
  auto bin_number = Id3{static_cast<int>(xLength / cut_off),static_cast<int>(yLength / cut_off),static_cast<int>(zLength / cut_off),};
  _para.SetParameter(PARA_BIN_NUMBER, bin_number);

  // 注意：这里如果是 verletlist方法 locator中 locator.SetRs(GetParameter<Real>(PARA_R_CORE));要修改合理
  if (Get<std::string>("type") == "RBL")
  {
    _para.SetParameter(PARA_R_CORE, Get<Real>("r_core"));
    _para.SetParameter(PARA_NEIGHBOR_SAMPLE_NUM, Get<Real>("neighbor_sample_num"));
  }
}
