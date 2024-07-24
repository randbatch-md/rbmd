#include "FieldName.h"
#include "Coulomb.h"
Coulomb::Coulomb(const Configuration& cfg)
  : HyperParameters(cfg)
{
}

void Coulomb::Execute()
{
  _para.SetParameter(PARA_COULOMB_TYPE, Get<std::string>("type"));
  // 注意：这里如果事 verletlist方法 locator中 locator.SetRs(GetParameter<Real>(PARA_R_CORE));要修改合理
  //if (Get<std::string>("type") == "RBE")
  //{
  _para.SetParameter(PARA_COULOMB_SAMPLE_NUM, Get<IdComponent>("coulomb_sample_num"));
  _para.SetParameter(PARA_ALPHA, Get<Real>("alpha"));
  _para.SetParameter(PARA_KMAX, Get<IdComponent>("kmax"));
  //}
}
