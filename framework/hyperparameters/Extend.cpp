#include "Extend.h"
#include "FieldName.h"
Extend::Extend(const Configuration& cfg)
  : HyperParameters(cfg)
{
}

void Extend::Execute()
{
  // 可以存在 fix_shake 和 special_bonds 或者不存在
  _para.SetParameter(PARA_FIX_SHAKE, GetJudge<std::string>("fix_shake"));
  auto a = _para.GetParameter<std::string>(PARA_FIX_SHAKE);
  _para.SetParameter(PARA_SPECIAL_BONDS, GetVectorValueJudge<int>("special_bonds"));
}
