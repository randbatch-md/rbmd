#include "Extend.h"
#include "FieldName.h"
Extend::Extend(const Configuration& cfg)
  : HyperParameters(cfg)
{
}

void Extend::Execute()
{
  SetParameters();
}

void Extend::SetParameters()
{
  _para.SetParameter(PARA_FIX_SHAKE, Get<std::string>("fix_shake"));
}
