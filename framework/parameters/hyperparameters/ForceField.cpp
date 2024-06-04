#include "ForceField.h"
#include "FieldName.h"
ForceField::ForceField(const Configuration& cfg)
  : HyperParameters(cfg)
{
}

void ForceField::Execute() 
{
  SetPara();
}

void ForceField::SetPara()
{
  _para.SetParameter(PARA_FORCE_FIELD_TYPE, Get<std::string>("type"));
  // 注意：这一块目前是在算EAM使用，有的案例比如H2O根本没有
  // 所以这一块 是否需要提出来在这里初始化、或者进行判断 看是否这只了是能文件
  if (Get<std::string>("type") == "EAM")
  {
    _para.SetParameter(PARA_POTENTIAL_FILE, Get<std::string>("potential_file")); 
  }
}
