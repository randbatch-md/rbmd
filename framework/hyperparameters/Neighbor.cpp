//==================================================================================
//  RBMD 1.0 is developed for random batch molecular dynamics calculation.
//
//  Copyright(C) 2024 SHANGHAI JIAOTONG UNIVERSITY CHONGQING RESEARCH INSTITUTE
//
//  This program is free software : you can redistribute it and /or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with this program.If not, see < https://www.gnu.org/licenses/>.
//
//  Contact Email : [your - email@example.com]
//==================================================================================

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
    _para.SetParameter(PARA_R_CORE, Get<Real>("cut_off")); 
  }
}
