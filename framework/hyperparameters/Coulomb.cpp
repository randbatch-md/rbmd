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
#include "Coulomb.h"
Coulomb::Coulomb(const Configuration& cfg)
  : HyperParameters(cfg)
{
}

void Coulomb::Execute()
{
  _para.SetParameter(PARA_COULOMB_TYPE, Get<std::string>("type"));
  // 注：这里如果是 verletlist方法 locator中 locator.SetRs(GetParameter<Real>(PARA_R_CORE));要修改合理
  _para.SetParameter(PARA_COULOMB_SAMPLE_NUM, Get<IdComponent>("coulomb_sample_num"));
  _para.SetParameter(PARA_ALPHA, Get<Real>("alpha"));
  _para.SetParameter(PARA_KMAX, Get<IdComponent>("kmax"));
}
