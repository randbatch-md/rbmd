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
