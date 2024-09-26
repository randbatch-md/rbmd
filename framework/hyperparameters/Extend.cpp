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
//  The post-processing data produced by VASPKIT may be used in any publications 
//  provided that its use is explicitly acknowledged. A suitable reference for VASPKIT is:
//  [1] Gao W, Zhao T, Guo Y, et al.RBMD: A molecular dynamics package enabling to simulate 
//  10 million all - atom particles in a single graphics processing unit[J].arXiv preprint arXiv : 2407.09315, 2024.
// 
//  Contact Email : [support_wz@sciai.com.cn]
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
