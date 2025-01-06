//==================================================================================
//  RBMD 2.2.0 is developed for random batch molecular dynamics calculation.
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

#include "ForceField.h"
#include "FieldName.h"
ForceField::ForceField(const Configuration& cfg)
  : HyperParameters(cfg)
{
}

void ForceField::Execute() 
{
  _para.SetParameter(PARA_FORCE_FIELD_TYPE, Get<std::string>("type"));
  _para.SetParameter(PARA_FILE_TYPE, (std::string) "NOT_EAM");
  // 注意：这一块目前是在算EAM使用，有的案例比如H2O根本没有
  // 所以这一块 是否需要提出来在这里初始化、或者进行判断 看是否这只了是能文件
  if (Get<std::string>("type") == "EAM")
  {
    _para.SetParameter(PARA_POTENTIAL_FILE, Get<std::string>("potential_file"));
    _para.SetParameter(PARA_FILE_TYPE, (std::string) "EAM");
  }
  else if ("CVFF" == Get<std::string>("type"))
  {
    _para.SetParameter(PARA_BOND_TYPE, Get<std::string>("bond_type"));
    _para.SetParameter(PARA_ANGLE_TYPE, Get<std::string>("angle_type"));
    _para.SetParameter(PARA_DIHEDRAL_TYPE, GetJudge<std::string>("dihedral_type"));
    _para.SetParameter(PARA_IMPROPER_TYPE, GetJudge<std::string>("improper_type"));
  }
}
