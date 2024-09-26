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
}
