//==================================================================================
//  RBMD 2.1.0 is developed for random batch molecular dynamics calculation.
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

#pragma once
#include <any>
#include <map>
#include <vtkm/cont/ArrayHandle.h>
#include <vtkm/cont/Field.h>

class Para
{
public:
  Para(){};

  ~Para(){};

  template<typename T>
  void AddField(const std::string& name, const vtkm::cont::ArrayHandle<T>& arrayHandle)
  {
    _field.insert(std::make_pair(
      name, vtkm::cont::Field(name, vtkm::cont::Field::Association::Any, arrayHandle)));
  }

  template<typename T>
  vtkm::cont::ArrayHandle<T> GetFieldAsArrayHandle(const std::string& name)
  {
    try
    {
      return _field[name].GetData().AsArrayHandle<vtkm::cont::ArrayHandle<T>>();
    }
    catch (const std::exception& e)
    {
      console::Error(e.what());
    }
  }

  template<typename T>
  void SetParameter(const std::string& name, const T& field)
  {
    _parameters[name] = field;
  }

  template<typename T>
  T GetParameter(const std::string& name)
  {
    try
    {
      return std::any_cast<T>(_parameters.at(name));
    }
    catch (const std::exception& e)
    {
      console::Error(e.what());
    }
  }

  bool HaveParameter(const std::string& name) const
  {
    return _parameters.find(name) != _parameters.end();
  }

protected:
  std::map<std::string, vtkm::cont::Field> _field;
  std::map<std::string, std::any> _parameters;
};