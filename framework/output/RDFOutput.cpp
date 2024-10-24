﻿//==================================================================================
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

#include "Executioner.h"
#include "RDFOutput.h"
#include "Application.h"
#include <vtkm/cont/Algorithm.h>
#include "output/worklet/OutPutWorklet.h"
#include "FieldName.h"

RDFOutput::RDFOutput(const Configuration& cfg)
  : FileOutput(cfg)
  , _RDF_file("rdf.rbmd")
  , _radius(Get<Real>("radius"))
  , _min_radius(0)
  , _dr(Get<Real>("dr"))
  , _statistics_rdf_steps(Get<Real>("statistics_rdf_steps"))
{
  _step_lower = _executioner.NumStep() - _statistics_rdf_steps;
  _step_upper = _executioner.NumStep() + 1;

 auto radius_num = (_radius - _min_radius) / _dr;

  _vRadius.resize(radius_num);
  std::generate(_vRadius.begin(),
                _vRadius.end(),
                [this](void) -> Real
                {
                  static Real min = _min_radius;
                  min += _dr;
                  return min;
                });
  auto init_way = _para.GetParameter<std::string>(PARA_INIT_WAY);
  if (init_way != "inbuild")
  {
    _atoms_pair = GetVectorOfVectorsValue<int>("atoms_pair");
    for (int i = 0; i < _atoms_pair.size(); i++)
    {
      vtkm::cont::ArrayHandle<Real> rdf_pair;
      rdf_pair.AllocateAndFill(_vRadius.size(), 0);
      _rdf.push_back(rdf_pair);
    }
    std::vector<int> atoms_pair_type;
    int i, type;
    int size = _atoms_pair.size() * 2 - 1;
    while (size >= 0)
    {
      i = 0;
      type = _atoms_pair[size / 2][size % 2];
      if (atoms_pair_type.size() == 0)
      {
        atoms_pair_type.push_back(type);
      }
      else
      {
        while (i < atoms_pair_type.size())
        {
          if (type != atoms_pair_type[i] && i == atoms_pair_type.size() - 1)
          {
            atoms_pair_type.push_back(type);
          }
          else if (type == atoms_pair_type[i])
          {
            break;
          }
          else
            i++;
        }
      }
      size--;
    }
    _para.SetParameter(PARA_ATOMS_PAIR_TYPE, atoms_pair_type);
  }
  else
  {
    vtkm::cont::ArrayHandle<Real> rdf_pair;
    rdf_pair.AllocateAndFill(_vRadius.size(), 0);
    _rdf.push_back(rdf_pair);
  }
}

void RDFOutput::Init() 
{
  _rdf_rho = _para.GetParameter<Real>(PARA_RDF_RHO);
  FileOutput::Init();
}

void RDFOutput::Execute()
{
  _rdfTimer.Start();
  ComputeRDF();
  if (ShouldOutput())
  {
    try
    {
      _RDF_file << "Radius ";
      auto init_way = _para.GetParameter<std::string>(PARA_INIT_WAY);
      if (init_way == "inbuild")
      {
        _RDF_file << "RDF" << std::endl;
      }
      else
      {
        for (int i = 0; i < _atoms_pair.size(); i++)
        {
          _RDF_file << "RDF:" << _atoms_pair[i][0] << "-" << _atoms_pair[i][1] << " ";
        }
        _RDF_file << std::endl;
      }
      
      auto r_num = _vRadius.size();
      
      for (auto i = 0; i < r_num; ++i)
      {
        auto r = _vRadius[i];
        _RDF_file << r ; 
        for (int j = 0; j < _rdf.size(); j++)
        {
          auto rdf_Protol = _rdf[j].ReadPortal();
          auto RDF = rdf_Protol.Get(i);
          Real rdf = (RDF / (_step_upper - _step_lower - 1));
          if (std::isnan(rdf))
          {
            _RDF_file << " ";
          }
          else
          {
            _RDF_file << " " << rdf;
          }
        }
        _RDF_file << std::endl;
      }
      _RDF_file.close();
    }
    catch (const std::exception& e)
    {
      _RDF_file.close();
      console::Error(e.what());
    }
  }
}

bool RDFOutput::ShouldOutput()
{
  if (_executioner.CurrentStep() == _executioner.NumStep() - 1)
  {
    return true;
  }
  return false;
}

void RDFOutput::ComputeRDF()
{
  if (_executioner.CurrentStep() >= _step_lower && _executioner.CurrentStep() < _step_upper)
  {
    auto molecule_id = _para.GetFieldAsArrayHandle<Id>(field::molecule_id);
    auto position = _para.GetFieldAsArrayHandle<Vec3f>(field::position);
    auto init_way = _para.GetParameter<std::string>(PARA_INIT_WAY);
    if (init_way == "inbuild")
    {
      auto center_position = _para.GetFieldAsArrayHandle<Vec3f>(field::center_position);
      auto target_position = _para.GetFieldAsArrayHandle<Vec3f>(field::target_position);
      auto num_center_pos = center_position.GetNumberOfValues();
      ContPointLocator locator;
      SetLocator(locator);

      locator.SetPosition(target_position);
      auto radius = vtkm::cont::make_ArrayHandle(_vRadius);
      OutPut::atoms::ComputeRDF(
        num_center_pos, _rdf_rho, radius, center_position, target_position, locator, _rdf[0]);
    }
    else
    {
      std::map<Id, ArrayHandle<Vec3f>> atom_pair_position;
      std::map<Id, std::vector<Id>> atom_pair_id;
      auto atom_pair_type = _para.GetParameter<std::vector<int>>(PARA_ATOMS_PAIR_TYPE);
      auto rdf_id = _para.GetFieldAsArrayHandle<Id>(field::atom_pair_id);
      auto rdf_position = _para.GetFieldAsArrayHandle<Vec3f>(field::atom_pair_position);
      auto atoms_pair_type_offsets =
        _para.GetFieldAsArrayHandle<Id>(field::atoms_pair_type_offsets);
      auto id_vec_arrayhandle =
        vtkm::cont::make_ArrayHandleGroupVecVariable(rdf_id, atoms_pair_type_offsets);
      auto position_vec_arrayhandle =
        vtkm::cont::make_ArrayHandleGroupVecVariable(rdf_position, atoms_pair_type_offsets);
      auto type = id_vec_arrayhandle.GetNumberOfValues();
      auto id_read_protol = id_vec_arrayhandle.GetPortalControl();
      auto position_read_protol = position_vec_arrayhandle.GetPortalControl();

      for (int i = 0; i < type; ++i)
      {
        auto id_vec = id_read_protol.Get(i);
        auto position_Vec = position_read_protol.Get(i);
        std::vector<vtkm::Id> id_data(id_vec.GetNumberOfComponents());
        std::vector<Vec3f> position_data(position_Vec.GetNumberOfComponents());
        for (vtkm::IdComponent j = 0; j < id_vec.GetNumberOfComponents(); ++j)
        {
          id_data[j] = id_vec[j];
          position_data[j] = position_Vec[j].Get();
        }
        //vtkm::cont::ArrayHandle<Id> id_array;
        //vtkm::cont::ArrayCopy(vtkm::cont::make_ArrayHandle(id_data), id_array);
        vtkm::cont::ArrayHandle<Vec3f> atom_type_position;
        vtkm::cont::ArrayCopy(vtkm::cont::make_ArrayHandle(position_data), atom_type_position);
        atom_pair_position[atom_pair_type[i]] = atom_type_position;
        atom_pair_id[atom_pair_type[i]] = id_data;
      }

      for (int i = 0; i < _atoms_pair.size(); i++)
      {
        auto center_position = atom_pair_position[_atoms_pair[i][0]];
        auto target_position = atom_pair_position[_atoms_pair[i][1]];

        vtkm::cont::ArrayHandle<Id> atom_id_center;
        vtkm::cont::ArrayCopy(vtkm::cont::make_ArrayHandle(atom_pair_id[_atoms_pair[i][0]]), atom_id_center);
        
        vtkm::cont::ArrayHandle<Id> atom_id_target;
        vtkm::cont::ArrayCopy(vtkm::cont::make_ArrayHandle(atom_pair_id[_atoms_pair[i][1]]), atom_id_target);

        auto num_center_pos = center_position.GetNumberOfValues();
        ContPointLocator locator;
        SetLocator(locator);
        locator.SetPosition(target_position, atom_id_target);

        _rdf_rho = target_position.GetNumberOfValues() / _para.GetParameter<Real>(PARA_VOLUME);
        auto radius = vtkm::cont::make_ArrayHandle(_vRadius);
        OutPut::ComputeRDF(num_center_pos,
                           _rdf_rho,
                           radius,
                           center_position,
                           position,
                           atom_id_center,
                           molecule_id,
                           locator,
                           _rdf[i]);
      }
    }
    }

}