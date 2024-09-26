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

#pragma once

#include "FileOutput.h"
#include <vtkm/cont/ArrayHandleGroupVecVariable.h>

class RDFOutput : public FileOutput
{
public:
  RDFOutput(const Configuration& cfg);
  virtual ~RDFOutput(){};

  void Init() override;
  void Execute() override;

protected:
  void ComputeRDF();
  bool ShouldOutput() override;

private:
  vtkm::IdComponent _step_lower;
  vtkm::IdComponent _step_upper;
  const Real _radius;
  const Real _min_radius;
  const Real _dr;
  std::vector<std::vector<int>> _atoms_pair;
  vtkm::IdComponent _statistics_rdf_steps;

  std::ofstream _RDF_file;
  vtkm::cont::Timer _rdfTimer;

  std::vector<Real> _vRadius;
  std::vector<vtkm::cont::ArrayHandle<Real>> _rdf;
  vtkm::cont::ArrayHandle<Real> rdf;

  Real _rdf_rho;
};