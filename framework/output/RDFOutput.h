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
  const std::vector<std::vector<int>> _atoms_pair;
  vtkm::IdComponent _statistics_rdf_steps;

  std::ofstream _RDF_file;
  vtkm::cont::Timer _rdfTimer;

  std::vector<Real> _vRadius;
  std::vector<vtkm::cont::ArrayHandle<Real>> _rdf;
  vtkm::cont::ArrayHandle<Real> rdf;

  Real _rdf_rho;
};