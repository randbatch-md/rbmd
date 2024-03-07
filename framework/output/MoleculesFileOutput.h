#pragma once

#include "FileOutput.h"
#include <vtkm/cont/ArrayHandleGroupVecVariable.h>

class MoleculesFileOutput : public FileOutput
{
public:
  MoleculesFileOutput(const Configuration& cfg);
  virtual ~MoleculesFileOutput(){};

  void Init() override;
  void Execute() override;

protected:
  void ComputeRDF();
  bool ShouldOutput() override;

private:
  bool _binary = false;
  bool _comput_RDF;
  bool _outpute_file;

  vtkm::IdComponent _step_lower;
  vtkm::IdComponent _step_upper;
  const Real _radius;
  const Real _min_radius;
  const Real _dr;
  const vtkm::IdComponent _center_type;
  const vtkm::IdComponent _target_type;
  vtkm::IdComponent _statistics_rdf_steps;

  std::ofstream _RDF_file;
  vtkm::cont::Timer _rdfTimer;

  std::vector<Real> _vRadius;
  vtkm::cont::ArrayHandle<Real> _rdf;
  vtkm::cont::ArrayHandle<Real> rdf;

  Real _rdf_rho;
};