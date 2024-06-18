#include "Executioner.h"
#include "RDFOutput.h"
#include "Application.h"
#include "LennardJonesSystem.h"
#include <vtkm/cont/Algorithm.h>
#include "output/worklet/OutPutWorklet.h"
#include "FieldName.h"

//RegisterObject(RDFOutput);

RDFOutput::RDFOutput(const Configuration& cfg)
  : FileOutput(cfg)
  , _RDF_file("rdf.rbmd")
  , _radius(Get<Real>("radius"))
  , _min_radius(0)
  , _dr(Get<Real>("dr"))
  , _statistics_rdf_steps(Get<Real>("statistics_rdf_steps"))
  , _atoms_pair(GetVectorOfVectorsValue<int>("atoms_pair"))
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

  _rdf.AllocateAndFill(radius_num, 0);
  _para.SetParameter(PARA_CENTER_TYPE, _atoms_pair[0][0]);
  _para.SetParameter(PARA_TARGET_TYPE, _atoms_pair[0][1]);
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
      _RDF_file << "Radius RDF" << std::endl;
      auto r_num = _vRadius.size();
      auto rdf_Protol = _rdf.ReadPortal();
      for (auto i = 0; i < r_num; ++i)
      {
        auto r = _vRadius[i];
        auto RDF = rdf_Protol.Get(i);
        Real rdf = (RDF / (_step_upper - _step_lower - 1));
        if (std::isnan(rdf))
        {
          _RDF_file << r << " " << std::endl;
        }
        else
        {
          _RDF_file << r << " " << rdf << std::endl;
        }
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
    auto center_position = _para.GetFieldAsArrayHandle<Vec3f>(field::center_position);
    auto target_position = _para.GetFieldAsArrayHandle<Vec3f>(field::target_position);
    auto atom_id_center = _para.GetFieldAsArrayHandle<Id>(field::atom_id_center);
    auto atom_id_target = _para.GetFieldAsArrayHandle<Id>(field::atom_id_target);
    auto molecule_id = _para.GetFieldAsArrayHandle<Id>(field::molecule_id);
    auto position = _para.GetFieldAsArrayHandle<Vec3f>(field::position);

    auto num_center_pos = center_position.GetNumberOfValues();
    ContPointLocator locator;
    SetLocator(locator);
    locator.SetPosition(target_position, atom_id_target);

    auto radius = vtkm::cont::make_ArrayHandle(_vRadius);
    OutPut::ComputeRDF(num_center_pos,
                   _rdf_rho,
                   radius,
                   center_position,
                   position,
                   atom_id_center,
                   molecule_id,
                   locator,
                   _rdf);
  }
}