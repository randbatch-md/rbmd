#include "Executioner.h"
#include "MoleculesFileOutput.h"
#include "Application.h"
#include "LennardJonesSystem.h"
#include <vtkm/cont/Algorithm.h>
#include "output/worklet/OutPutWorklet.h"
#include "FieldName.h"

//RegisterObject(MoleculesFileOutput);

MoleculesFileOutput::MoleculesFileOutput(const Configuration& cfg)
  : FileOutput(cfg)
  , _RDF_file(_file_base + ".csv")
  , _radius(Get<Real>("radius"))
  , _min_radius(Get<Real>("min_radius"))
  , _dr(Get<Real>("dr"))
  , _statistics_rdf_steps(Get<Real>("statistics_rdf_steps"))
  , _comput_RDF(Get<bool>("compute"))
  , _outpute_file(Get<bool>("output_file"))
  , _center_type(Get<vtkm::IdComponent>("center_type"))
  , _target_type(Get<vtkm::IdComponent>("target_type"))
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
  _system.SetParameter(PARA_CENTER_TYPE, Get<vtkm::IdComponent>("center_type"));
  _system.SetParameter(PARA_TARGET_TYPE, Get<vtkm::IdComponent>("target_type"));
}

void MoleculesFileOutput::Init() 
{
  if (_comput_RDF)
  {
    _rdf_rho = _system.GetParameter<Real>(PARA_RDF_RHO);
    FileOutput::Init();
  }
}

void MoleculesFileOutput::Execute()
{
  if (_comput_RDF)
  {
    _rdfTimer.Start();
    ComputeRDF();
    if (ShouldOutput() && _outpute_file)
    {
      try
      {
        _RDF_file << "RADIUS,RDF" << std::endl;
        auto r_num = _vRadius.size();
        auto rdf_Protol = _rdf.ReadPortal();
        for (auto i = 0; i < r_num; ++i)
        {
          auto r = _vRadius[i];
          auto RDF = rdf_Protol.Get(i);
          Real rdf = (RDF / (_step_upper - _step_lower - 1));
          if (std::isnan(rdf))
          {
            _RDF_file << r << "," << std::endl;
          }
          else
          {
            _RDF_file << r << "," << rdf << std::endl;
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
}

bool MoleculesFileOutput::ShouldOutput()
{
  if (_executioner.CurrentStep() == _executioner.NumStep() - 1)
  {
    return true;
  }
  return false;
}

void MoleculesFileOutput::ComputeRDF()
{
  if (_executioner.CurrentStep() >= _step_lower && _executioner.CurrentStep() < _step_upper)
  {
    auto center_position = _system.GetFieldAsArrayHandle<Vec3f>(field::center_position);
    auto target_position = _system.GetFieldAsArrayHandle<Vec3f>(field::target_position);
    auto atom_id_center = _system.GetFieldAsArrayHandle<Id>(field::atom_id_center);
    auto atom_id_target = _system.GetFieldAsArrayHandle<Id>(field::atom_id_target);
    auto molecule_id = _system.GetFieldAsArrayHandle<Id>(field::molecule_id);
    auto position = _system.GetFieldAsArrayHandle<Vec3f>(field::position);

    auto num_center_pos = center_position.GetNumberOfValues();
    ContPointLocator locator;
    SetLocator(locator);
    locator.SetPosition(target_position, atom_id_target);

    auto radius = vtkm::cont::make_ArrayHandle(_vRadius);
    OutPut::ComputeRDF(num_center_pos,
                       _rdf_rho,
                       radius,
                       center_position,
                       atom_id_center,
                       molecule_id,
                       locator,
                       _rdf);
  }
}