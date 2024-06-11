#include "MeshFreeSystem.h"
#include "FieldName.h"

MeshFreeSystem::MeshFreeSystem(const Configuration& cfg)
  : System(cfg)
  , _unit(Get<std::string>("unit"))
{
  InitPara();
  _para.SetParameter(PARA_UNIT, _unit);
}

void MeshFreeSystem::Init()
{
  System::Init();
  _position = _para.GetFieldAsArrayHandle<Vec3f>(field::position);
  InitPointLocator();
}

void MeshFreeSystem::InitField()
{
  _para.AddField(field::position, ArrayHandle<Vec3f>{});
}

void MeshFreeSystem::Evolve()
{
  _timer.Start();
  PreSolve();
  Solve();
  PostSolve();
}

void MeshFreeSystem::InitPara()
{
  //cut off
  auto& parser = _app.GetParser();
  if (nullptr != parser)
  {
    auto& locator_node = parser->GetJsonNode("locator");
    _para.SetParameter(PARA_CUTOFF, static_cast<Real>(locator_node["cut_off"].asFloat()));
    _para.SetParameter(PARA_R_CORE, static_cast<Real>(locator_node["rs"].asFloat()));
    _para.SetParameter(PARA_NEIGHBOR_SAMPLE_NUM, static_cast<Real>(locator_node["random_num"].asFloat()));  }
  else
  {
    console::Error("Parser is NULL!");
  }


  //unit factor
  if (_unit == "REAL")
  {
    _unit_factor._kB = 1.9872067 * vtkm::Pow(10.0, -3);
    _unit_factor._fmt2v = 4.186 * vtkm::Pow(10.0, -4);
    _unit_factor._mvv2e = 1.0 / (4.186 * vtkm::Pow(10.0, -4));
    _unit_factor._qqr2e = 332.06371;
  }
  else if (_unit == "LJ")
  {
    _unit_factor._kB = 1.0;
    _unit_factor._fmt2v = 1.0;
    _unit_factor._mvv2e = 1.0;
    _unit_factor._qqr2e = 1.0;
  }
  else if (_unit == "METAL")
  {
    _unit_factor._kB = 8.617343e-5;
    _unit_factor._fmt2v = 1.0 / 1.0364269e-4;
    _unit_factor._mvv2e = 1.0364269e-4;
    _unit_factor._qqr2e = 14.399645;
  }
  _para.SetParameter(PARA_UNIT_FACTOR, _unit_factor);
}

void MeshFreeSystem::InitPointLocator()
{
  vtkm::Vec<vtkm::Range, 3> range = _para.GetParameter<vtkm::Vec<vtkm::Range, 3>>(PARA_RANGE);
  vtkm::Vec3f left_bottom{ {
                             static_cast<vtkm::FloatDefault>(range[0].Min),
                           },
                           {
                             static_cast<vtkm::FloatDefault>(range[1].Min),
                           },
                           {
                             static_cast<vtkm::FloatDefault>(range[2].Min),
                           } };
  vtkm::Vec3f right_top{ {
                           static_cast<vtkm::FloatDefault>(range[0].Max),
                         },
                         {
                           static_cast<vtkm::FloatDefault>(range[1].Max),
                         },
                         {
                           static_cast<vtkm::FloatDefault>(range[2].Max),
                         } };
  _locator.SetRange(left_bottom, right_top);

  _locator.SetCutOff(_para.GetParameter<Real>(PARA_CUTOFF));

  _locator.SetRs(_para.GetParameter<Real>(PARA_R_CORE));

  _locator.SetPosition(_position);
}
