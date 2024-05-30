#include "MeshFreeSystem.h"
#include "FieldName.h"

MeshFreeSystem::MeshFreeSystem(const Configuration& cfg)
  : System(cfg)
  , _unit(Get<std::string>("unit"))
{
  InitPara();
  SetParameter(PARA_UNIT, _unit);
}

void MeshFreeSystem::Init()
{
  System::Init();
  _position = GetFieldAsArrayHandle<Vec3f>(field::position);
  InitPointLocator();
}

void MeshFreeSystem::InitField()
{
  AddField(field::position, ArrayHandle<Vec3f>{});
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
    SetParameter(PARA_CUTOFF, static_cast<Real>(locator_node["cut_off"].asFloat()));
    SetParameter(PARA_RS, static_cast<Real>(locator_node["rs"].asFloat()));
    SetParameter(PARA_RANDOM_NUM, static_cast<Real>(locator_node["random_num"].asFloat()));  }
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
  SetParameter(PARA_UNIT_FACTOR, _unit_factor);
}

void MeshFreeSystem::InitPointLocator()
{
  vtkm::Vec<vtkm::Range, 3> range = GetParameter<vtkm::Vec<vtkm::Range, 3>>(PARA_RANGE);
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

  _locator.SetCutOff(GetParameter<Real>(PARA_CUTOFF));

  _locator.SetRs(GetParameter<Real>(PARA_RS));

  _locator.SetPosition(_position);
}
