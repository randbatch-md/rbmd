
#pragma once

#include "Object.h"
#include "ContPointLocator.h"
#include "System.h"
#include "forceFunction/ContForceFunction.h"
#include "FieldName.h"
#include "topology/ContTopology.h"

class Application;
class Output : public Object
{
public:
  Output(const Configuration& cfg);

  virtual ~Output(){};

  virtual void Init();
  virtual void Execute() = 0;

protected:
  void SetTopology(ContTopology& topology)
  {
    auto pts_type = _system.GetFieldAsArrayHandle<Id>(field::pts_type);
    auto epsilon = _system.GetFieldAsArrayHandle<Real>(field::epsilon);
    auto sigma = _system.GetFieldAsArrayHandle<Real>(field::sigma);
    auto charge = _system.GetFieldAsArrayHandle<Real>(field::charge );
    auto molecule_id = _system.GetFieldAsArrayHandle<Id>(field::molecule_id);

    topology.SetAtomsType(pts_type);
    topology.SetEpsAndSigma(epsilon, sigma);
    topology.SetCharge(charge);
    topology.SetMolecularId(molecule_id);
  }

  void SetLocator(ContPointLocator& locator)
  {
      vtkm::Vec<vtkm::Range, 3> range = _system.GetParameter<vtkm::Vec<vtkm::Range, 3>>(PARA_RANGE);
      vtkm::Vec3f left_bottom{ {static_cast<vtkm::FloatDefault>(range[0].Min),},
                               {static_cast<vtkm::FloatDefault>(range[1].Min),},
                               {static_cast<vtkm::FloatDefault>(range[2].Min),} };
      vtkm::Vec3f right_top{ {static_cast<vtkm::FloatDefault>(range[0].Max),},
                             {static_cast<vtkm::FloatDefault>(range[1].Max),},
                             {static_cast<vtkm::FloatDefault>(range[2].Max),} };
      
      locator.SetRange(left_bottom, right_top);
      
      locator.SetCutOff(_system.GetParameter<Real>(PARA_CUTOFF));
      
      locator.SetRs(_system.GetParameter<Real>(PARA_RS));

      locator.SetPosition(_system.GetFieldAsArrayHandle<Vec3f>(field::position));
  }

  void SetForceFunction(ContForceFunction& force_function)
  {
      auto rhomax = _system.GetParameter<Real>(EAM_PARA_RHOMAX);
      auto nrho = _system.GetParameter<Id>(EAM_PARA_NRHO);
      auto drho = _system.GetParameter<Real>(EAM_PARA_DRHO);
      auto nr = _system.GetParameter<Id>(EAM_PARA_NR);
      auto dr = _system.GetParameter<Real>(EAM_PARA_DR);
      force_function.SetEAMParameters(rhomax, nrho, drho, nr, dr);
  }

protected:
  Application& _app;
  System& _system;
  std::shared_ptr<Executioner>& _executioner;
};
