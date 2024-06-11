
#pragma once

#include "Object.h"
#include "ContPointLocator.h"
//#include "System.h"
#include "forceFunction/ContForceFunction.h"
#include "FieldName.h"
#include "topology/ContTopology.h"
#include "Para.h"

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
    auto pts_type = _para.GetFieldAsArrayHandle<Id>(field::pts_type);
    auto epsilon = _para.GetFieldAsArrayHandle<Real>(field::epsilon);
    auto sigma = _para.GetFieldAsArrayHandle<Real>(field::sigma);
    auto charge = _para.GetFieldAsArrayHandle<Real>(field::charge );
    auto molecule_id = _para.GetFieldAsArrayHandle<Id>(field::molecule_id);

    topology.SetAtomsType(pts_type);
    topology.SetEpsAndSigma(epsilon, sigma);
    topology.SetCharge(charge);
    topology.SetMolecularId(molecule_id);
  }

  void SetLocator(ContPointLocator& locator)
  {
      vtkm::Vec<vtkm::Range, 3> range = _para.GetParameter<vtkm::Vec<vtkm::Range, 3>>(PARA_RANGE);
      vtkm::Vec3f left_bottom{ {static_cast<vtkm::FloatDefault>(range[0].Min),},
                               {static_cast<vtkm::FloatDefault>(range[1].Min),},
                               {static_cast<vtkm::FloatDefault>(range[2].Min),} };
      vtkm::Vec3f right_top{ {static_cast<vtkm::FloatDefault>(range[0].Max),},
                             {static_cast<vtkm::FloatDefault>(range[1].Max),},
                             {static_cast<vtkm::FloatDefault>(range[2].Max),} };
      
      locator.SetRange(left_bottom, right_top);
      
      locator.SetCutOff(_para.GetParameter<Real>(PARA_CUTOFF));
      
      locator.SetRs(_para.GetParameter<Real>(PARA_R_CORE));

      locator.SetPosition(_para.GetFieldAsArrayHandle<Vec3f>(field::position));
  }

  void SetForceFunction(ContForceFunction& force_function)
  {
      auto rhomax = _para.GetParameter<Real>(EAM_PARA_RHOMAX);
      auto nrho = _para.GetParameter<Id>(EAM_PARA_NRHO);
      auto drho = _para.GetParameter<Real>(EAM_PARA_DRHO);
      auto nr = _para.GetParameter<Id>(EAM_PARA_NR);
      auto dr = _para.GetParameter<Real>(EAM_PARA_DR);
      force_function.SetEAMParameters(rhomax, nrho, drho, nr, dr);
  }

protected:
  Application& _app;
  Para& _para;
  std::shared_ptr<Executioner>& _executioner;
};
