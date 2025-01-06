//==================================================================================
//  RBMD 2.2.0 is developed for random batch molecular dynamics calculation.
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


#pragma once

#include "Object.h"
#include "ContPointLocator.h"
#include "forceFunction/ContForceFunction.h"
#include "FieldName.h"
#include "topology/ContTopology.h"
#include "Para.h"
#include "vtkm/cont/ArrayHandle.h"
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
