#pragma once
#include "InitCondition.h"
#include "forceFunction/ContForceFunction.h"
#include "locator/ContPointLocator.h"
#include "topology/ContTopology.h"
#include "DataObject.h"
#include "FieldName.h"

class MeshFreeCondition: public InitCondition
{
public:
  MeshFreeCondition(const Configuration& cfg);
  ~MeshFreeCondition();
  void Execute() override;
  void UpdateField() override;
  void InitField() override;
  void SetParameters() override;

protected:
  void SetTopology(ContTopology& topology) 
  {
    //auto pts_type =  GetFieldAsArrayHandle<Id>(field::pts_type );
    //auto epsilon =  GetFieldAsArrayHandle<Real>(field::epsilon);
    //auto sigma =  GetFieldAsArrayHandle<Real>(field::sigma);
    //auto charge =  GetFieldAsArrayHandle<Real>(field::charge );
    //auto molecule_id =  GetFieldAsArrayHandle<Id>(field::molecule_id);
    auto pts_type = _para.GetFieldAsArrayHandle<Id>(field::pts_type);
    auto epsilon = _para.GetFieldAsArrayHandle<Real>(field::epsilon);
    auto sigma = _para.GetFieldAsArrayHandle<Real>(field::sigma);
    auto charge = _para.GetFieldAsArrayHandle<Real>(field::charge);
    auto molecule_id = _para.GetFieldAsArrayHandle<Id>(field::molecule_id);
    topology.SetAtomsType(pts_type);
    topology.SetEpsAndSigma(epsilon, sigma);
    topology.SetCharge(charge);
    topology.SetMolecularId(molecule_id);
  }

  void SetLocator(ContPointLocator& locator)
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

    locator.SetRange(left_bottom, right_top);

    locator.SetCutOff(_para.GetParameter<Real>(PARA_CUTOFF));

    locator.SetRs(_para.GetParameter<Real>(PARA_R_CORE));

    locator.SetPosition(_para.GetFieldAsArrayHandle<Vec3f>(field::position));
  }

private:
  void InitPosition();
  void InitMassAndVelocity();
  void InitPositionFlag();
  void InitId();
  void InitParameters();

protected:
  ArrayHandle<Vec3f> _position;
  std::vector<int> _dims;
  std::vector<Real> _x_range;
  std::vector<Real> _y_range;
  std::vector<Real> _z_range;
};