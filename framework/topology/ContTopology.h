#pragma once
#include <vtkm/cont/ExecutionObjectBase.h>
#include "topology/ExecTopology.h"
#include "FieldName.h"
#include "Types.h"
#include <vtkm/cont/ConvertNumComponentsToOffsets.h>
#include <vtkm/cont/ArrayHandleGroupVecVariable.h>


class ContTopology : vtkm::cont::ExecutionObjectBase
{
public:
  using GroupVecType = typename vtkm::cont::ArrayHandleGroupVecVariable<vtkm::cont::ArrayHandle<vtkm::Id>,
                                                                        vtkm::cont::ArrayHandle<vtkm::Id>>;

  void SetSourceAndOffsets(const vtkm::cont::ArrayHandle<vtkm::Id>& source_array,
                           const vtkm::cont::ArrayHandle<vtkm::Id>& offsets_array)
  {
    this->_group_vec_array = GroupVecType(source_array, offsets_array);
  }
  void SetMolecularId(const vtkm::cont::ArrayHandle<Id>& molecular_id)
  {
    this->_molecular_id = molecular_id;
  }

  void SetAtomsType(const vtkm::cont::ArrayHandle<Id>& atoms_type)
  {
    this->_atoms_type = atoms_type;
  }
  void SetCharge(const vtkm::cont::ArrayHandle<Real>& charge) 
  { 
      this->_charge = charge;
  }

  void SetEpsAndSigma(const vtkm::cont::ArrayHandle<Real>& epsilon,
                      const vtkm::cont::ArrayHandle<Real>& sigma)
  {
    this->_epsilon = epsilon;
    this->_sigma = sigma;
  }

  //void SetNeighbourIdAndNum(const vtkm::cont::ArrayHandle<vtkm::Id>& neighbour_j_id,
  //                         const vtkm::cont::ArrayHandle<vtkm::Id>& neighbour_j_num)
  //{
  //  this->_group_vec_neighbour = GroupVecType(neighbour_j_id, neighbour_j_num);
  //}

  VTKM_CONT
  ExecTopology PrepareForExecution(vtkm::cont::DeviceAdapterId device,
                                        vtkm::cont::Token& token) const;

private:

  GroupVecType _group_vec_array;
  //GroupVecType _group_vec_neighbour;
  vtkm::cont::ArrayHandle<vtkm::Id> _molecular_id;
  vtkm::cont::ArrayHandle<vtkm::Id> _atoms_type;
  vtkm::cont::ArrayHandle<Real> _charge;
  vtkm::cont::ArrayHandle<Real> _epsilon;
  vtkm::cont::ArrayHandle<Real> _sigma;
};