#pragma once
#include <vtkm/cont/ExecutionObjectBase.h>
#include "staticTable/ExecStaticTable.h"
#include "Types.h"


class ContStaticTable : vtkm::cont::ExecutionObjectBase
{
public:
  void SetTableIndex1(const vtkm::cont::ArrayHandle<IdComponent>& table_index1)
  {
    this->_table_index_1 = table_index1;
  }

  void SetTableRij1(const vtkm::cont::ArrayHandle<Real>& table_rij1)
  {
    this->_table_rij_1 = table_rij1;
  }

  void SetTabledRij1(const vtkm::cont::ArrayHandle<Real>& table_drij1)
  {
    this->_table_drij_1 = table_drij1;
  }

  void SetTableFunctionRij1(const vtkm::cont::ArrayHandle<Real>& table_function_rij1)
  {
    this->_table_function_rij_1 = table_function_rij1;
  }

  void SetTabledFunctionRij1(const vtkm::cont::ArrayHandle<Real>& table_dfunction_rij1)
  {
    this->_table_dfunction_rij_1 = table_dfunction_rij1;
  }

  void SetTableIndex2(const vtkm::cont::ArrayHandle<IdComponent>& table_index2)
  {
    this->_table_index_2 = table_index2;
  }

  void SetTableRij2(const vtkm::cont::ArrayHandle<Real>& table_rij2)
  {
    this->_table_rij_2 = table_rij2;
  }

  void SetTabledRij2(const vtkm::cont::ArrayHandle<Real>& table_drij2)
  {
    this->_table_drij_2 = table_drij2;
  }

  void SetTableFunctionRij2(const vtkm::cont::ArrayHandle<Real>& table_function_rij2)
  {
    this->_table_function_rij_2 = table_function_rij2;
  }

  void SetTabledFunctionRij2(const vtkm::cont::ArrayHandle<Real>& table_dfunction_rij2)
  {
    this->_table_dfunction_rij_2 = table_dfunction_rij2;
  }

  VTKM_CONT
  ExecStaticTable PrepareForExecution(vtkm::cont::DeviceAdapterId device,
                                       vtkm::cont::Token& token) const;

  private:
  vtkm::cont::ArrayHandle<IdComponent> _table_index_1;
  vtkm::cont::ArrayHandle<Real> _table_rij_1;
  vtkm::cont::ArrayHandle<Real> _table_drij_1;
  vtkm::cont::ArrayHandle<Real> _table_function_rij_1;
  vtkm::cont::ArrayHandle<Real> _table_dfunction_rij_1;
  vtkm::cont::ArrayHandle<IdComponent> _table_index_2;
  vtkm::cont::ArrayHandle<Real> _table_rij_2;
  vtkm::cont::ArrayHandle<Real> _table_drij_2;
  vtkm::cont::ArrayHandle<Real> _table_function_rij_2;
  vtkm::cont::ArrayHandle<Real> _table_dfunction_rij_2;


};