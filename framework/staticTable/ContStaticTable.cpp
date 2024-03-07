#include "staticTable/ContStaticTable.h"
#include <vtkm/cont/ArrayHandle.h>

VTKM_CONT ExecStaticTable ContStaticTable::PrepareForExecution(vtkm::cont::DeviceAdapterId device,
                                       vtkm::cont::Token& token) const
{
  return ExecStaticTable(this->_table_index_1.PrepareForInput(device, token),
                           this->_table_rij_1.PrepareForInput(device, token),
                           this->_table_drij_1.PrepareForInput(device, token),
                           this->_table_function_rij_1.PrepareForInput(device, token),
                           this->_table_dfunction_rij_1.PrepareForInput(device, token),
                           this->_table_index_2.PrepareForInput(device, token),
                           this->_table_rij_2.PrepareForInput(device, token),
                           this->_table_drij_2.PrepareForInput(device, token),
                           this->_table_function_rij_2.PrepareForInput(device, token),
                           this->_table_dfunction_rij_2.PrepareForInput(device, token)
                           );
}