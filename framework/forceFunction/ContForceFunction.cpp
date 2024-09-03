#include "forceFunction/ContForceFunction.h"
#include <vtkm/cont/ArrayHandle.h>

VTKM_CONT ExecForceFunction ContForceFunction::PrepareForExecution(vtkm::cont::DeviceAdapterId device,
                                       vtkm::cont::Token& token) const
{
  return ExecForceFunction(this->_cut_Off,
                           this->_alpha,
                           this->_volume,
                           this->_Vlength,
                           this->_box,
                           this->Kmax,
                           this->RBEP,
                           this->_rhomax,
                           this->_nrho,
                           this->_drho,
                           this->_nr,
                           this->_dr);
}
