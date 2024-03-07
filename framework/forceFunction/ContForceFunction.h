#pragma once
#include <vtkm/cont/ExecutionObjectBase.h>
#include "forceFunction/ExecForceFunction.h"
#include "Types.h"


class ContForceFunction : vtkm::cont::ExecutionObjectBase
{
public:

  void SetParameters(const Real& cut_off,
                     const Real& alpha,
                     const Real& volume,
                     const Real& vlength,
                     const IdComponent& kmax)
  {
    this->_cut_Off = cut_off;
    this->_alpha = alpha;
    this->_volume = volume;
    this->_Vlength = vlength;
    this->Kmax = kmax;
  }

  VTKM_CONT
  ExecForceFunction PrepareForExecution(vtkm::cont::DeviceAdapterId device,
                                       vtkm::cont::Token& token) const;

  private:
  Real _cut_Off;
  Real _alpha;
  Real _volume;
  Real _Vlength;
  IdComponent Kmax;
  IdComponent RBEP;
};