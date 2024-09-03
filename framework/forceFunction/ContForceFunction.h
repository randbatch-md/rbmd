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
                     const Vec3f& box,
                     const IdComponent& kmax)
  {
    this->_cut_Off = cut_off;
    this->_alpha = alpha;
    this->_volume = volume;
    this->_Vlength = vlength;
    this->_box = box;
    this->Kmax = kmax;
  }

  void SetEAMParameters(const Real& rhomax,
                        const Id& nrho,
                        const Real& drho,
                        const Id& nr,
                        const Real& dr)
  {
    this->_rhomax = rhomax;
    this->_nrho = nrho;
    this->_drho = drho;
    this->_nr = nr;
    this->_dr = dr;
  }

  VTKM_CONT
  ExecForceFunction PrepareForExecution(vtkm::cont::DeviceAdapterId device,
                                       vtkm::cont::Token& token) const;

  private:
  Real _cut_Off;
  Real _alpha;
  Real _volume;
  Real _Vlength;
  Vec3f _box;
  IdComponent Kmax;
  IdComponent RBEP;

  Real _rhomax;
  Id _nrho;
  Id _nr;
  Real _drho;
  Real _dr;
};