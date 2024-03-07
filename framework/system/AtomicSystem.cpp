﻿#include "AtomicSystem.h"
#include "FieldName.h"
#include "system/worklet/SystemWorklet.h"
#include <vtkm/cont/Algorithm.h>
AtomicSystem::AtomicSystem(const Configuration& cfg): MeshFreeSystem(cfg){}

void AtomicSystem::Init() 
{
  MeshFreeSystem::Init();
  SetForceFunction();
  SetTopology();

}

void AtomicSystem::InitialCondition() 
{
  SetCharge();
  _velocity = GetFieldAsArrayHandle<Vec3f>(field::velocity);
  _mass = GetFieldAsArrayHandle<Real>(field::mass);
  _atoms_id = GetFieldAsArrayHandle<Id>(field::atom_id);
}

std::vector<Vec2f> AtomicSystem::ComputeChargeStructureFactorRBE(Real& _Vlength,
                                                   ArrayHandle<Vec3f>& _psample)
{
  std::vector<Vec2f> rhok;
  auto p_number = _psample.GetNumberOfValues();
  ArrayHandle<Real> density_real;
  ArrayHandle<Real> density_image;
  for (Id i = 0; i < p_number; i++)
  {
    auto kl = _psample.ReadPortal().Get(i);
    kl = 2 * vtkm::Pi() * kl / _Vlength;
    SystemWorklet::ComputeChargeStructureFactorComponent(
      kl, _position, _charge, density_real, density_image);

    Real value_Re =
      vtkm::cont::Algorithm::Reduce(density_real, vtkm::TypeTraits<Real>::ZeroInitialization());
    Real value_Im =
      vtkm::cont::Algorithm::Reduce(density_image, vtkm::TypeTraits<Real>::ZeroInitialization());

    rhok.push_back({ value_Re, value_Im });
  }
  return rhok;
}
     
std::vector<Vec2f> AtomicSystem::ComputeChargeStructureFactorEwald(Real& _Vlength,
                                                                   IdComponent& Kmax)
{
  std::vector<Vec2f> rhok;
  ArrayHandle<Real> density_real;
  ArrayHandle<Real> density_image;
  for (Id i = -Kmax; i <= Kmax; i++)
  {
    for (Id j = -Kmax; j <= Kmax; j++)
    {
      for (Id k = -Kmax; k <= Kmax; k++)
      {
        if (!(i == 0 && j == 0 && k == 0))
        {
          Vec3f K = { Real(i), Real(j), Real(k) };
          K = 2 * vtkm::Pi() * K / _Vlength;
          SystemWorklet::ComputeChargeStructureFactorComponent(
            K, _position, _charge, density_real, density_image);
          Real value_Re = vtkm::cont::Algorithm::Reduce(
            density_real, vtkm::TypeTraits<Real>::ZeroInitialization());
          Real value_Im = vtkm::cont::Algorithm::Reduce(
            density_image, vtkm::TypeTraits<Real>::ZeroInitialization());

          rhok.push_back({ value_Re, value_Im });
        }
      }
    }
  }
  return rhok;
}

void AtomicSystem::ComputeRBEEleForce(ArrayHandle<Vec3f>& psample,
                                      IdComponent& RBE_P,
                                      ArrayHandle<Vec3f>& RBE_ele_force)
{
  auto Vlength = GetParameter<Real>(PARA_VLENGTH);
  auto rhok = ComputeChargeStructureFactorRBE(Vlength, psample);
  ArrayHandle<Vec2f> whole_rhok = vtkm::cont::make_ArrayHandle(rhok);
  SystemWorklet::ComputeNewRBEForce(RBE_P, _atoms_id, psample, whole_rhok, _force_function, _topology, _locator, RBE_ele_force);
}

void AtomicSystem::ComputeEwaldEleForce(IdComponent& Kmax,
                                        ArrayHandle<Vec3f>& Ewald_ele_force)
{
  auto Vlength = GetParameter<Real>(PARA_VLENGTH);
  auto rhok = ComputeChargeStructureFactorEwald(Vlength, Kmax);
  ArrayHandle<Vec2f> whole_rhok = vtkm::cont::make_ArrayHandle(rhok);

  whole_rhok = vtkm::cont::make_ArrayHandle(rhok);
      SystemWorklet::ComputeNewFarElectrostatics(
    Kmax, _atoms_id, whole_rhok, _force_function, _topology, _locator, Ewald_ele_force);
}

void AtomicSystem::InitField() 
{
  MeshFreeSystem::InitField();
  AddField(field::charge , ArrayHandle<Real>{});
  AddField(field::velocity, ArrayHandle<Vec3f>{});
  AddField(field::mass, ArrayHandle<Real>{});
  AddField(field::atom_id, ArrayHandle<Id>{});
  AddField(field::molecule_id, ArrayHandle<Id>{});
}
