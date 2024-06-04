#include "MolecularSystem.h"
#include "FieldName.h"
#include <vtkm/cont/Algorithm.h>
#include "system/worklet/SystemWorklet.h"
#include <vtkm/cont/ArrayHandleGroupVec.h>
#include "system/worklet/MolecularWorklet.h"

MolecularSystem::MolecularSystem(const Configuration& cfg)
  : MeshFreeSystem(cfg)
{
}

void MolecularSystem::Init() 
{
  MeshFreeSystem::Init();
  SetForceFunction();
  SetTopology();
}

void MolecularSystem::InitialCondition() 
{
	// 物理场初始化
  _charge = _para.GetFieldAsArrayHandle<Real>(field::charge);
  _topology.SetCharge(_charge);
  _velocity = _para.GetFieldAsArrayHandle<Vec3f>(field::velocity);
  _mass = _para.GetFieldAsArrayHandle<Real>(field::mass);

  _molecule_id = _para.GetFieldAsArrayHandle<Id>(field::molecule_id);
  _atoms_id = _para.GetFieldAsArrayHandle<Id>(field::atom_id);
}

void MolecularSystem::SetForceFunction() {}

void MolecularSystem::SetTopology() {}

void MolecularSystem::InitField() 
{
  MeshFreeSystem::InitField();
  _para.AddField(field::charge, ArrayHandle<Real>{});
  _para.AddField(field::velocity, ArrayHandle<Vec3f>{});
  _para.AddField(field::mass, ArrayHandle<Real>{});
  _para.AddField(field::molecule_id, ArrayHandle<Id>{});
  _para.AddField(field::atom_id, ArrayHandle<Id>{});
}

std::vector<Vec2f> MolecularSystem::ComputeChargeStructureFactorRBE(Real& _Vlength,
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

std::vector<Vec2f> MolecularSystem::ComputeChargeStructureFactorEwald(Real& _Vlength,
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


void MolecularSystem::ComputeRBEEleForce(ArrayHandle<Vec3f>& psample,
                                         IdComponent& RBE_P,
                                         ArrayHandle<Vec3f>& RBE_ele_force)
{
  auto Vlength = _para.GetParameter<Real>(PARA_VLENGTH);
  auto rhok = ComputeChargeStructureFactorRBE(Vlength, psample);
  ArrayHandle<Vec2f> whole_rhok = vtkm::cont::make_ArrayHandle(rhok);
  SystemWorklet::ComputeNewRBEForce(RBE_P, _atoms_id, psample, whole_rhok, _force_function,  _topology,  _locator, RBE_ele_force);


}

void MolecularSystem::ComputeEwaldEleForce(IdComponent& Kmax,
                                           ArrayHandle<Vec3f>& Ewald_ele_force)
{
  auto Vlength = _para.GetParameter<Real>(PARA_VLENGTH);
  auto rhok = ComputeChargeStructureFactorEwald(Vlength, Kmax);
  ArrayHandle<Vec2f>  whole_rhok = vtkm::cont::make_ArrayHandle(rhok);
    SystemWorklet::ComputeNewFarElectrostatics(Kmax,
                                               _atoms_id,
                                               whole_rhok,
                                               _force_function,
                                               _topology,
                                               _locator,
                                               Ewald_ele_force); 
}