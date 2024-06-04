#include "MadelungSystem.h"
#include "Application.h"
#include "Executioner.h"
#include <vtkm/VectorAnalysis.h>
#include <vtkm/cont/Algorithm.h>
#include <vtkm/worklet/WorkletCellNeighborhood.h>
#include <vtkm/worklet/WorkletMapField.h>
#include <vtkm/worklet/WorkletPointNeighborhood.h>
#include "FieldName.h"

struct MadelungSystem::ComputeCoeffient : vtkm::worklet::WorkletMapField
{
  using ControlSignature = void(FieldIn pts, FieldIn charge, FieldOut force);
  using ExecutionSignature = void(_1, _2, _3);

  template<typename CoordType>
  VTKM_EXEC void operator()(const CoordType& pts, const Real& charge, Real& force) const
  {
    auto p = pts - _pts;
    Real dis = vtkm::Magnitude(p);
    if (dis < 0.0001)
      force = 0;
    else
      force = charge / dis;
  }
  Vec3f _pts{ 51, 51, 51 };
};


//RegisterObject(MadelungSystem);
MadelungSystem::MadelungSystem(const Configuration& cfg)
  : System(cfg)
  , _madelung_result(0.0f)
{

}

void MadelungSystem::Init()
{
  System::Init();
}

void MadelungSystem::Evolve() 
{
  auto&& charge = _para.GetFieldAsArrayHandle<Real>(field::charge);
  auto&& position = _para.GetFieldAsArrayHandle<Vec3f>(field::position);

  ArrayHandle<Real> force;
  vtkm::cont::Invoker{}(MadelungSystem::ComputeCoeffient{}, position, charge, force);

  _madelung_result = vtkm::cont::Algorithm::Reduce(force, vtkm::TypeTraits<Real>::ZeroInitialization());
  _para.SetParameter(gtest::madelung, _madelung_result);
  console::Info("Madelung常数：", _madelung_result);
  console::Info("MadelungSystem::Evolve运行完成，", "计算时间: ", _timer.GetElapsedTime(), " s");
}

void MadelungSystem::InitField() 
{
  _para.AddField(field::position, ArrayHandle<Vec3f>{});
  _para.AddField(field::charge, ArrayHandle<Real>{});
}