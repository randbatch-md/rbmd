#include "MadelungInitCondition.h"
#include "FieldName.h"
#include <vtkm/cont/DataSetBuilderUniform.h>
#include <vtkm/worklet/WorkletPointNeighborhood.h>
#include <vtkm/cont/Invoker.h>
#include <vtkm/cont/ArrayCopy.h>

struct MadelungInitCondition::SetChargeWorklet : vtkm::worklet::WorkletPointNeighborhood
{
  using ControlSignature = void(CellSetIn, FieldOut charge);
  using ExecutionSignature = void(Boundary, _2);

  template<typename T>
  VTKM_EXEC void operator()(const vtkm::exec::BoundaryState& boundary, T& charge) const
  {
    if ((boundary.IJK[0] + boundary.IJK[1] + boundary.IJK[2]) % 2 == 0)
      charge = -1;
    else
    {
      charge = 1;
    }
  }
};

MadelungInitCondition::MadelungInitCondition(const Configuration& cfg)
  : InitCondition(cfg)
  , _dims(GetVectorValue<int>("dims"))
  , _x_range(GetVectorValue<Real>("x_range"))
  , _y_range(GetVectorValue<Real>("y_range"))
  , _z_range(GetVectorValue<Real>("z_range"))
{

}

void MadelungInitCondition::Execute() 
{
  UpdateField();
}

void MadelungInitCondition::UpdateField() 
{
  //update position
  auto position = _para.GetFieldAsArrayHandle<Vec3f>(field::position);
  Vec3f start{ _x_range[0], _y_range[0], _z_range[0] };
  Vec3f end{ _x_range[1], _y_range[1], _z_range[1] };
  Vec3f space = end - start;
  Id3 point_dim = { _dims[0], _dims[1], _dims[2]};
  for (size_t i = 0; i < 3; i++)
  {
    space[i] = (end[i] - start[i]) / _dims[i];
  }
  auto data_set = vtkm::cont::DataSetBuilderUniform::Create(point_dim, start, space);
  vtkm::cont::ArrayCopy(data_set.GetCoordinateSystem().GetData(), position);
  
  //update charge
  auto charge = _para.GetFieldAsArrayHandle<Real>(field::charge);
  vtkm::cont::Invoker{}(SetChargeWorklet{}, data_set.GetCellSet(),charge);
}
