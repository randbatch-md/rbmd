#pragma once
#include <vtkm/worklet/DispatcherMapTopology.h>
#include <vtkm/worklet/Keys.h>
#include <vtkm/worklet/ScatterCounting.h>
#include <vtkm/worklet/WorkletMapField.h>
#include <vtkm/worklet/WorkletMapTopology.h>
#include <vtkm/worklet/WorkletReduceByKey.h>

//namespace galois
//{
using vtkm::worklet::WorkletMapField;
using vtkm::worklet::WorkletVisitCellsWithPoints;
using vtkm::worklet::DispatcherMapTopology;
using vtkm::worklet::WorkletVisitCellsWithPoints;
using vtkm::worklet::Keys;
using vtkm::worklet::ScatterCounting;
using vtkm::worklet::WorkletReduceByKey;
//}
