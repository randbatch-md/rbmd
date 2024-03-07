#include "SetupDeviceAction.h"
#include "Application.h"

#define VTKM_NO_ERROR_ON_MIXED_CUDA_CXX_TAG

#include <vtkm/cont/DeviceAdapter.h>
//#include <vtkm/cont/RuntimeDeviceInformation.h>
#include <vtkm/cont/RuntimeDeviceTracker.h>
#ifdef VTKM_ENABLE_TBB
#include "tbb/task_scheduler_init.h"
#endif //  VTKM_ENABLE_TBB

  SetupDeviceAction::SetupDeviceAction(Application& app)
  : Action(app)
  {

  }

void SetupDeviceAction::Execute()
{
  auto& parallelNode =  _parser.GetJsonNode("parallel");

  auto& device = _app.DeviceTag();
  auto type = parallelNode["type"].asString();
  if (type == "serial")
  {
    device.reset(new vtkm::cont::DeviceAdapterTagSerial);
  }
  else if (type == "tbb")
  {
    device.reset(new vtkm::cont::DeviceAdapterTagTBB);
  }
  else if (type == "cuda")
  {
    device.reset(new vtkm::cont::DeviceAdapterTagCuda);
  }
  else if (type == "hip")
  {
    device.reset(new vtkm::cont::DeviceAdapterTagKokkos);
  }
  else if (type == "openmp")
  {
    device.reset(new vtkm::cont::DeviceAdapterTagOpenMP);
  }
  else
  {
    console::Warning( "未知并行模式:", type, ", 切换为串行模式....");
    device.reset(new vtkm::cont::DeviceAdapterTagSerial);
  }

#ifdef VTKM_ENABLE_TBB
  int num_threads = parallelNode["num_threads"].asInt();
  int max_threads = tbb::task_scheduler_init::default_num_threads();
  if (num_threads > max_threads)
  {
    console::Warning( "系统核数 < 设置核数");
    num_threads = max_threads;
  }
  std::cout << " 系统最大核数：" << max_threads << " 当前运行核数：" << num_threads << std::endl;

  tbb::task_scheduler_init init(num_threads);
#endif

  #ifdef VTKM_ENABLE_OPENMP
  int num_threads = parallelNode["num_threads"].asInt();
  int max_threads = omp_get_max_threads();
  if (num_threads > max_threads)
  {
    console::Warning("系统核数 < 设置核数");
    num_threads = max_threads;
  }
  std::cout << " 系统最大核数：" << max_threads << " 当前运行核数：" << num_threads << std::endl;

  omp_set_num_threads(num_threads);
#endif

  //console::Success("SetupDeviceAction: ", device->GetName());
}

#undef VTKM_NO_ERROR_ON_MIXED_CUDA_CXX_TAG
