#include "MDApplication.h"
#include <iostream>
#include <string>

#ifdef VTKM_HIP
#include <Kokkos_Core.hpp>
#endif // VTKM_HIP

int main(int argc, char* argv[])
{
#ifdef VTKM_HIP
#include <Kokkos_Core.hpp>
#endif // VTKM_HIP


  std::shared_ptr<Application> app = std::make_shared<MDApplication>(argc, argv);
  app->Run();

  return 0;
};