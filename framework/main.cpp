#include "MDApplication.h"
#include <iostream>
#include <string>

const std::string RBMD_VERSION = "1.0.0";
int main(int argc, char* argv[])
{
  if (std::string(argv[1]) == "--version")
  {
    std::cout << RBMD_VERSION << std::endl;
    return 0;
  }

  std::shared_ptr<Application> app = std::make_shared<MDApplication>(argc, argv);
  app->Run();

  return 0;
};