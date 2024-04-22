#include "RBMD.h"

RBMD::RBMD(int argc, char* argv[]) 
{
  int iargc=1;
  if (argc==2)
  {
    if (std::string(argv[iargc]) == "--version")
    {
      OutputVersion();
    }
    else if (std::string(argv[iargc]) == "-help")
    {
      std::cout << "---Reference website: https://www.randbatch.com/guide/CASE_STUDIES.html---" << std::endl;
    }
    else
    {
      std::cout << "---Please enter the configuration command---" << std::endl
                << "---You can enter '-help' to query---" << std::endl;
    }
  }
  else if (argc > 2)
  {
    if (std::string(argv[iargc]) == "-j" && argc > 3)
    {
      std::cout << "--- The formate is: '-j'+ '*.json'---" << std::endl;
    }
    else
    {
      std::cout << "---Please enter the configuration command---" << std::endl;
    }
  }
  else
  {
    std::cout << "---Please enter the configuration command---" << std::endl;
  }
}

void RBMD::OutputVersion() 
{
  const std::string RBMD_VERSION = "1.0.0";
  std::cout << RBMD_VERSION << std::endl;
  exit(0);
}
