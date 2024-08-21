#pragma once
#include <string>
class InitGlobal
{
public:
  InitGlobal(const std::string& device,int argc, char** argv);
  void KokkosInit(int argc, char** argv);
  void MPIInit(int argc, char** argv);
  ~InitGlobal() = default;

protected:
};
