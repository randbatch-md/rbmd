#pragma once

class InitGlobal
{
public:
  InitGlobal(int argc, char** argv);
  void KokkosInit(int argc, char** argv);
  ~InitGlobal() = default;

protected:
};