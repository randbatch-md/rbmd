﻿#pragma once 

#include "Output.h"
#include <string>

class System;
class Executioner;
class FileOutput : public Output
{
public:
  FileOutput(const Configuration& cfg);
  virtual ~FileOutput(){};

  void Init() override;
  void Execute() override;

protected:
  virtual bool ShouldOutput();
protected:
  System& _system;
  Executioner& _executioner;
  int _interval;
};

