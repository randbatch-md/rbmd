#include "FieldName.h"
#include "ExecutionNPT.h"

ExecutionNPT::ExecutionNPT(const Configuration& cfg)
  : ExecutionMD(cfg)
{
  ExecutionMD::SetParameters();
}

void ExecutionNPT::Init() {}

void ExecutionNPT::Execute()
{
  //_timer.Start();
  
}

void ExecutionNPT::SetParameters()
{
}