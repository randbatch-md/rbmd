#include "FieldName.h"
#include "ExecutionNVT.h"

ExecutionNVT::ExecutionNVT(const Configuration& cfg)
  : ExecutionMD(cfg)
{
  Execution::SetParameters();
}

void ExecutionNVT::Init() 
{
  _timer.Start();

}

void ExecutionNVT::Execute()
{

}

void ExecutionNVT::SetParameters()
{
}