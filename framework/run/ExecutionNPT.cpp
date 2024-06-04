#include "FieldName.h"
#include "ExecutionNPT.h"

ExecutionNPT::ExecutionNPT(const Configuration& cfg)
  : Execution(cfg)
{
  Execution::SetParameters();
}

void ExecutionNPT::Execute()
{
}

void ExecutionNPT::SetParameters()
{
}