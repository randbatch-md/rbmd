#include "FieldName.h"
#include "ExecutionNVE.h"

ExecutionNVE::ExecutionNVE(const Configuration& cfg)
  : ExecutionMD(cfg)
{
  Execution::SetParameters();
}

void ExecutionNVE::Init() {}

void ExecutionNVE::Execute()
{

}

void ExecutionNVE::SetParameters()
{
}