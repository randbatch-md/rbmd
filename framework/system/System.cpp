#include "System.h"
#include "Application.h"
#include "InitCondition.h"

System::System(const Configuration& cfg)
  : Object(cfg)
  , _app(*(Get<Application*>("_app")))
{
}

void System::Init()
{
  _timer.Start();

  //add variable
  InitField();

  for (const auto& init : _app.GetInitConditionWarehouse())
  {
    init->Execute();
  }
}