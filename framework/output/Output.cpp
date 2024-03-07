#include "Output.h"
#include "Application.h"
#include "FieldName.h"

Output::Output(const Configuration& cfg)
  : Object(cfg)
  , _app(*(Get<Application*>("_app")))
  , _system(*(_app.GetSystem()))
  , _executioner(_app.GetExecutioner())
{

}

void Output::Init()
{
}
