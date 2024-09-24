#include "Output.h"
#include "Application.h"
#include "FieldName.h"

Output::Output(const Configuration& cfg)
  : Object(cfg)
  , _app(*(Get<Application*>("_app")))
  , _para(*(_app.GetParameter()))
  , _executioner(_app.GetExecutioner())
{

}

void Output::Init()
{
}
