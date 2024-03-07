#include "Action.h"
#include "Application.h"
#include "JsonParser.h"

Action::Action(Application& app)
  : _app(app)
  , _parser(*(app.GetParser()))
{

}
