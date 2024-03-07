#pragma once

class Application;
class JsonParser;

class Action
{
public:
  Action(Application& app);
  virtual ~Action() = default;

  virtual void Execute() = 0;

protected:
  Application& _app;
  JsonParser& _parser;
};
