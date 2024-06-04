#pragma once
#include "Object.h"
#include "DataObject.h"
#include "Para.h"

class InitCondition : public Object
{
public:
  InitCondition(const Configuration& cfg)
    : Object(cfg)
    , _app(*Get<Application*>("app"))
    , _para(*(_app.GetParameter())) // 这里取出来 system为空，因为已经没有system一说了
  {  
  };
  ~InitCondition() = default;

  virtual void Execute() = 0;
  //virtual void InitParameter()=0;
  virtual void InitField() = 0;
  virtual void SetParameters() = 0;

protected:
  virtual void UpdateField()=0;

protected:
  Application& _app;
  Para& _para;
};