//==================================================================================
//  RBMD 1.0 is developed for random batch molecular dynamics calculation.
//
//  Copyright(C) 2024 SHANGHAI JIAOTONG UNIVERSITY CHONGQING RESEARCH INSTITUTE
//
//  This program is free software : you can redistribute it and /or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with this program.If not, see < https://www.gnu.org/licenses/>.
//
//  Contact Email : [your - email@example.com]
//==================================================================================

#pragma once
#include "Object.h"
#include "DataObject.h"
#include "Para.h"

class InitCondition : public Object
{
public:
  InitCondition(const Configuration& cfg)
    : Object(cfg)
    , _app(*Get<Application*>("_app"))
    , _para(*(_app.GetParameter()))
  {
  };

  ~InitCondition() = default;

  virtual void Execute() = 0;
  virtual void InitField() = 0;
  virtual void SetParameters() = 0;

protected:
  virtual void UpdateField()=0;

protected:
  Application& _app;
  Para& _para;
};