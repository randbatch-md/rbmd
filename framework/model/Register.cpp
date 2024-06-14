﻿#include "Register.h"
#include "Factory.h"
#include "LJInitCondition.h"
#include "MadelungInitCondition.h"
#include "ModelFileInitCondition.h"
#include "AtomsFileOutput.h"
#include "AtomsTableOutput.h"
#include "RDFOutput.h"
#include "ThermoOutput.h"
#include "TempOutput.h"
#include "MSDOutput.h"
#include "Progress.h"
#include "H2OSystem.h"
#include "LennardJonesSystem.h"
#include "MadelungSystem.h"
#include "SalineSolutionSystem.h"
#include "NaClSystem.h"
#include "ERFInitCondition.h"
#include "TrajectoryOutput.h"
#include "VACFOutput.h"
#include "PEOSystem.h"
#include "EAMSystem.h"
void RegisterObjectGlobal() 
{
  RegisterObject(LJInitCondition);
  RegisterObject(MadelungInitCondition);
  RegisterObject(ModelFileInitCondition);
  RegisterObject(AtomsFileOutput);
  RegisterObject(AtomsTableOutput);
  RegisterObject(ConsoleOutput);
  RegisterObject(RDFOutput);
  RegisterObject(ThermoOutput);
  RegisterObject(TempOutput);
  RegisterObject(MSDOutput);
  RegisterObject(TrajectoryOutput);
  RegisterObject(Progress);
  RegisterObject(H2OSystem);
  RegisterObject(LennardJonesSystem);
  RegisterObject(MadelungSystem);
  RegisterObject(SalineSolutionSystem);
  RegisterObject(NaClSystem);
  //RegisterObject(ERFInitCondition);
  RegisterObject(VACFOutput);
  RegisterObject(PEOSystem);
  RegisterObject(EAMSystem);  
}