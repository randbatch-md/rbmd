#include "MeshFreeFileInitCondition.h"

MeshFreeFileInitCondition::MeshFreeFileInitCondition(const Configuration& cfg)
  : InitCondition(cfg)
  , _file(Get<std::string>("file"))
{

}
