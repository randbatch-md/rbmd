#include "MeshFreeFileInitCondition.h"

MeshFreeFileInitCondition::MeshFreeFileInitCondition(const Configuration& cfg)
  : InitCondition(cfg)
  , _file(Get<std::string>("file"))
{
  _file.open(_para.GetParameter<std::string>(PARA_READ_FILE));
}

//void MeshFreeFileInitCondition::InitParameter() 
//{
//}

void MeshFreeFileInitCondition::InitField() 
{
  _para.AddField(field::charge, ArrayHandle<Real>{});
  _para.AddField(field::velocity, ArrayHandle<Vec3f>{});
  _para.AddField(field::mass, ArrayHandle<Real>{});
  _para.AddField(field::molecule_id, ArrayHandle<Id>{});
  _para.AddField(field::atom_id, ArrayHandle<Id>{});
  _para.AddField(field::position, ArrayHandle<Vec3f>{});
}

void MeshFreeFileInitCondition::SetParameters() {}
