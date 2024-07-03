#include "MeshFreeFileInitCondition.h"
#include "FieldName.h"
MeshFreeFileInitCondition::MeshFreeFileInitCondition(const Configuration& cfg)
  : InitCondition(cfg)
  , _file(Get<std::string>("file"))
{
}

void MeshFreeFileInitCondition::InitField()
{
  _para.AddField(field::charge, ArrayHandle<Real>{});
  _para.AddField(field::velocity, ArrayHandle<Vec3f>{});
  _para.AddField(field::mass, ArrayHandle<Real>{});
  _para.AddField(field::molecule_id, ArrayHandle<Id>{});
  _para.AddField(field::atom_id, ArrayHandle<Id>{});
  _para.AddField(field::position, ArrayHandle<Vec3f>{});
  _para.AddField(field::special_offsets, ArrayHandle<Id>{});
  _para.AddField(field::special_weights, ArrayHandle<Real>{});
  _para.AddField(field::special_ids, ArrayHandle<Id>{});
}

void MeshFreeFileInitCondition::SetParameters() {}