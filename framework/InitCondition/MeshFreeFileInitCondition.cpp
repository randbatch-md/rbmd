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