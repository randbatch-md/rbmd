#include "FieldName.h"
#include "ReadData.h"
ReadData::ReadData(const Configuration& cfg)
  : InitConfiguration(cfg) 
{
}

void ReadData::Execute()
{
  InitField(); 
  SetPara();
}

void ReadData::SetPara()
{
  SetParameter(gtest::velocity_type, Get<std::string>("velocity_type"));
  SetParameter(ATOM_STYLE, Get<std::string>("atom_style"));
  SetParameter(PARA_UNIT, Get<std::string>("unit"));
  SetParameter(PARA_READ_FILE, Get<std::string>("file"));
}

void ReadData::InitField()
{
  AddField(field::charge, ArrayHandle<Real>{});
  AddField(field::velocity, ArrayHandle<Vec3f>{});
  AddField(field::mass, ArrayHandle<Real>{});
  AddField(field::molecule_id, ArrayHandle<Id>{});
  AddField(field::atom_id, ArrayHandle<Id>{});
  AddField(field::position, ArrayHandle<Vec3f>{});

  AddField(field::position_flag, ArrayHandle<Id3>{});
  AddField(field::bond_atom_id, ArrayHandle<Id>{});
  AddField(field::bond_type, ArrayHandle<Id>{});
  AddField(field::bond_coeffs_k, ArrayHandle<Real>{});
  AddField(field::bond_coeffs_equilibrium, ArrayHandle<Real>{});
  AddField(field::angle_atom_id, ArrayHandle<Id>{});
  AddField(field::angle_type, ArrayHandle<Id>{});
  AddField(field::angle_coeffs_k, ArrayHandle<Real>{});
  AddField(field::angle_coeffs_equilibrium, ArrayHandle<Real>{});
  AddField(field::atom_id_center, ArrayHandle<Id>{});
  AddField(field::atom_id_target, ArrayHandle<Id>{});
  AddField(field::pts_type, ArrayHandle<Id>{});
  AddField(field::center_position, ArrayHandle<Vec3f>{});
  AddField(field::target_position, ArrayHandle<Vec3f>{});
  AddField(field::epsilon, ArrayHandle<Real>{});
  AddField(field::sigma, ArrayHandle<Real>{});
  AddField(field::dihedrals_atom_id, ArrayHandle<Id>{});
  AddField(field::dihedrals_type, ArrayHandle<Id>{});
  AddField(field::signal_atoms_id, ArrayHandle<Id>{});
  AddField(field::special_source_array, ArrayHandle<Id>{});
  AddField(field::special_offsets_array, ArrayHandle<Id>{});
}
