#include "ModelFileInitCondition.h"
#include <vtkm/cont/DataSetBuilderUniform.h>
#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <algorithm>
#include <vtkm/cont/ArrayCopy.h>
#include "FieldName.h"
#include <vtkm/cont/ConvertNumComponentsToOffsets.h>
#include <vtkm/cont/ArrayHandleGroupVecVariable.h>

//RegisterObject(ModelFileInitCondition);
ModelFileInitCondition::ModelFileInitCondition(const Configuration& cfg)
  : MeshFreeFileInitCondition(cfg)
{
  group_id_init.insert(std::make_pair("wat", MolecularType::H20));
  group_id_init.insert(std::make_pair("nacl", MolecularType::NACL));
  SetParameters();
}

void ModelFileInitCondition::Execute() 
{
  try
  {
    //SetParameters();
    InitField();
    ReadDataFile(_file);
    _file.close();
    UpdateField();
  }
  catch (const std::exception& e)
  {
    _file.close();
    std::cout << e.what() << std::endl;
  }
}

void ModelFileInitCondition::ReadDataFile(std::ifstream& file)
{

  if (file.is_open())
  {
    Header(file);
    Parser(file);
    Velocity();
  }
  else
  {
    std::cout << "Unable to open file";
  }
}

void ModelFileInitCondition::SetParameters()
{
  _para.SetParameter(gtest::velocity_type, Get<std::string>("velocity_type"));
  _para.SetParameter(ATOM_STYLE, Get<std::string>("atom_style"));
  _para.SetParameter(PARA_UNIT, Get<std::string>("unit"));
  _para.SetParameter(PARA_FILE_DIHEDRALS, false);
  _para.SetParameter(PARA_FILE_ANGLES, false);
}

void ModelFileInitCondition::InitField()
{
  MeshFreeFileInitCondition::InitField();
  _para.AddField(field::position_flag, ArrayHandle<Id3>{});
  _para.AddField(field::bond_atom_id, ArrayHandle<Id>{});
  _para.AddField(field::bond_type, ArrayHandle<Id>{});
  _para.AddField(field::bond_coeffs_k, ArrayHandle<Real>{});
  _para.AddField(field::bond_coeffs_equilibrium, ArrayHandle<Real>{});
  _para.AddField(field::angle_atom_id, ArrayHandle<Id>{});
  _para.AddField(field::angle_type, ArrayHandle<Id>{});
  _para.AddField(field::angle_coeffs_k, ArrayHandle<Real>{});
  _para.AddField(field::angle_coeffs_equilibrium, ArrayHandle<Real>{});
  _para.AddField(field::atom_id_center, ArrayHandle<Id>{});
  _para.AddField(field::atom_id_target, ArrayHandle<Id>{});
  _para.AddField(field::atom_pair_id, ArrayHandle<Id>{});
  _para.AddField(field::atoms_pair_type_offsets, ArrayHandle<Id>{});
  _para.AddField(field::atom_pair_position, ArrayHandle<Vec3f>{});
  _para.AddField(field::pts_type, ArrayHandle<Id>{});
  _para.AddField(field::center_position, ArrayHandle<Vec3f>{});
  _para.AddField(field::target_position, ArrayHandle<Vec3f>{});
  _para.AddField(field::epsilon, ArrayHandle<Real>{});
  _para.AddField(field::sigma, ArrayHandle<Real>{});
  _para.AddField(field::signal_atoms_id, ArrayHandle<Id>{});
  _para.AddField(field::special_source_array, ArrayHandle<Id>{});
  _para.AddField(field::special_offsets_array, ArrayHandle<Id>{});

  _para.AddField(field::dihedrals_atom_id, ArrayHandle<Id>{});
  _para.AddField(field::dihedrals_type, ArrayHandle<Id>{});
  _para.AddField(field::dihedrals_coeffs_k, ArrayHandle<Real>{});
  _para.AddField(field::dihedrals_coeffs_sign, ArrayHandle<vtkm::IdComponent>{});
  _para.AddField(field::dihedrals_coeffs_multiplicity, ArrayHandle<vtkm::IdComponent>{});
}

void ModelFileInitCondition::Header(std::ifstream& file)
{
  std::string line;
  Real x_low, x_high, y_low, y_high, z_low, z_high;
  while (getline(file, line))
  {
    if (!line.empty() && line != "\r")
    {
      if (line.find("atoms") != std::string::npos)
      {
        std::istringstream iss(line);
        iss >> _header._num_atoms;
      }
      else if (line.find("bonds") != std::string::npos)
      {
        std::istringstream iss(line);
        iss >> _header._num_bonds;
      }
      else if (line.find("angles") != std::string::npos)
      {
        std::istringstream iss(line);
        iss >> _header._num_angles;
      }
      else if (line.find("dihedrals") != std::string::npos)
      {
        std::istringstream iss(line);
        iss >> _header._num_dihedrals;
      }
      else if (line.find("impropers") != std::string::npos)
      {
        std::istringstream iss(line);
        iss >> _header._num_impropers;
      }
      else if (line.find("atom types") != std::string::npos)
      {
        std::istringstream iss(line);
        iss >> _header._num_atoms_type;
      }
      else if (line.find("bond types") != std::string::npos)
      {
        std::istringstream iss(line);
        iss >> _header._num_bound_type;
      }
      else if (line.find("angle types") != std::string::npos)
      {
        std::istringstream iss(line);
        iss >> _header._num_angle_type;
      }
      else if (line.find("dihedral types") != std::string::npos)
      {
        std::istringstream iss(line);
        iss >> _header._num_dihedrals_type;
      }
      else if (line.find("xlo xhi") != std::string::npos)
      {
        std::istringstream iss(line);
        iss >> x_low >> x_high;
      }
      else if (line.find("ylo yhi") != std::string::npos)
      {
        std::istringstream iss(line);
        iss >> y_low >> y_high;
      }
      else if (line.find("zlo zhi") != std::string::npos)
      {
        std::istringstream iss(line);
        iss >> z_low >> z_high;
        break;
      }
    }
  }

  auto xLength = x_high - x_low;
  auto yLength = y_high - y_low;
  auto zLength = z_high - z_low;
  _para.SetParameter(PARA_VLENGTH, xLength);
  _para.SetParameter(PARA_VOLUME, xLength * yLength * zLength);
  _para.SetParameter(PARA_RHO, _header._num_atoms / (xLength * yLength * zLength));
  _para.SetParameter(NTYPES, _header._num_atoms_type);
  auto cut_off = _para.GetParameter<Real>(PARA_CUTOFF);
  auto bin_number = Id3{
    static_cast<int>(xLength / cut_off),
    static_cast<int>(yLength / cut_off),
    static_cast<int>(zLength / cut_off),
  };
  auto range = vtkm::Vec<vtkm::Range, 3>{
    { x_low, x_high },
    { y_low, y_high },
    { z_low, z_high }};
  _para.SetParameter(PARA_BIN_NUMBER, bin_number);
  _para.SetParameter(PARA_RANGE, range);
}

void ModelFileInitCondition::Parser(std::ifstream& file) 
{
  std::string line;
  while (getline(file, line))
  {
    if (!line.empty() && line != "\r")
    {
      if (line.find("Masses") != std::string::npos)
      {
        Masses(file, line);
      }
      else if (line.find("Pair Coeffs") != std::string::npos)
      {
        PairCoeffs(file, line);
      }
      else if (line.find("Bond Coeffs") != std::string::npos)
      {
        BondCoeffs(file, line);
      }
      else if (line.find("Angle Coeffs") != std::string::npos)
      {
        AngleCoeffs(file, line);
      }
      else if (line.find("Dihedral Coeffs") != std::string::npos)
      {
        DihedralsCoeffs(file, line);
      }
      else if (line.find("group") != std::string::npos)
      {
        Group(line);
      }
      else if (line.find("Atoms") != std::string::npos)
      {
        Atoms(file, line);
      }
      else if (line.find("Bonds") != std::string::npos)
      {
        Bonds(file, line);
      }
      else if (line.find("Angles") != std::string::npos)
      {
        Angles(file, line);
        _para.SetParameter(PARA_FILE_ANGLES, true);
      }
      else if (line.find("Dihedrals") != std::string::npos)
      {
        Dihedrals(file, line);
        _para.SetParameter(PARA_FILE_DIHEDRALS, true);
      }
      else if (line.find("Velocities") != std::string::npos)
      {
        VelocityRead(file, line);
      }
    }
  }
}

void ModelFileInitCondition::Masses(std::ifstream& file, std::string line)
{
  _mass.resize(_header._num_atoms_type);
  std::vector<Real> atoms_id_masses(_header._num_atoms_type);
  for (int i = 0; i < _header._num_atoms_type && getline(file, line); i)
  {
    if (!line.empty() && line != "\r")
    {
      std::istringstream iss(line);
      iss >> atoms_id_masses[i] >> _mass[i];
      i++;
    }
  }

  auto unit = _para.GetParameter<std::string>(PARA_UNIT);
  if (unit == "LJ")
  {
    _mass = std::vector<Real>(_header._num_atoms_type, 1.0);
  }
}

void ModelFileInitCondition::PairCoeffs(std::ifstream& file, std::string& line)
{
  std::vector<Real> atoms_id_paircoeffs(_header._num_atoms_type);
  _sigma_atom.resize(_header._num_atoms_type);
  _eps_atom.resize(_header._num_atoms_type);

  for (int i = 0; i < _header._num_atoms_type && getline(file, line); i)
  {
    if (!line.empty() && line != "\r")
    {
      std::istringstream iss(line);
      iss >> atoms_id_paircoeffs[i] >> _eps_atom[i] >> _sigma_atom[i];
      i++;
    }
  }

  // EPS SIGMA Pts_type
  auto epsilon = _para.GetFieldAsArrayHandle<Real>(field::epsilon);
  vtkm::cont::ArrayCopy(vtkm::cont::make_ArrayHandle(_eps_atom), epsilon);

  auto sigma = _para.GetFieldAsArrayHandle<Real>(field::sigma);
  vtkm::cont::ArrayCopy(vtkm::cont::make_ArrayHandle(_sigma_atom), sigma);
}

void ModelFileInitCondition::BondCoeffs(std::ifstream& file, std::string& line)
{
  std::vector<Real> bond_id_bondcoeffs(_header._num_bound_type);
  _bond_coeffs_k.resize(_header._num_bound_type);
  _bond_coeffs_equilibrium.resize(_header._num_bound_type);

  for (int i = 0; i < _header._num_bound_type && getline(file, line); i)
  {
    if (!line.empty() && line != "\r")
    {
      std::istringstream iss(line);
      iss >> bond_id_bondcoeffs[i] >> _bond_coeffs_k[i] >> _bond_coeffs_equilibrium[i];
      i++;
    }
  }
}

void ModelFileInitCondition::AngleCoeffs(std::ifstream& file, std::string& line)
{
  std::vector<Real> bond_id_anglecoeffs(_header._num_angle_type);
  _angle_coeffs_k.resize(_header._num_angle_type);
  _angle_coeffs_equilibrium.resize(_header._num_angle_type);

  for (int i = 0; i < _header._num_angle_type && getline(file, line); i)
  {
    if (!line.empty() && line != "\r")
    {
      std::istringstream iss(line);
      iss >> bond_id_anglecoeffs[i] >> _angle_coeffs_k[i] >> _angle_coeffs_equilibrium[i];
      i++;
    }
  }
}

void ModelFileInitCondition::DihedralsCoeffs(std::ifstream& file, std::string& line)
{
  std::vector<Real> bond_id_dihedralscoeffs(_header._num_dihedrals_type);
  _dihedrals_coeffs_k.resize(_header._num_dihedrals_type);
  _dihedrals_coeffs_sign.resize(_header._num_dihedrals_type);
  _dihedrals_coeffs_multiplicity.resize(_header._num_dihedrals_type);

  for (int i = 0; i < _header._num_dihedrals_type && getline(file, line); i)
  {
    if (!line.empty() && line != "\r")
    {
      std::istringstream iss(line);
      iss >> bond_id_dihedralscoeffs[i] >> _dihedrals_coeffs_k[i] >> _dihedrals_coeffs_sign[i] >>
        _dihedrals_coeffs_multiplicity[i];
      i++;
    }
  }
}

void ModelFileInitCondition::Group(std::string& line)
{
  std::istringstream iss(line);
  std::string keyword;
  iss >> keyword;

  if ("group" == keyword)
  {
    std::string groupID;
    std::string type;
    std::vector<Id> atomtypes;

    iss >> groupID >> type;
    vtkm::Id atomType;
    if (group_id_init.find(groupID) != group_id_init.end())
    {
      while (iss >> atomType)
      {
        _group_info[atomType] = group_id_init[groupID];
      }
    }
    else
    {
      std::cout << "cannot find " << groupID << std::endl;
    }
  }
}
 
void ModelFileInitCondition::Atoms(std::ifstream& file, std::string& line)
{
  _atoms_id.resize(_header._num_atoms);
  _molecules_id.resize(_header._num_atoms);
  _atoms_type.resize(_header._num_atoms);
  _charge.resize(_header._num_atoms);
  _position.resize(_header._num_atoms);
  _molecular_type.resize(_header._num_atoms);
  
  vtkm::Id temp_type;
  vtkm::Id temp_id;
  vtkm::Id temp_molecule_id;
  auto atom_style = _para.GetParameter<std::string>(ATOM_STYLE);
  //
  if (atom_style == "atomic")
  {
    for (int i = 0; i < _header._num_atoms && getline(file, line); i)
    {
      if (!line.empty() && line != "\r")
      {
        std::istringstream iss(line);
        iss >> temp_id >> temp_type >> _position[temp_id - 1][0] >> _position[temp_id - 1][1] >>
          _position[temp_id - 1][2];
        _atoms_type[temp_id - 1] = temp_type - 1;
        _atoms_id[temp_id - 1] = temp_id - 1;
        AtomsMapInsert(_atoms_type[temp_id - 1], _atoms_id[temp_id - 1]);
        i++;
      }
    }
    SetAtomsFieldEAM();
  }
  //
  else if (atom_style == "full")
  {
    for (int i = 0; i < _header._num_atoms && getline(file, line); i)
    {
      if (!line.empty() && line != "\r")
      {
        std::istringstream iss(line);
        iss >> temp_id >> temp_molecule_id >> temp_type >> _charge[temp_id - 1] >>
          _position[temp_id - 1][0] >> _position[temp_id - 1][1] >> _position[temp_id - 1][2];
        _atoms_type[temp_id - 1] = temp_type - 1;
        _molecular_type[temp_id - 1] = _group_info[temp_type];
        _atoms_id[temp_id - 1] = temp_id - 1;
        _molecules_id[temp_id - 1] = temp_molecule_id - 1;
        MolecularMapInsert(_molecules_id[temp_id - 1], _atoms_id[temp_id - 1]);
        AtomsMapInsert(_atoms_type[temp_id - 1], _atoms_id[temp_id - 1]);
        AtomstoMolecular(_atoms_id[temp_id - 1], _molecules_id[temp_id - 1]);
        i++;
      }
    }
    SetAtomsField();
    SetMolecularGroup();
  }
  //
  else if (atom_style == "charge")
  {
    for (int i = 0; i < _header._num_atoms && getline(file, line); i)
    {
      if (!line.empty() && line != "\r")
      {
        std::istringstream iss(line);
        iss >> temp_id >> temp_type >> _charge[temp_id - 1] >> _position[temp_id - 1][0] >>
          _position[temp_id - 1][1] >> _position[temp_id - 1][2];
        _atoms_type[temp_id - 1] = temp_type - 1;
        _atoms_id[temp_id - 1] = temp_id - 1;
        AtomsMapInsert(_atoms_type[temp_id - 1], _atoms_id[temp_id - 1]);
        i++;
      }
    }
    SetAtomsField();
  }
  /*
  if (_header._num_bound_type!=0)
  {
    for (int i = 0; i < _header._num_atoms && getline(file, line); i)
    {
      if (!line.empty() && line != "\r")
      {
        std::istringstream iss(line);
        iss >> temp_id >> temp_molecule_id >> temp_type >> _charge[temp_id - 1] >>
          _position[temp_id - 1][0] >> _position[temp_id - 1][1] >> _position[temp_id - 1][2];
        _atoms_type[temp_id - 1] = temp_type - 1;
        _molecular_type[temp_id - 1] = _group_info[temp_type];
        _atoms_id[temp_id - 1] = temp_id - 1;
        _molecules_id[temp_id - 1] = temp_molecule_id - 1;
        MolecularMapInsert(_molecules_id[temp_id - 1], _atoms_id[temp_id - 1]);
        AtomsMapInsert(_atoms_type[temp_id - 1], _atoms_id[temp_id - 1]);
        AtomstoMolecular(_atoms_id[temp_id - 1], _molecules_id[temp_id - 1]);
        i++;
      }
    }
    SetMolecularGroup();
  }
  else if (_header._num_bound_type == 0)
  {
    for (int i = 0; i < _header._num_atoms && getline(file, line); i)
    {
      if (!line.empty() && line != "\r")
      {
        std::istringstream iss(line);
        iss >> temp_id >> temp_type >> _charge[temp_id - 1] >>
        _position[temp_id - 1][0] >> _position[temp_id - 1][1] >> _position[temp_id - 1][2];
        _atoms_type[temp_id - 1] = temp_type - 1;
        _atoms_id[temp_id - 1] = temp_id - 1;
        AtomsMapInsert(_atoms_type[temp_id - 1], _atoms_id[temp_id - 1]);
        i++;
      }
    }
  }
  SetAtomsField();
  */

}

void ModelFileInitCondition::Bonds(std::ifstream& file, std::string& line)
{
  std::unordered_set<vtkm::Id> bond_atoms_set;
  std::vector<vtkm::Id> bond_id_vec(_header._num_bonds);
  _bond_type.resize(_header._num_bonds);
  _bond_atoms_id.resize(_header._num_bonds * 2);
  int index = 0;
  vtkm::Id bond_type;
  vtkm::Id temp_id1;
  vtkm::Id temp_id2;
  for (int i = 0; i < _header._num_bonds && getline(file, line); i)
  {
    if (!line.empty() && line != "\r")
    {
      std::istringstream iss(line);
      iss >> bond_id_vec[i] >> bond_type >> temp_id1 >> temp_id2;
      _bond_type[i] = bond_type - 1;

      auto bond_id1 = temp_id1 - 1;
      auto bond_id2 = temp_id2 - 1;

      _bond_atoms_id[index++] = bond_id1;
      _bond_atoms_id[index++] = bond_id2;
      i++;

      bond_atoms_set.insert(bond_id1);
      bond_atoms_set.insert(bond_id2);

      //special
      _special_map.insert(std::make_pair(bond_id1, bond_id2));
      _special_map.insert(std::make_pair(bond_id2, bond_id1));
    }
  }
  for (const auto value : _atoms_id)
  {
    if (bond_atoms_set.find(value) == bond_atoms_set.end())
    {
      _signal_atoms_id.push_back(value);
    }
  }
  SetBondField();
}

void ModelFileInitCondition::Angles(std::ifstream& file, std::string& line)
{
  std::vector<vtkm::Id> angle_id_vec(_header._num_angles);
  _angle_type.resize(_header._num_angles);
  _angle_atoms_id.resize(_header._num_angles * 3);

  int index = 0;
  vtkm::Id angle_type;
  vtkm::Id temp_id1;
  vtkm::Id temp_id2;
  vtkm::Id temp_id3;
  for (int i = 0; i < _header._num_angles && getline(file, line); i)
  {
    if (!line.empty() && line != "\r")
    {
      std::istringstream iss(line);
      iss >> angle_id_vec[i] >> angle_type >> temp_id1 >> temp_id2 >> temp_id3;
      _angle_type[i] = angle_type - 1;
      _angle_atoms_id[index++] = temp_id1 -1;
      _angle_atoms_id[index++] = temp_id2 -1;
      _angle_atoms_id[index++] = temp_id3 -1;
      i++;
    }
  }
  SetAngleField();  
}

void ModelFileInitCondition::Dihedrals(std::ifstream& file, std::string& line)
{
  std::vector<vtkm::Id> dihedrals_id_vec(_header._num_dihedrals);
  _dihedrals_type.resize(_header._num_dihedrals);
  _dihedrals_atoms_id.resize(_header._num_dihedrals * 4);

  int index = 0;
  vtkm::Id dihedrals_type;
  vtkm::Id temp_id1;
  vtkm::Id temp_id2;
  vtkm::Id temp_id3;
  vtkm::Id temp_id4;
  for (int i = 0; i < _header._num_dihedrals && getline(file, line); i)
  {
    if (!line.empty() && line != "\r")
    {
      std::istringstream iss(line);
      iss >> dihedrals_id_vec[i] >> dihedrals_type >> temp_id1 >> temp_id2 >> temp_id3 >> temp_id4;
      _dihedrals_type[i] = dihedrals_type - 1;
      _dihedrals_atoms_id[index++] = temp_id1 - 1;
      _dihedrals_atoms_id[index++] = temp_id2 - 1;
      _dihedrals_atoms_id[index++] = temp_id3 - 1;
      _dihedrals_atoms_id[index++] = temp_id4 - 1;
      i++;
    }
    SetDihedralsField();
  }

  auto dihedrals_atoms_id = _para.GetFieldAsArrayHandle<Id>(field::dihedrals_atom_id);
  vtkm::cont::ArrayCopy(vtkm::cont::make_ArrayHandle(_dihedrals_atoms_id), dihedrals_atoms_id);

  auto dihedrals_type_array = _para.GetFieldAsArrayHandle<Id>(field::dihedrals_type);
  vtkm::cont::ArrayCopy(vtkm::cont::make_ArrayHandle(_dihedrals_type), dihedrals_type_array);
}

void ModelFileInitCondition::VelocityRead(std::ifstream& file, std::string& line)
{
  _velocity.resize(_header._num_atoms);
  vtkm::Id temp_id;
  for (int i = 0; i < _header._num_atoms && getline(file, line); i)
  {
    if (!line.empty() && line != "\r")
    {
      std::istringstream iss(line);
      iss >> temp_id >> _velocity[temp_id - 1][0] >> _velocity[temp_id - 1][1] >>
        _velocity[temp_id - 1][2]; 
      i++;
    }
  }

  auto velocity = _para.GetFieldAsArrayHandle<Vec3f>(field::velocity);
  vtkm::cont::ArrayCopy(vtkm::cont::make_ArrayHandle(_velocity), velocity);
}

void ModelFileInitCondition::Velocity()
{
  if (_velocity.size() == NULL)
  {
    _velocity.resize(_header._num_atoms);
    auto velocity_type = _para.GetParameter<std::string>(gtest::velocity_type);

    if (velocity_type == "GAUSS")
    {
      for (auto i = 0; i < _header._num_atoms; ++i)
      {
        Real u10 = RandomValue<Real>(0.0, 1.0);
        Real u20 = RandomValue<Real>(0.0, 1.0);
        Real z00 = vtkm::Sqrt(-2.0 * vtkm::Log(u10)) * vtkm::Cos(2.0 * vtkm::Pif() * u20);
        Real u11 = RandomValue<Real>(0.0, 1.0);
        Real u21 = RandomValue<Real>(0.0, 1.0);
        Real z01 = vtkm::Sqrt(-2.0 * vtkm::Log(u11)) * vtkm::Cos(2.0 * vtkm::Pif() * u21);
        Real u12 = RandomValue<Real>(0.0, 1.0);
        Real u22 = RandomValue<Real>(0.0, 1.0);
        Real z02 = vtkm::Sqrt(-2.0 * vtkm::Log(u12)) * vtkm::Cos(2.0 * vtkm::Pif() * u22);
        //Real sigmaz = 1;
        auto kB = 1.9872067 * vtkm::Pow(10.0, -3);
        Real sigmaz;
        auto unit = _para.GetParameter<std::string>(PARA_UNIT);
        if (unit == "REAL")
        {
          sigmaz = vtkm::Sqrt(kB * 298.0 / _mass[_atoms_type[i]]);
        }
        else
        {
          sigmaz = vtkm::Sqrt(1.0 * 1.0 / _mass[_atoms_type[i]]);
        }
        Vec3f gaussianr = Vec3f(z00 * sigmaz, z01 * sigmaz, z02 * sigmaz);
        _velocity[i] = gaussianr;
      }
    }
    else if (velocity_type == "ZERO" || "TEST")
    {
      for (auto i = 0; i < _header._num_atoms; ++i)
      {
        _velocity[i] = Vec3f(0.0, 0.0, 0.0);
      }
    }
    auto velocity = _para.GetFieldAsArrayHandle<Vec3f>(field::velocity);
    vtkm::cont::ArrayCopy(vtkm::cont::make_ArrayHandle(_velocity), velocity);
  }
  else
  {
    return;
  }
}

void ModelFileInitCondition::MolecularMapInsert(vtkm::Id& key, vtkm::Id& value)
{
  auto it = _molecular_map.find(key);
  if (it != _molecular_map.end())
  {
    it->second.push_back(value);
  }
  else
  {
    _molecular_map.insert(std::make_pair(key,std::vector<Id>{ value }));
  }
}

void ModelFileInitCondition::AtomsMapInsert(vtkm::Id& key, vtkm::Id& value)
{
  auto it = _atoms_map.find(key);
  if (it != _atoms_map.end())
  {
    it->second.push_back(value);
  }
  else
  {
    _atoms_map.insert(std::make_pair(key, std::vector<Id>{ value }));
  }
}

void ModelFileInitCondition::AtomstoMolecular(vtkm::Id& key, vtkm::Id& value)
{
  auto it = atom_to_molecular_map.find(key);
  if (it != atom_to_molecular_map.end())
  {
    it->second = value;
  }
  else
  {
    atom_to_molecular_map.insert(std::make_pair(key, value));
  }
}

void ModelFileInitCondition::SetMolecularGroup()
{
  //
  std::vector<std::vector<vtkm::Id>> atoms_gro;
  for (const auto atom_id : _atoms_id)
  {
    auto molecular_id = atom_to_molecular_map[atom_id];
    auto atoms_vec = _molecular_map[molecular_id];
    atoms_gro.emplace_back(atoms_vec);
  }

  std::vector<Id> atoms_vec_gro;
  std::vector<vtkm::Id> countVector;
  for (const std::vector<vtkm::Id>& innerVector : atoms_gro)
  {
    //扁平化
    atoms_vec_gro.insert(atoms_vec_gro.end(), innerVector.begin(), innerVector.end());
    //计数数组
    countVector.push_back(innerVector.size());
  }

  //sourceArray
  vtkm::cont::ArrayHandle<vtkm::Id> sourceArray = 
           vtkm::cont::make_ArrayHandle<vtkm::Id>(atoms_vec_gro);

  //offsetArray
  vtkm::cont::ArrayHandle<vtkm::Id> countArray =
           vtkm::cont::make_ArrayHandle<vtkm::Id>(countVector);

  vtkm::Id sourceArraySize;
  vtkm::cont::ArrayHandle<vtkm::Id> offsetsArray =
           vtkm::cont::ConvertNumComponentsToOffsets(countArray, sourceArraySize);

  auto source_array = _para.GetFieldAsArrayHandle<Id>(field::special_source_array);
  vtkm::cont::ArrayCopy(sourceArray, source_array);
  auto offsets_array = _para.GetFieldAsArrayHandle<Id>(field::special_offsets_array);
  vtkm::cont::ArrayCopy(offsetsArray, offsets_array);

  //
  auto molecule_id = _para.GetFieldAsArrayHandle<Id>(field::molecule_id);
  vtkm::cont::ArrayCopy(vtkm::cont::make_ArrayHandle(_molecules_id), molecule_id);
}

void ModelFileInitCondition::SetAtomsField() 
{
  // init Id
  auto atom_id = _para.GetFieldAsArrayHandle<Id>(field::atom_id);
  vtkm::cont::ArrayCopy(vtkm::cont::make_ArrayHandle(_atoms_id), atom_id);

  // init position
  ArrayHandle<Vec3f> position = _para.GetFieldAsArrayHandle<Vec3f>(field::position);
  vtkm::cont::ArrayCopy(vtkm::cont::make_ArrayHandle(_position), position);

  // init mass
  auto mass = _para.GetFieldAsArrayHandle<Real>(field::mass);
  mass.Allocate(_header._num_atoms);
  auto&& mass_protol = mass.WritePortal();
  for (int i = 0; i < _header._num_atoms; ++i)
  {
    mass_protol.Set(i, _mass[_atoms_type[i]]);
  }

  // init charge
  auto charge = _para.GetFieldAsArrayHandle<Real>(field::charge);
  vtkm::cont::ArrayCopy(vtkm::cont::make_ArrayHandle(_charge), charge);

  auto pts_type = _para.GetFieldAsArrayHandle<Id>(field::pts_type);
  vtkm::cont::ArrayCopy(vtkm::cont::make_ArrayHandle(_atoms_type), pts_type);

  auto init_way = _para.GetParameter<std::string>(PARA_INIT_WAY);
  if (init_way == "inbuild")
  {
    auto center_type = _para.GetParameter<IdComponent>(PARA_CENTER_TYPE);
    auto atoms_center = _atoms_map[center_type];
    auto atoms_id_center = _para.GetFieldAsArrayHandle<Id>(field::atom_id_center);
    vtkm::cont::ArrayCopy(vtkm::cont::make_ArrayHandle(atoms_center), atoms_id_center);
    
    auto target_type = _para.GetParameter<IdComponent>(PARA_TARGET_TYPE);
    auto atoms_target = _atoms_map[target_type];
    auto atoms_id_target = _para.GetFieldAsArrayHandle<Id>(field::atom_id_target);
    vtkm::cont::ArrayCopy(vtkm::cont::make_ArrayHandle(atoms_target), atoms_id_target);
  }
  else
  {
    auto atom_pair_type = _para.GetParameter<std::vector<int>>(PARA_ATOMS_PAIR_TYPE);
    std::vector<Id> atom_pair_id_temp;
    std::vector<Id> offsets;
    vtkm::Id total_size = 0;

    for (int i = 0; i < atom_pair_type.size(); i++)
    {
      auto atoms_type_id = _atoms_map[atom_pair_type[i]];
      atom_pair_id_temp.insert(atom_pair_id_temp.end(), atoms_type_id.begin(), atoms_type_id.end());

      offsets.push_back(total_size);
      total_size += atoms_type_id.size();
    }
    offsets.push_back(total_size);
    auto atom_pair_id = _para.GetFieldAsArrayHandle<Id>(field::atom_pair_id);
    vtkm::cont::ArrayCopy(vtkm::cont::make_ArrayHandle(atom_pair_id_temp), atom_pair_id);
    auto atoms_pair_type_offsets = _para.GetFieldAsArrayHandle<Id>(field::atoms_pair_type_offsets);
    vtkm::cont::ArrayCopy(vtkm::cont::make_ArrayHandle(offsets), atoms_pair_type_offsets);
  }


  // init position_flag
  auto position_flag = _para.GetFieldAsArrayHandle<Id3>(field::position_flag);
  position_flag.AllocateAndFill(_header._num_atoms, { 0, 0, 0 });
}

void ModelFileInitCondition::SetBondField() 
{
  auto signal_atoms_id = _para.GetFieldAsArrayHandle<Id>(field::signal_atoms_id);
  vtkm::cont::ArrayCopy(vtkm::cont::make_ArrayHandle(_signal_atoms_id), signal_atoms_id);

  auto bond_atom_id = _para.GetFieldAsArrayHandle<Id>(field::bond_atom_id);
  vtkm::cont::ArrayCopy(vtkm::cont::make_ArrayHandle(_bond_atoms_id), bond_atom_id);

  auto bond_type_handle = _para.GetFieldAsArrayHandle<Id>(field::bond_type);
  vtkm::cont::ArrayCopy(vtkm::cont::make_ArrayHandle(_bond_type), bond_type_handle);

  std::vector<Real> bond_coeffs_k_temp(_header._num_bonds);
  std::vector<Real> bond_coeffs_equilibrium_temp(_header._num_bonds);
  for (size_t i = 0; i < _header._num_bonds; i++)
  {
    bond_coeffs_k_temp[i] = _bond_coeffs_k[_bond_type[i]];
    bond_coeffs_equilibrium_temp[i] = _bond_coeffs_equilibrium[_bond_type[i]];
  }

  auto bond_coeffs_k = _para.GetFieldAsArrayHandle<Real>(field::bond_coeffs_k);
  vtkm::cont::ArrayCopy(vtkm::cont::make_ArrayHandle(bond_coeffs_k_temp), bond_coeffs_k);

  auto bond_coeffs_equilibrium = _para.GetFieldAsArrayHandle<Real>(field::bond_coeffs_equilibrium);
  vtkm::cont::ArrayCopy(vtkm::cont::make_ArrayHandle(bond_coeffs_equilibrium_temp),
                        bond_coeffs_equilibrium);

  if (_para.GetParameter<bool>(PARA_DIHEDRALS_FORCE))
  {
    _special_bonds = _para.GetParameter<std::vector<int>>(PARA_SPECIAL_BONDS);
    SetSpecialBonds();
  }
}

void ModelFileInitCondition::SetSpecialBonds() 
{
  std::vector<Real> special_weights;
  std::vector<Id> special_ids;
  std::vector<Id> special_offsets;

  for (const auto& atoms_id : _atoms_id)
  {
    //没有成键的部分
    if (_special_map.find(atoms_id) == _special_map.end())
    {
      special_weights.push_back(1.0);
      special_ids.push_back(atoms_id);
      special_offsets.push_back(1);
      continue;
    }

    //成键的部分
    Id offset = 0;
    auto link_0 = _special_map.equal_range(atoms_id);
    for (auto it0 = link_0.first; it0 != link_0.second; ++it0)
    {
      //一级连接
      int key_1 = it0->second;
      special_weights.push_back(_special_bonds[0]);
      special_ids.push_back(key_1);
      offset++;

      if (_special_map.find(key_1) == _special_map.end())
        continue;

      //二级链接
      auto link_1 = _special_map.equal_range(key_1);
      for (auto it1 = link_1.first; it1 != link_1.second; ++it1)
      {
        auto key_2 = it1->second;
        if (atoms_id == key_2)
          continue;

        special_weights.push_back(_special_bonds[1]);
        special_ids.push_back(key_2);
        offset++;

        if (_special_map.find(key_2) == _special_map.end())
          continue;

        //三级连接
        auto link_2 = _special_map.equal_range(key_2);
        for (auto it2 = link_2.first; it2 != link_2.second; ++it2)
        {
          auto key_3 = it2->second;
          if (key_1 == key_3)
            continue;

          special_weights.push_back(_special_bonds[2]);
          special_ids.push_back(key_3);
          offset++;
        }
      }
    }

    special_offsets.push_back(offset);
  }

  vtkm::Id offsetSize;
  vtkm::cont::ArrayHandle<vtkm::Id> offsetsArray = vtkm::cont::ConvertNumComponentsToOffsets(
    vtkm::cont::make_ArrayHandle(special_offsets), offsetSize);

  auto special_offsets_array = _para.GetFieldAsArrayHandle<Id>(field::special_offsets);
  vtkm::cont::ArrayCopy(offsetsArray, special_offsets_array);

  auto special_ids_array = _para.GetFieldAsArrayHandle<Id>(field::special_ids);
  vtkm::cont::ArrayCopy(vtkm::cont::make_ArrayHandle(special_ids), special_ids_array);

  auto special_weights_array = _para.GetFieldAsArrayHandle<Real>(field::special_weights);
  vtkm::cont::ArrayCopy(vtkm::cont::make_ArrayHandle(special_weights), special_weights_array);

  //vtkm::Id offsetSize;
  //vtkm::cont::ArrayHandle<vtkm::Id> offsetsArray = vtkm::cont::ConvertNumComponentsToOffsets(
  //  vtkm::cont::make_ArrayHandle(special_offsets), offsetSize);

  //auto groupSpecialIds = vtkm::cont::make_ArrayHandleGroupVecVariable(
  //  vtkm::cont::make_ArrayHandle(special_ids), offsetsArray);

  //auto groupSpecialWeights = vtkm::cont::make_ArrayHandleGroupVecVariable(
  //  vtkm::cont::make_ArrayHandle(special_weights), offsetsArray);

  //if (groupSpecialIds.GetNumberOfValues() != groupSpecialWeights.GetNumberOfValues())
  //  std::cout << "Error: Parse weight failed!" << std::endl;

  //for (int i = 0; i < groupSpecialIds.GetNumberOfValues(); ++i)
  //{
  //  std::cout << "id: " << i << "  ";
  //  auto item_ids = groupSpecialIds.ReadPortal().Get(i);
  //  auto item_weights = groupSpecialWeights.ReadPortal().Get(i);
  //  if (item_ids.GetNumberOfComponents() != item_weights.GetNumberOfComponents())
  //    std::cout << "Error: Parse weight failed!" << std::endl;
  //  for (int i = 0; i < item_ids.GetNumberOfComponents(); ++i)
  //  {
  //    std::cout << item_ids[i] << "-" << item_weights[i] << "  ";
  //  }

  //  std::cout << std::endl;
  //}


  //auto special_weights = _para.GetFieldAsArrayHandle<Real>(field::special_weights);
  //vtkm::cont::ArrayCopy(vtkm::cont::make_ArrayHandle(_special_weights), special_weights);
}

void ModelFileInitCondition::SetAngleField() 
{
  auto angle_atom_id = _para.GetFieldAsArrayHandle<Id>(field::angle_atom_id);
  vtkm::cont::ArrayCopy(vtkm::cont::make_ArrayHandle(_angle_atoms_id), angle_atom_id);

  auto angle_type_array = _para.GetFieldAsArrayHandle<Id>(field::angle_type);
  vtkm::cont::ArrayCopy(vtkm::cont::make_ArrayHandle(_angle_type), angle_type_array);

  std::vector<Real> angle_coeffs_k_temp(_header._num_angles);
  std::vector<Real> angle_coeffs_equilibrium_temp(_header._num_angles);
  for (size_t i = 0; i < _header._num_angles; i++)
  {
    angle_coeffs_k_temp[i] = _angle_coeffs_k[_angle_type[i]];
    angle_coeffs_equilibrium_temp[i] = _angle_coeffs_equilibrium[_angle_type[i]];
  }

  auto angle_coeffs_k = _para.GetFieldAsArrayHandle<Real>(field::angle_coeffs_k);
  vtkm::cont::ArrayCopy(vtkm::cont::make_ArrayHandle(angle_coeffs_k_temp), angle_coeffs_k);

  auto angle_coeffs_equilibrium =
    _para.GetFieldAsArrayHandle<Real>(field::angle_coeffs_equilibrium);
  vtkm::cont::ArrayCopy(vtkm::cont::make_ArrayHandle(angle_coeffs_equilibrium_temp),
                        angle_coeffs_equilibrium);
}

void ModelFileInitCondition::SetDihedralsField()
{
  auto dihedrals_atom_id = _para.GetFieldAsArrayHandle<Id>(field::dihedrals_atom_id);
  vtkm::cont::ArrayCopy(vtkm::cont::make_ArrayHandle(_dihedrals_atoms_id), dihedrals_atom_id);

  auto dihedrals_type_array = _para.GetFieldAsArrayHandle<Id>(field::dihedrals_type);
  vtkm::cont::ArrayCopy(vtkm::cont::make_ArrayHandle(_dihedrals_type), dihedrals_type_array);

  std::vector<Real> dihedrals_coeffs_k_temp(_header._num_dihedrals);
  std::vector<vtkm::IdComponent> dihedrals_coeffs_sign_temp(_header._num_dihedrals);
  std::vector<vtkm::IdComponent> dihedrals_coeffs_multiplicity_temp(_header._num_dihedrals);
  for (size_t i = 0; i < _header._num_dihedrals; i++)
  {
    dihedrals_coeffs_k_temp[i] = _dihedrals_coeffs_k[_dihedrals_type[i]];
    dihedrals_coeffs_sign_temp[i] = _dihedrals_coeffs_sign[_dihedrals_type[i]];
    dihedrals_coeffs_multiplicity_temp[i] = _dihedrals_coeffs_multiplicity[_dihedrals_type[i]];
  }

  auto dihedrals_coeffs_k = _para.GetFieldAsArrayHandle<Real>(field::dihedrals_coeffs_k);
  vtkm::cont::ArrayCopy(vtkm::cont::make_ArrayHandle(dihedrals_coeffs_k_temp), dihedrals_coeffs_k);

  auto dihedrals_coeffs_sign =
    _para.GetFieldAsArrayHandle<vtkm::IdComponent>(field::dihedrals_coeffs_sign);
  vtkm::cont::ArrayCopy(vtkm::cont::make_ArrayHandle(dihedrals_coeffs_sign_temp),
                        dihedrals_coeffs_sign);

  auto dihedrals_coeffs_multiplicity =
    _para.GetFieldAsArrayHandle<vtkm::IdComponent>(field::dihedrals_coeffs_multiplicity);
  vtkm::cont::ArrayCopy(vtkm::cont::make_ArrayHandle(dihedrals_coeffs_multiplicity_temp),
                        dihedrals_coeffs_multiplicity);
}

void ModelFileInitCondition::SetAtomsFieldEAM()
{
  // init Id
  auto atom_id = _para.GetFieldAsArrayHandle<Id>(field::atom_id);
  vtkm::cont::ArrayCopy(vtkm::cont::make_ArrayHandle(_atoms_id), atom_id);

  // init position
  ArrayHandle<Vec3f> position = _para.GetFieldAsArrayHandle<Vec3f>(field::position);
  vtkm::cont::ArrayCopy(vtkm::cont::make_ArrayHandle(_position), position);

  // init mass
  auto mass = _para.GetFieldAsArrayHandle<Real>(field::mass);
  mass.Allocate(_header._num_atoms);
  auto&& mass_protol = mass.WritePortal();
  for (int i = 0; i < _header._num_atoms; ++i)
  {
    mass_protol.Set(i, _mass[_atoms_type[i]]);
  }

  auto pts_type = _para.GetFieldAsArrayHandle<Id>(field::pts_type);
  vtkm::cont::ArrayCopy(vtkm::cont::make_ArrayHandle(_atoms_type), pts_type);

    // init position_flag
  auto position_flag = _para.GetFieldAsArrayHandle<Id3>(field::position_flag);
  position_flag.AllocateAndFill(_header._num_atoms, { 0, 0, 0 });
}
