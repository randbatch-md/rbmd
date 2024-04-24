﻿#include "ModelFileInitCondition.h"
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
  _system.SetParameter(gtest::velocity_type, Get<std::string>("velocity_type"));
  _system.SetParameter(ATOM_STYLE, Get<std::string>("atom_style"));
}

void ModelFileInitCondition::Execute() 
{
  try
  {
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
  _system.SetParameter(PARA_VLENGTH, xLength);
  _system.SetParameter(PARA_VOLUME, xLength * yLength * zLength);
  _system.SetParameter(PARA_RHO, _header._num_atoms / (xLength * yLength * zLength));
  _system.SetParameter(NTYPES, _header._num_atoms_type);
  auto cut_off = _system.GetParameter<Real>(PARA_CUTOFF);
  auto bin_number = Id3{
    static_cast<int>(xLength / cut_off),
    static_cast<int>(yLength / cut_off),
    static_cast<int>(zLength / cut_off),
  };
  auto range = vtkm::Vec<vtkm::Range, 3>{
    { x_low, x_high },
    { y_low, y_high },
    { z_low, z_high }};
  _system.SetParameter(PARA_BIN_NUMBER, bin_number);
  _system.SetParameter(PARA_RANGE, range);
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
      else if (line.find("BondBond") != std::string::npos)
      {
        BondBondCoeffs(file, line);
      }
      else if (line.find("BondAngle") != std::string::npos)
      {
        BondAngleCoeffs(file, line);
      }
      else if (line.find("Bond Coeffs") != std::string::npos)
      {
        BondCoeffs(file, line);
      }
      else if (line.find("Angle Coeffs") != std::string::npos)
      {
        AngleCoeffs(file, line);
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
      }
      else if (line.find("Velocities") != std::string::npos)
      {
        VelocityRead(file, line);
      }
      else if (line.find("Dihedrals") != std::string::npos)
      {
        Dihedrals(file, line);
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

  auto unit = _system.GetParameter<std::string>(PARA_UNIT);
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
  auto epsilon = _system.GetFieldAsArrayHandle<Real>(field::epsilon);
  vtkm::cont::ArrayCopy(vtkm::cont::make_ArrayHandle(_eps_atom), epsilon);

  auto sigma = _system.GetFieldAsArrayHandle<Real>(field::sigma);
  vtkm::cont::ArrayCopy(vtkm::cont::make_ArrayHandle(_sigma_atom), sigma);
}

void ModelFileInitCondition::BondCoeffs(std::ifstream& file, std::string& line)
{
  std::vector<Real> bond_id_bondcoeffs(_header._num_bound_type);
  _bond_coeffs_k.resize(_header._num_bound_type);
  _bond_coeffs_k3.resize(_header._num_bound_type);
  _bond_coeffs_k4.resize(_header._num_bound_type);
  _bond_coeffs_equilibrium.resize(_header._num_bound_type);

  if (1)
  {
    for (int i = 0; i < _header._num_bound_type && getline(file, line); i)
    {
      if (!line.empty() && line != "\r")
      {
        std::istringstream iss(line);
        iss >> bond_id_bondcoeffs[i] >> _bond_coeffs_equilibrium[i] >> _bond_coeffs_k[i] >>
          _bond_coeffs_k3[i] >> _bond_coeffs_k4[i];
        i++;
      }
    }
  }
  else
  {
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
}

void ModelFileInitCondition::AngleCoeffs(std::ifstream& file, std::string& line)
{
  std::vector<Real> angle_id_anglecoeffs(_header._num_angle_type);
  _angle_coeffs_k.resize(_header._num_angle_type);
  _angle_coeffs_k3.resize(_header._num_angle_type);
  _angle_coeffs_k4.resize(_header._num_angle_type);
  _angle_coeffs_equilibrium.resize(_header._num_angle_type);

  if (1)
  {
    for (int i = 0; i < _header._num_angle_type && getline(file, line); i)
    {
      if (!line.empty() && line != "\r")
      {
        std::istringstream iss(line);
        iss >> angle_id_anglecoeffs[i] >> _angle_coeffs_equilibrium[i] >> _angle_coeffs_k[i] >>
          _angle_coeffs_k3[i] >> _angle_coeffs_k4[i];
        i++;
      }
    }
  }
  else
  {
    for (int i = 0; i < _header._num_angle_type && getline(file, line); i)
    {
      if (!line.empty() && line != "\r")
      {
        std::istringstream iss(line);
        iss >> angle_id_anglecoeffs[i] >> _angle_coeffs_k[i] >> _angle_coeffs_equilibrium[i];
        i++;
      }
    }
  }
}

void ModelFileInitCondition::BondBondCoeffs(std::ifstream& file, std::string& line)
{
  std::vector<Real> bond_id_bbcoeffs(_header._num_angle_type);
  _bb_m.resize(_header._num_angle_type);
  _bb_r1.resize(_header._num_angle_type);
  _bb_r2.resize(_header._num_angle_type);

  for (int i = 0; i < _header._num_angle_type && getline(file, line); i)
  {
    if (!line.empty() && line != "\r")
    {
      std::istringstream iss(line);
      iss >> bond_id_bbcoeffs[i] >> _bb_m[i] >> _bb_r1[i] >> _bb_r2[i];
      i++;
    }
  }
}

void ModelFileInitCondition::BondAngleCoeffs(std::ifstream& file, std::string& line)
{
  std::vector<Real> bond_id_bacoeffs(_header._num_angle_type);
  _ba_n1.resize(_header._num_angle_type);
  _ba_n2.resize(_header._num_angle_type);
  _ba_r1.resize(_header._num_angle_type);
  _ba_r2.resize(_header._num_angle_type);

  for (int i = 0; i < _header._num_angle_type && getline(file, line); i)
  {
    if (!line.empty() && line != "\r")
    {
      std::istringstream iss(line);
      iss >> bond_id_bacoeffs[i] >> _ba_n1[i] >> _ba_n2[i] >> _ba_r1[i] >> _ba_r2[i];
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
  auto atom_style = _system.GetParameter<std::string>(ATOM_STYLE);
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
      _bond_atoms_id[index++] = temp_id1 - 1;
      _bond_atoms_id[index++] = temp_id2 - 1;
      i++;

      bond_atoms_set.insert(temp_id1 - 1);
      bond_atoms_set.insert(temp_id2 - 1);
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
  }

  auto dihedrals_atoms_id = _system.GetFieldAsArrayHandle<Id>(field::dihedrals_atom_id);
  vtkm::cont::ArrayCopy(vtkm::cont::make_ArrayHandle(_dihedrals_atoms_id), dihedrals_atoms_id);

  auto dihedrals_type_array = _system.GetFieldAsArrayHandle<Id>(field::dihedrals_type);
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

  auto velocity = _system.GetFieldAsArrayHandle<Vec3f>(field::velocity);
  vtkm::cont::ArrayCopy(vtkm::cont::make_ArrayHandle(_velocity), velocity);
}

void ModelFileInitCondition::Velocity()
{
  if (_velocity.size() == NULL)
  {
    _velocity.resize(_header._num_atoms);
    auto velocity_type = _system.GetParameter<std::string>(gtest::velocity_type);

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
        auto unit = _system.GetParameter<std::string>(PARA_UNIT);
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
    auto velocity = _system.GetFieldAsArrayHandle<Vec3f>(field::velocity);
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

  auto source_array = _system.GetFieldAsArrayHandle<Id>(field::special_source_array);
  vtkm::cont::ArrayCopy(sourceArray, source_array);
  auto offsets_array = _system.GetFieldAsArrayHandle<Id>(field::special_offsets_array);
  vtkm::cont::ArrayCopy(offsetsArray, offsets_array);

  //
  auto molecule_id = _system.GetFieldAsArrayHandle<Id>(field::molecule_id);
  vtkm::cont::ArrayCopy(vtkm::cont::make_ArrayHandle(_molecules_id), molecule_id);
}

void ModelFileInitCondition::SetAtomsField() 
{
  // init Id
  auto atom_id = _system.GetFieldAsArrayHandle<Id>(field::atom_id);
  vtkm::cont::ArrayCopy(vtkm::cont::make_ArrayHandle(_atoms_id), atom_id);

  // init position
  ArrayHandle<Vec3f> position = _system.GetFieldAsArrayHandle<Vec3f>(field::position);
  vtkm::cont::ArrayCopy(vtkm::cont::make_ArrayHandle(_position), position);

  // init mass
  auto mass = _system.GetFieldAsArrayHandle<Real>(field::mass);
  mass.Allocate(_header._num_atoms);
  auto&& mass_protol = mass.WritePortal();
  for (int i = 0; i < _header._num_atoms; ++i)
  {
    mass_protol.Set(i, _mass[_atoms_type[i]]);
  }

  // init charge
  auto charge = _system.GetFieldAsArrayHandle<Real>(field::charge);
  vtkm::cont::ArrayCopy(vtkm::cont::make_ArrayHandle(_charge), charge);

  auto pts_type = _system.GetFieldAsArrayHandle<Id>(field::pts_type);
  vtkm::cont::ArrayCopy(vtkm::cont::make_ArrayHandle(_atoms_type), pts_type);

  // Init Center_Target Position
  auto center_type = _system.GetParameter<IdComponent>(PARA_CENTER_TYPE);
  auto atoms_center = _atoms_map[center_type];
  auto atoms_id_center = _system.GetFieldAsArrayHandle<Id>(field::atom_id_center);
  vtkm::cont::ArrayCopy(vtkm::cont::make_ArrayHandle(atoms_center), atoms_id_center);

  auto target_type = _system.GetParameter<IdComponent>(PARA_TARGET_TYPE);
  auto atoms_target = _atoms_map[target_type];
  auto atoms_id_target = _system.GetFieldAsArrayHandle<Id>(field::atom_id_target);
  vtkm::cont::ArrayCopy(vtkm::cont::make_ArrayHandle(atoms_target), atoms_id_target);

  // init position_flag
  auto position_flag = _system.GetFieldAsArrayHandle<Id3>(field::position_flag);
  position_flag.AllocateAndFill(_header._num_atoms, { 0, 0, 0 });
}

void ModelFileInitCondition::SetBondField() 
{
  auto signal_atoms_id = _system.GetFieldAsArrayHandle<Id>(field::signal_atoms_id);
  vtkm::cont::ArrayCopy(vtkm::cont::make_ArrayHandle(_signal_atoms_id), signal_atoms_id);

  auto bond_atom_id = _system.GetFieldAsArrayHandle<Id>(field::bond_atom_id);
  vtkm::cont::ArrayCopy(vtkm::cont::make_ArrayHandle(_bond_atoms_id), bond_atom_id);

  auto bond_type_handle = _system.GetFieldAsArrayHandle<Id>(field::bond_type);
  vtkm::cont::ArrayCopy(vtkm::cont::make_ArrayHandle(_bond_type), bond_type_handle);

  std::vector<Real> bond_coeffs_k_temp(_header._num_bonds);
  std::vector<Real> bond_coeffs_k3_temp(_header._num_bonds);
  std::vector<Real> bond_coeffs_k4_temp(_header._num_bonds);
  std::vector<Real> bond_coeffs_equilibrium_temp(_header._num_bonds);
  if (1)
  {
    for (size_t i = 0; i < _header._num_bonds; i++)
    {
      bond_coeffs_k_temp[i] = _bond_coeffs_k[_bond_type[i]];
      bond_coeffs_k3_temp[i] = _bond_coeffs_k3[_bond_type[i]];
      bond_coeffs_k4_temp[i] = _bond_coeffs_k4[_bond_type[i]];
      bond_coeffs_equilibrium_temp[i] = _bond_coeffs_equilibrium[_bond_type[i]];
    }

    auto bond_coeffs_k3 = _system.GetFieldAsArrayHandle<Real>(field::bond_coeffs_k3);
    vtkm::cont::ArrayCopy(vtkm::cont::make_ArrayHandle(bond_coeffs_k3_temp), bond_coeffs_k3);

    auto bond_coeffs_k4 = _system.GetFieldAsArrayHandle<Real>(field::bond_coeffs_k4);
    vtkm::cont::ArrayCopy(vtkm::cont::make_ArrayHandle(bond_coeffs_k4_temp), bond_coeffs_k4);
  }
  else
  {
    for (size_t i = 0; i < _header._num_bonds; i++)
    {
      bond_coeffs_k_temp[i] = _bond_coeffs_k[_bond_type[i]];
      bond_coeffs_equilibrium_temp[i] = _bond_coeffs_equilibrium[_bond_type[i]];
    }
  }

  auto bond_coeffs_k = _system.GetFieldAsArrayHandle<Real>(field::bond_coeffs_k);
  vtkm::cont::ArrayCopy(vtkm::cont::make_ArrayHandle(bond_coeffs_k_temp), bond_coeffs_k);

  auto bond_coeffs_equilibrium = _system.GetFieldAsArrayHandle<Real>(field::bond_coeffs_equilibrium);
  vtkm::cont::ArrayCopy(vtkm::cont::make_ArrayHandle(bond_coeffs_equilibrium_temp),
                        bond_coeffs_equilibrium);
}

void ModelFileInitCondition::SetAngleField() 
{
  auto angle_atom_id = _system.GetFieldAsArrayHandle<Id>(field::angle_atom_id);
  vtkm::cont::ArrayCopy(vtkm::cont::make_ArrayHandle(_angle_atoms_id), angle_atom_id);

  auto angle_type_array = _system.GetFieldAsArrayHandle<Id>(field::angle_type);
  vtkm::cont::ArrayCopy(vtkm::cont::make_ArrayHandle(_angle_type), angle_type_array);

  std::vector<Real> angle_coeffs_k_temp(_header._num_angles);
  std::vector<Real> angle_coeffs_k3_temp(_header._num_angles);
  std::vector<Real> angle_coeffs_k4_temp(_header._num_angles);
  std::vector<Real> bb_m_temp(_header._num_angles);
  std::vector<Real> bb_r1_temp(_header._num_angles);
  std::vector<Real> bb_r2_temp(_header._num_angles);
  std::vector<Real> ba_n1_temp(_header._num_angles);
  std::vector<Real> ba_n2_temp(_header._num_angles);
  std::vector<Real> ba_r1_temp(_header._num_angles);
  std::vector<Real> ba_r2_temp(_header._num_angles);
  std::vector<Real> angle_coeffs_equilibrium_temp(_header._num_angles);

  if (1)
  {
    for (size_t i = 0; i < _header._num_angles; i++)
    {
      angle_coeffs_k_temp[i] = _angle_coeffs_k[_angle_type[i]];
      angle_coeffs_k3_temp[i] = _angle_coeffs_k3[_angle_type[i]];
      angle_coeffs_k4_temp[i] = _angle_coeffs_k4[_angle_type[i]];
      angle_coeffs_equilibrium_temp[i] = _angle_coeffs_equilibrium[_angle_type[i]];
      bb_m_temp[i] = _bb_m[_angle_type[i]];
      bb_r1_temp[i] = _bb_r1[_angle_type[i]];
      bb_r2_temp[i] = _bb_r2[_angle_type[i]];
      ba_n1_temp[i] = _ba_n1[_angle_type[i]];
      ba_n2_temp[i] = _ba_n2[_angle_type[i]];
      ba_r1_temp[i] = _ba_r1[_angle_type[i]];
      ba_r2_temp[i] = _ba_r2[_angle_type[i]];
      auto angle_coeffs_k3 = _system.GetFieldAsArrayHandle<Real>(field::angle_coeffs_k3);
      auto angle_coeffs_k4 = _system.GetFieldAsArrayHandle<Real>(field::angle_coeffs_k4);
      auto bb_m = _system.GetFieldAsArrayHandle<Real>(field::bb_m);
      auto bb_r1 = _system.GetFieldAsArrayHandle<Real>(field::bb_r1);
      auto bb_r2 = _system.GetFieldAsArrayHandle<Real>(field::bb_r2);
      auto ba_n1 = _system.GetFieldAsArrayHandle<Real>(field::ba_n1);
      auto ba_n2 = _system.GetFieldAsArrayHandle<Real>(field::ba_n2);
      auto ba_r1 = _system.GetFieldAsArrayHandle<Real>(field::ba_r1);
      auto ba_r2 = _system.GetFieldAsArrayHandle<Real>(field::ba_r2);
      vtkm::cont::ArrayCopy(vtkm::cont::make_ArrayHandle(angle_coeffs_k3_temp), angle_coeffs_k3);
      vtkm::cont::ArrayCopy(vtkm::cont::make_ArrayHandle(angle_coeffs_k4_temp), angle_coeffs_k4);
      vtkm::cont::ArrayCopy(vtkm::cont::make_ArrayHandle(bb_m_temp), bb_m);
      vtkm::cont::ArrayCopy(vtkm::cont::make_ArrayHandle(bb_r1_temp), bb_r1);
      vtkm::cont::ArrayCopy(vtkm::cont::make_ArrayHandle(bb_r2_temp), bb_r2);
      vtkm::cont::ArrayCopy(vtkm::cont::make_ArrayHandle(ba_n1_temp), ba_n1);
      vtkm::cont::ArrayCopy(vtkm::cont::make_ArrayHandle(ba_n2_temp), ba_n2);
      vtkm::cont::ArrayCopy(vtkm::cont::make_ArrayHandle(ba_r1_temp), ba_r1);
      vtkm::cont::ArrayCopy(vtkm::cont::make_ArrayHandle(ba_r2_temp), ba_r2);
    }
  }
  else
  {
    for (size_t i = 0; i < _header._num_angles; i++)
    {
      angle_coeffs_k_temp[i] = _angle_coeffs_k[_angle_type[i]];
      angle_coeffs_equilibrium_temp[i] = _angle_coeffs_equilibrium[_angle_type[i]];
    }
  }

  auto angle_coeffs_k = _system.GetFieldAsArrayHandle<Real>(field::angle_coeffs_k);
  vtkm::cont::ArrayCopy(vtkm::cont::make_ArrayHandle(angle_coeffs_k_temp), angle_coeffs_k);

  auto angle_coeffs_equilibrium =
    _system.GetFieldAsArrayHandle<Real>(field::angle_coeffs_equilibrium);
  vtkm::cont::ArrayCopy(vtkm::cont::make_ArrayHandle(angle_coeffs_equilibrium_temp),
                        angle_coeffs_equilibrium);
}

void ModelFileInitCondition::SetAtomsFieldEAM()
{
  // init Id
  auto atom_id = _system.GetFieldAsArrayHandle<Id>(field::atom_id);
  vtkm::cont::ArrayCopy(vtkm::cont::make_ArrayHandle(_atoms_id), atom_id);

  // init position
  ArrayHandle<Vec3f> position = _system.GetFieldAsArrayHandle<Vec3f>(field::position);
  vtkm::cont::ArrayCopy(vtkm::cont::make_ArrayHandle(_position), position);

  // init mass
  auto mass = _system.GetFieldAsArrayHandle<Real>(field::mass);
  mass.Allocate(_header._num_atoms);
  auto&& mass_protol = mass.WritePortal();
  for (int i = 0; i < _header._num_atoms; ++i)
  {
    mass_protol.Set(i, _mass[_atoms_type[i]]);
  }

  auto pts_type = _system.GetFieldAsArrayHandle<Id>(field::pts_type);
  vtkm::cont::ArrayCopy(vtkm::cont::make_ArrayHandle(_atoms_type), pts_type);

    // init position_flag
  auto position_flag = _system.GetFieldAsArrayHandle<Id3>(field::position_flag);
  position_flag.AllocateAndFill(_header._num_atoms, { 0, 0, 0 });
}
