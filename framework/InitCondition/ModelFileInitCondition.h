//==================================================================================
//  RBMD 2.2.0 is developed for random batch molecular dynamics calculation.
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
//  The post-processing data produced by VASPKIT may be used in any publications 
//  provided that its use is explicitly acknowledged. A suitable reference for VASPKIT is:
//  [1] Gao W, Zhao T, Guo Y, et al.RBMD: A molecular dynamics package enabling to simulate 
//  10 million all - atom particles in a single graphics processing unit[J].arXiv preprint arXiv : 2407.09315, 2024.
// 
//  Contact Email : [support_wz@sciai.com.cn]
//==================================================================================

#pragma once
#include "MeshFreeFileInitCondition.h"
#include "FieldName.h"
#include "math/Utiles.h"
#include <map>
#include <vector>
struct ModelHeader
{
  vtkm::Id _num_atoms;
  vtkm::Id _num_bonds;
  vtkm::Id _num_angles;
  vtkm::Id _num_dihedrals;
  vtkm::Id _num_impropers;
  vtkm::Id _num_atoms_type;
  vtkm::Id _num_bound_type;
  vtkm::Id _num_angle_type;
  vtkm::Id _num_dihedrals_type;
};

enum MolecularType
{
	H20 = 0,
	NACL,
};

class ModelFileInitCondition : public MeshFreeFileInitCondition
{
public:
  ModelFileInitCondition(const Configuration& cfg);
  virtual ~ModelFileInitCondition() = default;

  void Execute() override;
  void ReadDataFile(std::ifstream&file);
  void InitField() override;
  void SetParameters() override;

private:
  void Header(std::ifstream& file);
  void Parser(std::ifstream& file);

  void Masses(std::ifstream& file, std::string line);
  void PairCoeffs(std::ifstream& file, std::string& line);
  void BondCoeffs(std::ifstream& file, std::string& line);
  void AngleCoeffs(std::ifstream& file, std::string& line);
  void DihedralsCoeffs(std::ifstream& file, std::string& line);
  void Group(std::string& line);
  void Atoms(std::ifstream& file, std::string& line);
  void Bonds(std::ifstream& file, std::string& line);
  void Angles(std::ifstream& file, std::string& line);
  void Dihedrals(std::ifstream& file, std::string& line);
  void VelocityRead(std::ifstream& file, std::string& line);
  void Velocity();

  void MolecularMapInsert(vtkm::Id& key, vtkm::Id& value);
  void AtomsMapInsert(vtkm::Id& key, vtkm::Id& value);
  void AtomstoMolecular(vtkm::Id& key, vtkm::Id& value);
  void SetMolecularVec();
  void SetMolecularGroup();
  void SetAtomsField();
  void SetBondField();
  void SetSpecialBonds();
  void SetAngleField();
  void SetAtomsFieldEAM();
  void SetDihedralsField();


private:
  ModelHeader _header;
  std::vector<Real> _mass;
  std::vector<Real> _sigma_atom;
  std::vector<Real> _eps_atom;
  std::vector<Real> _bond_coeffs_k;
  std::vector<Real> _bond_coeffs_equilibrium;
  std::vector<Real> _angle_coeffs_k;
  std::vector<Real> _angle_coeffs_equilibrium;
  std::vector<Id> _atoms_id;
  std::vector<Id> _signal_atoms_id;
  std::vector<Id> _molecules_id;
  std::vector<Id> _atoms_type;
  std::vector<Real> _charge;
  std::vector<Vec3f> _position;
  std::vector<Id> _bond_type;
  std::vector<Id> _bond_atoms_id;
  std::vector<Id> _angle_type;
  std::vector<Id> _angle_atoms_id;
  std::vector<Vec3f> _velocity;
  std::vector<Id> _dihedrals_type;
  std::vector<Id> _dihedrals_atoms_id;

  std::vector<Real> _dihedrals_coeffs_k;
  std::vector<vtkm::IdComponent> _dihedrals_coeffs_sign;
  std::vector<vtkm::IdComponent> _dihedrals_coeffs_multiplicity;

  std::map<Id, std::vector<Id>> _molecular_map;
  std::map<Id, Id> atom_to_molecular_map;
  std::map<Id, std::vector<Id>> _atoms_map;

  std::map<std::string, Id> group_id_init;
  std::map<Id, Id> _group_info;
  std::vector<Id> _molecular_type;
  std::vector<int> _special_bonds;
  std::multimap<Id, Id> _special_map;
};
