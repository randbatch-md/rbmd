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

private:
  void Header(std::ifstream& file);
  void Parser(std::ifstream& file);

  void Masses(std::ifstream& file, std::string line);
  void PairCoeffs(std::ifstream& file, std::string& line);
  void BondCoeffs(std::ifstream& file, std::string& line);
  void AngleCoeffs(std::ifstream& file, std::string& line);
  void BondBondCoeffs(std::ifstream& file, std::string& line);
  void BondAngleCoeffs(std::ifstream& file, std::string& line);
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
  void SetAngleField();
  void SetAtomsFieldEAM();

private:
  ModelHeader _header;
  std::vector<Real> _mass;
  std::vector<Real> _sigma_atom;
  std::vector<Real> _eps_atom;
  std::vector<Real> _bond_coeffs_k;
  std::vector<Real> _bond_coeffs_k3;
  std::vector<Real> _bond_coeffs_k4;
  std::vector<Real> _bond_coeffs_equilibrium;
  std::vector<Real> _angle_coeffs_k;
  std::vector<Real> _angle_coeffs_k3;
  std::vector<Real> _angle_coeffs_k4;
  std::vector<Real> _angle_coeffs_equilibrium;
  std::vector<Real> _bb_m;
  std::vector<Real> _bb_r1;
  std::vector<Real> _bb_r2;
  std::vector<Real> _ba_n1;
  std::vector<Real> _ba_n2;
  std::vector<Real> _ba_r1;
  std::vector<Real> _ba_r2;
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

  std::map<Id, std::vector<Id>> _molecular_map;
  std::map<Id, Id> atom_to_molecular_map;
  std::map<Id, std::vector<Id>> _atoms_map;

  std::map<std::string, Id> group_id_init;
  std::map<Id, Id> _group_info;
  std::vector<Id> _molecular_type;

};
