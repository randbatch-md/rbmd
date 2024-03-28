﻿#pragma once
#include <string>

namespace field
{
	const std::string position = "position";
	const std::string velocity = "velocity";
	const std::string mass = "mass";
	const std::string charge = "charge";
	const std::string epsilon = "epsilon";
	const std::string sigma = "sigma";
	const std::string pts_type = "pts_type";
	const std::string molecule_id = "molecule_id";
    const std::string special_source_array = "special_source_array";
    const std::string special_offsets_array = "special_offsets_array";
	const std::string molecule_type = "molecule_type";
	const std::string atom_id = "atom_id";
	const std::string center_position = "center_position";
	const std::string target_position = "target_position";
	const std::string atom_id_center = "atom_id_center";
	const std::string atom_id_target = "atom_id_target";
	const std::string position_flag = "position_flag";
	const std::string bond_atom_id = "bond_atom_id";
	const std::string bond_type = "bond_type";
	const std::string bond_coeffs_k = "bond_coeffs_k";
	const std::string bond_coeffs_equilibrium = "bond_coeffs_equilibrium";
	const std::string angle_atom_id = "angle_atom_id";
	const std::string angle_type = "angle_type";
	const std::string angle_coeffs_k = "angle_coeffs_k";
	const std::string angle_coeffs_equilibrium = "angle_coeffs_equilibrium";
	const std::string dihedrals_atom_id = "dihedrals_atom_id";
	const std::string dihedrals_type = "dihedrals_type";
	const std::string signal_atoms_id = "signal_atoms_id";
    const std::string frho_spline = "frho_spline";
    const std::string rhor_spline = "rhor_spline";
    const std::string z2r_spline = "z2r_spline";
}

namespace gtest
{
	const std::string madelung = "madelung";
	const std::string ave_temp_t = "ave_temp_t";
	const std::string ave_potential_energy = "ave_potential_energy";
	const std::string ave_kin_energy = "ave_kin_energy";
    const std::string velocity_type = "velocity_type";
	const std::string test = "test";
}

const std::string PARA_CUTOFF = "cut_off";
const std::string PARA_RS = "rs";
const std::string PARA_RANDOM_NUM = "random_num";
const std::string PARA_VLENGTH = "vLength";
const std::string PARA_VOLUME = "volume";
const std::string PARA_BIN_NUMBER = "bin_number";
const std::string PARA_RANGE = "range";
const std::string PARA_UNIT = "unit";
const std::string PARA_UNIT_FACTOR = "unit_factor";
const std::string PARA_RHO = "rho";
const std::string PARA_RDF_RHO = "rdf_rho";
const std::string PARA_ALPHA = "alpha";
const std::string PARA_KMAX = "kmax";
const std::string PARA_BOND_ENERGY = "bond_energy";
const std::string PARA_ANGLE_ENERGY = "angle_energy";
const std::string PARA_TEMPT_SUM = "tempt_sum";
const std::string PARA_TEMPT = "temp_t";
const std::string PARA_RBE_P = "rbeP";
const std::string PARA_CENTER_TYPE = "center_type";
const std::string PARA_TARGET_TYPE = "target_type";
const std::string PARA_LOCATOR = "locator";
const std::string NTYPES = "num_atom_ntypes";
const std::string ATOM_STYLE = "atom_style";
const std::string EAM_PARA_CUTOFF = "eam_cut_off";
const std::string EAM_PARA_RHOMAX = "eam_rhomax";
const std::string EAM_PARA_NRHO = "eam_nrho";
const std::string EAM_PARA_DRHO = "eam_drho";
const std::string EAM_PARA_NR = "eam_nr";
const std::string EAM_PARA_DR = "eam_dr";