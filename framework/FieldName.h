﻿//==================================================================================
//  RBMD 2.1.0 is developed for random batch molecular dynamics calculation.
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
	const std::string atom_pair_id = "atom_pair_id";
	const std::string atom_pair_position = "atom_pair_position";
    const std::string atoms_pair_type_offsets = "atoms_pair_type_offsets";
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
    const std::string special_weights = "special_weights";
    const std::string special_ids = "special_ids";
    const std::string special_offsets = "special_offsets";
    const std::string dihedrals_coeffs_k = "dihedrals_coeffs_k";
    const std::string dihedrals_coeffs_sign = "dihedrals_coeffs_sign";
    const std::string dihedrals_coeffs_multiplicity = "dihedrals_coeffs_multiplicity";
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
const std::string PARA_R_CORE = "r_core";
const std::string PARA_NEIGHBOR_SAMPLE_NUM = "neighbor_sample_num";
const std::string PARA_COULOMB_SAMPLE_NUM = "coulomb_sample_num";
const std::string PARA_VLENGTH = "vLength";
const std::string PARA_BOX = "box";
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
const std::string PARA_ATOMS_PAIR_TYPE = "atoms_pair_type";


const std::string PARA_LOCATOR = "locator";
const std::string NTYPES = "num_atom_ntypes";
const std::string ATOM_STYLE = "atom_style";
const std::string EAM_PARA_CUTOFF = "eam_cut_off";
const std::string EAM_PARA_RHOMAX = "eam_rhomax";
const std::string EAM_PARA_NRHO = "eam_nrho";
const std::string EAM_PARA_DRHO = "eam_drho";
const std::string EAM_PARA_NR = "eam_nr";
const std::string EAM_PARA_DR = "eam_dr";
const std::string PARA_DIHEDRAL_ENERGY = "dihedral_energy";

const std::string PARA_ATOM_MAP = "atom_map";
const std::string PARA_FORCE_FIELD_TYPE = "force_field_type";
const std::string PARA_POTENTIAL_FILE = "potential_file";
const std::string PARA_NEIGHBOR_TYPE = "neighbor_type";
const std::string PARA_FIX_SHAKE = "fix_shake";
const std::string PARA_COULOMB_TYPE = "coulomb_type";
const std::string PARA_ENSEMBLE = "ensemble";
const std::string PARA_TEMPERATURE = "temperature";
const std::string PARA_PRESSURE = "pressure";
const std::string PARA_TEMP_CTRL_TYPE = "temp_ctrl_type";
const std::string PARA_PRESS_CTRL_TYPE = "press_ctrl_type";
const std::string PARA_TIMESTEP = "timestep";
const std::string PARA_NUM_STEPS = "num_steps";
const std::string PARA_READ_FILE = "read_file";

const std::string PARA_FAR_FORCE = "far_force";
const std::string PARA_DIHEDRALS_FORCE = "dihedrals_force";
const std::string PARA_SPECIAL_BONDS = "special_bonds";

const std::string PARA_INIT_WAY = "init_way";
const std::string PARA_FILE_TYPE = "file_type";

const std::string PARA_FILE_DIHEDRALS = "file_dihedrals";
const std::string PARA_FILE_ANGLES = "file_angles";
