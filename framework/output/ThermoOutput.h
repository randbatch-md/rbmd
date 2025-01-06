//==================================================================================
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
#include <vtkm/cont/ArrayHandleGroupVecVariable.h>
#include "ConsoleOutput.h"

enum STATUS
{
  POTENTIAL_ENERGY,
  TOTAL_ENERGY,
  KIN_ENERGY,
  TEMPERATURE
};

class ThermoOutput : public ConsoleOutput
{
public:
  ThermoOutput(const Configuration& cfg);
  virtual ~ThermoOutput();

  void Init() override;
  void Execute() override;

private:
  void ComputePotentialEnergy();
  void Residual();
  void AddDataToTable();
  void WriteToFile();
  bool ShouldOutput();
  void StatisticalStatus();
  void ComputeEAMPotentialEnergy();

  void PostData();
  void PostExecute();

private:
  std::ofstream _file;
  std::ofstream _system_state;

  int _interval;
  Real _spec_far_ele_potential_energy_avr;
  Real _cut_off;

  Vec3f _box;
  Real _volume;

  Real _bond_energy;
  Real _angle_energy;
  Real _rho;
  Real _tempT_sum;
  Real _tempT;
  Real _potential_energy_avr_old;
  Real _residual;
  Real _potential_energy_avr; //(non bond)
  Real _potential_energy;
  Real _kinteic_energy;
  Real _total_energy;
  Real _lj_potential_energy_avr;
  Real _near_ele_potential_energy_avr;
  Real _far_ele_potential_energy_avr;
  Real _self_potential_energy_avr;
  Real _temperature;

  Real _dihedrals_energy;

  std::map<STATUS, std::vector<Real>> _status;
};
