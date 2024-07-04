#pragma once
#include <vtkm/cont/ArrayHandleGroupVecVariable.h>
#include "ConsoleOutput.h"
#include "MeshFreeSystem.h"

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
  void SpecialFarCoulEnergy();

private:
  std::ofstream _file;
  std::ofstream _system_state;
  IdComponent _Kmax;

  int _interval;
  Real _spec_far_ele_potential_energy_avr;
  Real _cut_off;
  Real _Vlength;
  Vec3f _box;
  Real _volume;
  Real _alpha;
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
