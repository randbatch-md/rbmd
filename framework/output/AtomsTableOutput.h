#pragma once
#include <vtkm/cont/ArrayHandleGroupVecVariable.h>
#include "ConsoleOutput.h"
#include "MeshFreeSystem.h"

class AtomsTableOutput : public ConsoleOutput
{
public:
  AtomsTableOutput(const Configuration& cfg);
  virtual ~AtomsTableOutput();

  void Init() override;
  void Execute() override;

private:
  void ComputePotentialEnergy();
  void Residual();
  void AddDataToTable();
  void WriteToFile();
  bool ShouldOutput();
  void ComputeEAMPotentialEnergy();

private:
  vtkm::cont::Timer _PotentialTimer;
  std::ofstream _file;
  IdComponent _Kmax;

  int _interval;
  bool _out_initial;
  bool _binary = false;
  bool _output_file;
  bool _compute;
  Real _spec_far_ele_potential_energy_avr;
  Real _cut_off;
  Real _Vlength;
  Real _volume;
  Real _alpha;
  Real _rho;
  Real _tempT_sum;
  Real _tempT;
  Real _potential_energy_avr_old;
  Real _residual;
  Real _potential_energy;
  Real _kinteic_energy;
  Real _total_energy;
  Real _lj_potential_energy_avr;
  Real _near_ele_potential_energy_avr;
  Real _far_ele_potential_energy_avr;
  Real _self_potential_energy_avr;
  Real _temperature;
  Real _pressure;
};
