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

#include "Executioner.h"
#include "ThermoOutput.h"
#include <vtkm/cont/Algorithm.h>
#include "Application.h"
#include "ConsoleOutput.h"
#include "math/Math.h"
#include "output/worklet/OutPutWorklet.h"
#include "run/worklet/RunWorklet.h"
#include "FieldName.h"
#include <ctime>
#include "forceFunction/ContForceFunction.h"
#include "UnitFactor.h"

const std::string HEADER_POTENTIAL_ENERGY_NAME = "Potential_Energy";
const std::string HEADER_RESIDUAL_NAME = "Residual";
const std::string HEADER_KINETIC_ENERGY_NAME = "Kinetic_Energy";
const std::string HEADER_NON_BOND_ENERGY_NAME = "Non_Bond_Energy";
const std::string HEADER_TOTAL_ENERGY_NAME = "Total_Energy";
const std::string HEADER_BOND_ENERGY_NAME = "Bond_Energy";
const std::string HEADER_ANGLE_ENERGY_NAME = "Angle_Energy";
const std::string HEADER_KBT_NAME = "KBT";
const std::string HEADER_CUMULATIVE_TIME_NAME = "Cumulative_Time";
const std::string HEADER_DIHEDRAL_ENERGY_NAME = "Dihedral_Energy";

ThermoOutput::ThermoOutput(const Configuration& cfg)
  : ConsoleOutput(cfg)
  , _file("energy.rbmd")
  , _interval(Get<int>("interval", 1))
  //, _system_state("SystemState.csv")
  , _cut_off(0.0)
  , _volume(0.0)
  , _bond_energy(0.0)
  , _angle_energy(0.0)
  , _dihedrals_energy(0.0)
  , _rho(0.0)
  , _tempT_sum(0.0)
  , _tempT(0.0)
{

}

ThermoOutput::~ThermoOutput()
{
  _file.close();
  _system_state.close();
}

void ThermoOutput::Init()
{
  ConsoleOutput::Init();
  AddHeader(HEADER_POTENTIAL_ENERGY_NAME);
  AddHeader(HEADER_RESIDUAL_NAME);
  AddHeader(HEADER_KBT_NAME);
  AddHeader(HEADER_KINETIC_ENERGY_NAME);
  AddHeader(HEADER_TOTAL_ENERGY_NAME);
  AddHeader(HEADER_CUMULATIVE_TIME_NAME);
  
  if (_para.GetParameter<std::string>(PARA_INIT_WAY) != "inbuild")
  {
    AddHeader(HEADER_NON_BOND_ENERGY_NAME);
    AddHeader(HEADER_BOND_ENERGY_NAME);
    AddHeader(HEADER_ANGLE_ENERGY_NAME);
    AddHeader(HEADER_DIHEDRAL_ENERGY_NAME);
  }

  _cut_off = _para.GetParameter<Real>(PARA_CUTOFF);
  _volume = _para.GetParameter<Real>(PARA_VOLUME);
  _box = _para.GetParameter<Vec3f>(PARA_BOX);
  _rho = _para.GetParameter<Real>(PARA_RHO);
}

void ThermoOutput::Execute()
{
  if (_para.HaveParameter(PARA_BOND_ENERGY))
  {
    _bond_energy = _para.GetParameter<Real>(PARA_BOND_ENERGY);
  }
  if (_para.HaveParameter(PARA_ANGLE_ENERGY))
  {
    _angle_energy = _para.GetParameter<Real>(PARA_ANGLE_ENERGY);
  }
  if (_para.HaveParameter(PARA_DIHEDRAL_ENERGY))
  {
    _dihedrals_energy = _para.GetParameter<Real>(PARA_DIHEDRAL_ENERGY);
  }
  _tempT_sum = _para.GetParameter<Real>(PARA_TEMPT_SUM);
  _tempT = _para.GetParameter<Real>(PARA_TEMPT);

  if (_para.GetParameter<std::string>(PARA_FILE_TYPE) == "EAM")
  {
    ComputeEAMPotentialEnergy();
  }
  else
  {
    ComputePotentialEnergy();
  }

  Residual();

  AddDataToTable();

  WriteToFile();

  //StatisticalStatus();
  ConsoleOutput::Execute();
  //PostData();
}

bool ThermoOutput::ShouldOutput()
{
  auto current_step = _executioner->CurrentStep();
  if (current_step == _executioner->NumStep())
    return true;
  else
    return current_step % _interval == 0;
}

void ThermoOutput::ComputePotentialEnergy()
{ 
  auto position = _para.GetFieldAsArrayHandle<Vec3f>(field::position);
  auto N = position.GetNumberOfValues();

  _lj_potential_energy_avr = 0.0;
  _near_ele_potential_energy_avr = 0.0;
  _self_potential_energy_avr = 0.0;
  _far_ele_potential_energy_avr = 0.0;
  

  auto force_field = _para.GetParameter<std::string>(PARA_FORCE_FIELD_TYPE);
  if ("LJ/CUT" == force_field)
  {
      _lj_potential_energy_avr = _para.GetParameter<Real>(PARA_LJ_ENERGY);
      _potential_energy = _lj_potential_energy_avr;
  }
  else if ("LJ/CUT/COUL/LONG" == force_field)
  {
      _lj_potential_energy_avr = _para.GetParameter<Real>(PARA_LJ_ENERGY);
      _near_ele_potential_energy_avr = _para.GetParameter<Real>(PARA_COUL_ENERGY);
      _far_ele_potential_energy_avr = _para.GetParameter<Real>(PARA_EWALD_LONG_ENERGY);
      _potential_energy = _lj_potential_energy_avr + _far_ele_potential_energy_avr + _near_ele_potential_energy_avr;
  }
  else if ("CVFF" == force_field)
  {
     _lj_potential_energy_avr = _para.GetParameter<Real>(PARA_LJ_ENERGY);
     _near_ele_potential_energy_avr = _para.GetParameter<Real>(PARA_COUL_ENERGY);
     _far_ele_potential_energy_avr = _para.GetParameter<Real>(PARA_EWALD_LONG_ENERGY);

     _potential_energy_avr = _lj_potential_energy_avr + _far_ele_potential_energy_avr + _near_ele_potential_energy_avr;
     _potential_energy = _potential_energy_avr + _bond_energy + _angle_energy + _dihedrals_energy; //
  }

  //TOTAL:

  _kinteic_energy = _tempT_sum / N * 0.5;

  _total_energy = _potential_energy + _kinteic_energy;

  _temperature = _tempT;
 
}

void ThermoOutput::Residual()
{
  if (_executioner->CurrentStep() <= 0)
  {
    _residual = 0;
  }
  if (_executioner->CurrentStep() > 0)
  {
    _residual = _potential_energy_avr - _potential_energy_avr_old;
  }
  _potential_energy_avr_old = _potential_energy_avr;
}

void ThermoOutput::AddDataToTable()
{
  _output_screen = true;
  if (_output_screen)
  {
    AddData(HEADER_POTENTIAL_ENERGY_NAME, _potential_energy);
    AddData(HEADER_RESIDUAL_NAME, _residual);
    AddData(HEADER_TOTAL_ENERGY_NAME, _total_energy);
    AddData(HEADER_KINETIC_ENERGY_NAME, _kinteic_energy);
    AddData(HEADER_KBT_NAME, _tempT);
    AddData(HEADER_CUMULATIVE_TIME_NAME, _cumulative_time + _timer.GetElapsedTime());
    if (_para.GetParameter<std::string>(PARA_INIT_WAY) != "inbuild")
    {
      AddData(HEADER_NON_BOND_ENERGY_NAME, _potential_energy_avr);
      AddData(HEADER_BOND_ENERGY_NAME, _bond_energy);
      AddData(HEADER_ANGLE_ENERGY_NAME, _angle_energy);
      AddData(HEADER_DIHEDRAL_ENERGY_NAME, _dihedrals_energy); 
    }
    
  }
}

void ThermoOutput::WriteToFile()
{
  if (_executioner->CurrentStep() == 0)
  {
    _file << "Step"
          << " "
          << "Time"
          << " "
          << "VanderWaalsEnergy" // lj energy
          << " "
          << "NearCoulombicEnergy" // add near energy
          << " "                   //
          << "FarCoulombicEnergy"  // add far energy
          << " "                   //
          << "Residual"
          << " "
          << "KinticEnergy"
          << " "
          << "PotentialEnergy"
          << " "
          << "TotalEnergy"
          << " "
          << "CumulativeTime";
    if (_para.GetParameter<std::string>(PARA_INIT_WAY) != "inbuild")
    {
      _file << " "
            << "NonBonEnergy"
            << " "
            << "BondEnergy"
            << " "
            << "AngleEnergy"
            << " "
            << "DihedralEnergy"
            << " "
            << "SpecialEnergy" ;
    }
    _file << std::endl;
  }
  if (ShouldOutput())
  {
    try
    {
        _file << _executioner->CurrentStep() << " " 
            << _time<< " "
            << _lj_potential_energy_avr << " "
            << _near_ele_potential_energy_avr << " "
            << _far_ele_potential_energy_avr << " " 
            << _residual << " "
            << _kinteic_energy << " "
            << _potential_energy << " "
            << _total_energy << " "
            << _cumulative_time + _timer.GetElapsedTime() << " ";
      if (_para.GetParameter<std::string>(PARA_INIT_WAY) != "inbuild")
        {
        _file << _potential_energy_avr << " "
              << _bond_energy << " "
              << _angle_energy << " " << _dihedrals_energy << " "
              << _spec_far_ele_potential_energy_avr;
        }
        _file << std::endl;
      
    }
    catch (const std::exception& e)
    {
      _file.close();
      console::Error(e.what());
    }
  }
}

void ThermoOutput::ComputeEAMPotentialEnergy()
{
  auto atoms_id = _para.GetFieldAsArrayHandle<Id>(field::atom_id);
  auto N = atoms_id.GetNumberOfValues();
  //
  auto cut_off = _para.GetParameter<Real>(EAM_PARA_CUTOFF);
  auto Vlength = _para.GetParameter<Real>(PARA_VLENGTH);
  auto box = _para.GetParameter<Vec3f>(PARA_BOX);

  auto rhor_spline = _para.GetFieldAsArrayHandle<Vec7f>(field::rhor_spline);
  auto frho_spline = _para.GetFieldAsArrayHandle<Vec7f>(field::frho_spline);
  auto z2r_spline = _para.GetFieldAsArrayHandle<Vec7f>(field::z2r_spline);

  //
  ContForceFunction force_function;
  SetForceFunction(force_function);

  ContTopology topology;
  SetTopology(topology);

  ContPointLocator locator;
  SetLocator(locator);

  //1:compute EAM_rho;
  ArrayHandle<Real> EAM_rho;
  OutPut::EAM_rho(cut_off, box, atoms_id, rhor_spline, locator, topology, force_function, EAM_rho);

  //2: embedding_energy_atom
  ArrayHandle<Real> embedding_energy_atom;
  OutPut::EAM_EmbeddingEnergy(
    atoms_id, EAM_rho, frho_spline, locator, topology, force_function, embedding_energy_atom);

  auto embedding_energy_total = vtkm::cont::Algorithm::Reduce(
    embedding_energy_atom, vtkm::TypeTraits<Real>::ZeroInitialization());

  //3: pair_energy_atom
  ArrayHandle<Real> pair_energy_atom;
  OutPut::EAM_PairEnergy(
    cut_off, box, atoms_id, z2r_spline, locator, topology, force_function, pair_energy_atom);

  auto pair_energy_atom_total =
    vtkm::cont::Algorithm::Reduce(pair_energy_atom, vtkm::TypeTraits<Real>::ZeroInitialization());

  //4:_potential_energy = embedding_energy_total +  pair_energy_atom_total
  _potential_energy = (embedding_energy_total + pair_energy_atom_total) / N; // ;

  //5._kinteic_energy
  _kinteic_energy = _tempT_sum / N * 0.5;

  _total_energy = _potential_energy + _kinteic_energy;

  //6:temperature
  _temperature = _tempT;

  //7.all_potential_atom

  //ArrayHandle<Real> all_potential_atom;
  //OutPut::SumEAM_Pe_Atom(embedding_energy_atom, pair_energy_atom, all_potential_atom);

  //for (Id i = 0; i < all_potential_atom.GetNumberOfValues(); i++)
  //{
  //  std::cout << "id" << i << "，all_potential_atom= " << all_potential_atom.ReadPortal().Get(i) << std::endl;
  //}
}

void ThermoOutput::StatisticalStatus()
{
  if (_executioner->CurrentStep() >= 1)
  {
    _status[POTENTIAL_ENERGY].push_back(_potential_energy);
    _status[TOTAL_ENERGY].push_back(_total_energy);
    _status[KIN_ENERGY].push_back(_kinteic_energy);
    _status[TEMPERATURE].push_back(_temperature);
  }
}

void ThermoOutput::PostData() 
{
  auto current_step = _executioner->CurrentStep();
  if (current_step == 1)
  {
    _system_state << "   "
               << " , "
               << "Initial"
               << " , "
               << "Final"
               << " , "
               << "Average"
               << " , "
               << "Std.Dev." 
               << " , "
               << "state" 
               << " , "
               << "finish time" 
               << " , "
               << "GPU time" 
               << " , "
               << "CPU time" << std::endl;
  }
  if (current_step == _executioner->NumStep())
  {
    PostExecute();
  }
}

void ThermoOutput::PostExecute() 
{
  try
  {
    //status
    //total energy
    auto init_total_energy = *((_status[TOTAL_ENERGY]).begin());
    auto final_total_energy = *((_status[TOTAL_ENERGY]).end() - 1);
    auto average_total_energy = std::accumulate(_status[TOTAL_ENERGY].begin(), _status[TOTAL_ENERGY].end(), (Real)0.0) 
        / _status[TOTAL_ENERGY].size();

    Real sum_std_total_energy = 0.0;
    std::for_each(_status[TOTAL_ENERGY].begin(), _status[TOTAL_ENERGY].end(), [&](const auto& energy) -> void{ 
        sum_std_total_energy += std::pow(energy - average_total_energy, 2); });
    auto std_total_energy = std::sqrt(sum_std_total_energy / (_status[TOTAL_ENERGY]).size());

    //potential energy
    auto init_potential_energy = *((_status[POTENTIAL_ENERGY]).begin());
    auto final_potential_energy = *((_status[POTENTIAL_ENERGY]).end() - 1);
    auto average_potential_energy = std::accumulate(_status[POTENTIAL_ENERGY].begin(),_status[POTENTIAL_ENERGY].end(),(Real)0.0) 
        / _status[POTENTIAL_ENERGY].size();
    _para.SetParameter(gtest::ave_potential_energy, average_potential_energy);
    Real sum_std_potential_energy = 0.0;
    std::for_each(_status[POTENTIAL_ENERGY].begin(), _status[POTENTIAL_ENERGY].end(), [&](const auto& energy) -> void
                  { sum_std_potential_energy += std::pow(energy - average_potential_energy, 2); });
    auto std_potential_energy = std::sqrt(sum_std_potential_energy / (_status[POTENTIAL_ENERGY]).size());

    //kin energy
    auto init_kin_energy = *((_status[KIN_ENERGY]).begin());
    auto final_kin_energy = *((_status[KIN_ENERGY]).end() - 1);
    auto average_kin_energy = std::accumulate(_status[KIN_ENERGY].begin(),_status[KIN_ENERGY].end(),(Real)0.0) 
        / _status[KIN_ENERGY].size();
    _para.SetParameter(gtest::ave_kin_energy, average_kin_energy);
    Real sum_std_kin_energy = 0.0;
    std::for_each(_status[KIN_ENERGY].begin(), _status[KIN_ENERGY].end(), [&](const auto& energy) -> void
                  { sum_std_kin_energy += std::pow(energy - average_kin_energy, 2); });
    auto std_kin_energy = std::sqrt(sum_std_kin_energy / (_status[KIN_ENERGY]).size());

    //temperature
    auto init_temperature = *((_status[TEMPERATURE]).begin());
    auto final_temperature = *((_status[TEMPERATURE]).end() - 1);
    auto average_temperature = std::accumulate(_status[TEMPERATURE].begin(), _status[TEMPERATURE].end(), (Real)0.0) 
        / _status[TEMPERATURE].size();
    _para.SetParameter(gtest::ave_temp_t ,average_temperature);
    Real sum_std_temperature = 0.0;
    std::for_each(_status[TEMPERATURE].begin(),_status[TEMPERATURE].end(), [&](const auto& energy) -> void
                  { sum_std_temperature += std::pow(energy - average_temperature, 2); });
    auto std_temperature = std::sqrt(sum_std_temperature / (_status[TEMPERATURE]).size());

    std::time_t currentTime = std::time(nullptr);
    std::tm* localTime = std::localtime(&currentTime);
    char buffer[80];
    std::strftime(buffer, sizeof(buffer), "%Y/%m/%d   %H:%M:%S", localTime);

      _system_state << "Tot.energy (kcal/mol)"   << " , " << init_total_energy << " , " << final_total_energy << " , " << average_total_energy << " , " << std_total_energy <<" , "
                    <<"finish"<< " , " <<buffer<< " , "<<0<< " , "<< _cumulative_time / 60 << "minites  "<<"("<<_cumulative_time<<"s"<<")" << std::endl
                    << "Pot.energy (kcal/mol)"   << " , " << init_potential_energy << " , " << final_potential_energy << " , " << average_potential_energy << " , " << std_potential_energy << std::endl
                    << "Kin.energy (kcal/mol)"   << " , " << init_kin_energy << " , " << final_kin_energy << " , " << average_kin_energy << " , " << std_kin_energy << std::endl
                    << "Temperature (K)"         << " , " << init_temperature << " , " <<final_temperature << " , " <<average_temperature << " , " <<std_temperature << std::endl
                    << "Volume (ANG^3)"          << " , " << _volume<< " , " << _volume << " , " << _volume << " , " << 0 <<std::endl
                    << "Density (1/ANG)"         << " , " << _rho << " , " << _rho << " , " << _rho << " , " << 0<<std::endl;

  }
  catch (const std::exception& e)
    {
      _system_state.close();
      console::Error(e.what());
    }
}
