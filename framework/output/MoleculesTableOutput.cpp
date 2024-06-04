#include "Executioner.h"
#include "MoleculesTableOutput.h"
#include <vtkm/cont/Algorithm.h>
#include "Application.h"
#include "ConsoleOutput.h"
#include "math/Math.h"
#include "output/worklet/OutPutWorklet.h"
#include "FieldName.h"
#include <ctime>
#include "forceFunction/ContForceFunction.h"
#include "UnitFactor.h"

const std::string HEADER_POTENTIAL_ENERGY_NAME = "POTENTIAL_ENERGY";
const std::string HEADER_RESIDUAL_NAME = "RESIDUAL";
const std::string HEADER_KINETIC_ENERGY_NAME = "KINETIC_ENERGY";
const std::string HEADER_NON_BOND_ENERGY_NAME = "NON_BOND_ENERGY";
const std::string HEADER_TOTAL_ENERGY_NAME = "TOTAL_ENERGY";
const std::string HEADER_BOND_ENERGY_NAME = "BOND_ENERGY";
const std::string HEADER_ANGLE_ENERGY_NAME = "ANGLE_ENERGY";
const std::string HEADER_KBT_NAME = "KBT";
const std::string HEADER_CUMULATIVE_TIME_NAME = "CUMULATIVE_TIME";

//RegisterObject(MoleculesTableOutput);




MoleculesTableOutput::MoleculesTableOutput(const Configuration& cfg)
  : ConsoleOutput(cfg)
  , _binary(Get<bool>("binary", false))
  , _file(_name + ".csv")
  , _interval(Get<int>("interval", 1))
  , _out_initial(Get<bool>("out_initial", false))
  , _output_file(Get<bool>("output_file"))
  , _compute(Get<bool>("compute"))
  , _system_state("SystemState.csv")
  , _cut_off(0.0)
  , _Vlength(0.0)
  , _volume(0.0)
  , _Kmax(0)
  , _alpha(0.0)
  , _bond_energy(0.0)
  , _angle_energy(0.0)
  , _rho(0.0)
  , _tempT_sum(0.0)
  , _tempT(0.0)
{

}

MoleculesTableOutput::~MoleculesTableOutput()
{
  _file.close();
  _system_state.close();
}

void MoleculesTableOutput::Init()
{
  if ( _compute)  //_output_screen &&
  {
    ConsoleOutput::Init();
    AddHeader(HEADER_POTENTIAL_ENERGY_NAME);
    AddHeader(HEADER_RESIDUAL_NAME);
    AddHeader(HEADER_KBT_NAME);
    AddHeader(HEADER_KINETIC_ENERGY_NAME);
    AddHeader(HEADER_NON_BOND_ENERGY_NAME);
    AddHeader(HEADER_TOTAL_ENERGY_NAME);
    AddHeader(HEADER_BOND_ENERGY_NAME);
    AddHeader(HEADER_ANGLE_ENERGY_NAME);
    AddHeader(HEADER_CUMULATIVE_TIME_NAME);

    _cut_off = _para.GetParameter<Real>(PARA_CUTOFF);
    _volume = _para.GetParameter<Real>(PARA_VOLUME);
    _Vlength = _para.GetParameter<Real>(PARA_VLENGTH);
    _Kmax = _para.GetParameter<IdComponent>(PARA_KMAX);
    _alpha = _para.GetParameter<Real>(PARA_ALPHA);
    _rho = _para.GetParameter<Real>(PARA_RHO);
  }
}

void MoleculesTableOutput::Execute()
{
  if (_compute)
  {
    if (_para.HaveParameter(PARA_BOND_ENERGY))
    {
      _bond_energy = _para.GetParameter<Real>(PARA_BOND_ENERGY);
      
    }
    if (_para.HaveParameter(PARA_ANGLE_ENERGY))
    {
      _angle_energy = _para.GetParameter<Real>(PARA_ANGLE_ENERGY);
    }

    _tempT_sum = _para.GetParameter<Real>(PARA_TEMPT_SUM);
    _tempT = _para.GetParameter<Real>(PARA_TEMPT);

    ComputePotentialEnergy();

    Residual();

    AddDataToTable();

    WriteToFile();

    StatisticalStatus();
  }
  ConsoleOutput::Execute();
  PostData();
}

bool MoleculesTableOutput::ShouldOutput()
{
  auto current_step = _executioner->CurrentStep();
  if (current_step == 0)
    return _out_initial;
  else if (current_step == _executioner->NumStep())
    return true;
  else
    return current_step % _interval == 0;
}

void MoleculesTableOutput::ComputePotentialEnergy()
{ 
  ArrayHandle<Real> lj_potential_energy;
  auto atoms_id = _para.GetFieldAsArrayHandle<Id>(field::atom_id);
  auto position = _para.GetFieldAsArrayHandle<Vec3f>(field::position);

  ContForceFunction force_function;

  ContTopology topology;
  SetTopology(topology);

  ContPointLocator locator;
  SetLocator(locator);

  OutPut::ComputePotentialEnergy( _cut_off, atoms_id, locator, topology, force_function, lj_potential_energy);

  auto lj_potential_energy_total =
    vtkm::cont::Algorithm::Reduce(lj_potential_energy, vtkm::TypeTraits<Real>::ZeroInitialization());
  _lj_potential_energy_avr = lj_potential_energy_total / position.GetNumberOfValues();

  //TODO: turn to parameter
  auto N = position.GetNumberOfValues();

  auto unit_factor = _para.GetParameter<UnitFactor>(PARA_UNIT_FACTOR);
  // self potential energy
  ArrayHandle<Real> _self_energy;
  auto charge = _para.GetFieldAsArrayHandle<Real>(field::charge);
  OutPut::ComputeSqCharge(charge, _self_energy);
  auto self_potential_energy_total = -vtkm::Sqrt(_alpha / vtkm::Pi()) *
    vtkm::cont::Algorithm::Reduce(_self_energy, vtkm::TypeTraits<Real>::ZeroInitialization());
  _self_potential_energy_avr = self_potential_energy_total / N;
  _self_potential_energy_avr = _self_potential_energy_avr * unit_factor._qqr2e;

  //_PotentialTimer.Start();
  ArrayHandle<Real> near_ele_potential_energy;
  OutPut::ComputeNearElePotential(
    _cut_off, _alpha, atoms_id, locator, topology, force_function, near_ele_potential_energy);

  auto near_ele_potential_energy_total = vtkm::cont::Algorithm::Reduce(
    near_ele_potential_energy, vtkm::TypeTraits<Real>::ZeroInitialization());
  _near_ele_potential_energy_avr =
    near_ele_potential_energy_total / position.GetNumberOfValues();
  _near_ele_potential_energy_avr = _near_ele_potential_energy_avr * unit_factor._qqr2e;

  _spec_far_ele_potential_energy_avr = 0.0;

  SpecialFarCoulEnergy();

  _near_ele_potential_energy_avr = _near_ele_potential_energy_avr - _spec_far_ele_potential_energy_avr;

  Real Volume = vtkm::Pow(_Vlength, 3);
  ArrayHandle<Real> Density_Real;
  ArrayHandle<Real> Density_Image;
  Real far_ele_potential_energy_total = 0.0;
  for (Id i = -_Kmax; i <= _Kmax; i++)
  {
    for (Id j = -_Kmax; j <= _Kmax; j++)
    {
      for (Id k = -_Kmax; k <= _Kmax; k++)
      {
        if (!(i == 0 && j == 0 && k == 0))
        {
          Vec3f K = { Real(i), Real(j), Real(k) };
          K = 2 * vtkm::Pi() * K / _Vlength;
          Real Range_K = vtkm::Magnitude(K);
          OutPut::ComputeDensity(
            K, position, charge, Density_Real, Density_Image);
          Real Value_Re = vtkm::cont::Algorithm::Reduce(
            Density_Real, vtkm::TypeTraits<Real>::ZeroInitialization());
          Real Value_Im = vtkm::cont::Algorithm::Reduce(
            Density_Image, vtkm::TypeTraits<Real>::ZeroInitialization());
          Real Range_density2 = vtkm::Pow(Value_Re, 2) + vtkm::Pow(Value_Im, 2);

          far_ele_potential_energy_total +=
            vtkm::Exp(-Range_K * Range_K / (4 * _alpha)) * Range_density2 /
            (Range_K * Range_K);
        }
      }
    }
  }

  far_ele_potential_energy_total = far_ele_potential_energy_total * (2 * vtkm::Pi() / Volume);
  _far_ele_potential_energy_avr = far_ele_potential_energy_total / N;
  _far_ele_potential_energy_avr = _far_ele_potential_energy_avr * unit_factor._qqr2e;

  _far_ele_potential_energy_avr = _far_ele_potential_energy_avr + _self_potential_energy_avr;

  _potential_energy_avr = _lj_potential_energy_avr + _far_ele_potential_energy_avr +
    _near_ele_potential_energy_avr;

  _potential_energy = _potential_energy_avr + _bond_energy + _angle_energy;

  _kinteic_energy = _tempT_sum / N * 0.5;

  _total_energy = _potential_energy + _kinteic_energy;

  _temperature = _tempT;
}

void MoleculesTableOutput::SpecialFarCoulEnergy()
{
  auto position = _para.GetFieldAsArrayHandle<Vec3f>(field::position);
  auto N = position.GetNumberOfValues();

  Real Volume = vtkm::Pow(_Vlength, 3);

  ArrayHandle<Real> Spec_far_coul_energy;
  Real _spec_far_ele_potential_energy_total = 0.0;
  auto atoms_id = _para.GetFieldAsArrayHandle<Id>(field::atom_id);
  auto unit_factor = _para.GetParameter<UnitFactor>(PARA_UNIT_FACTOR);

  auto source_array = _para.GetFieldAsArrayHandle<Id>(field::special_source_array);
  auto offsets_array = _para.GetFieldAsArrayHandle<Id>(field::special_offsets_array);
  auto groupVecArray = vtkm::cont::make_ArrayHandleGroupVecVariable(source_array, offsets_array);

  ContForceFunction force_function;

  ContTopology topology;
  SetTopology(topology);

  ContPointLocator locator;
  SetLocator(locator);

  OutPut::ComputeSpecialFarCoul(_Vlength, atoms_id, groupVecArray, locator, topology, force_function, Spec_far_coul_energy); 

  _spec_far_ele_potential_energy_total = vtkm::cont::Algorithm::Reduce(
    Spec_far_coul_energy, vtkm::TypeTraits<Real>::ZeroInitialization());
  _spec_far_ele_potential_energy_avr =
    _spec_far_ele_potential_energy_total / position.GetNumberOfValues();
  _spec_far_ele_potential_energy_avr =
    0.5 * _spec_far_ele_potential_energy_avr * unit_factor._qqr2e;
}

void MoleculesTableOutput::Residual()
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

void MoleculesTableOutput::AddDataToTable()
{
  if (_output_screen)
  {
    AddData(HEADER_POTENTIAL_ENERGY_NAME, _potential_energy);
    AddData(HEADER_RESIDUAL_NAME, _residual);
    AddData(HEADER_NON_BOND_ENERGY_NAME, _potential_energy_avr);
    AddData(HEADER_TOTAL_ENERGY_NAME, _total_energy);
    AddData(HEADER_KINETIC_ENERGY_NAME, _kinteic_energy);
    AddData(HEADER_KBT_NAME, _tempT);
    AddData(HEADER_BOND_ENERGY_NAME, _bond_energy);
    AddData(HEADER_ANGLE_ENERGY_NAME, _angle_energy);
    AddData(HEADER_CUMULATIVE_TIME_NAME, _cumulative_time + _timer.GetElapsedTime());
  }
}

void MoleculesTableOutput::WriteToFile()
{
  if (ShouldOutput()&&_output_file)
  {
    if (_executioner->CurrentStep() == 1)
    {
      _file << "Step"
            << " , "
            << "Time"
            << " , "
            << "potentialLJEnergyAvr"        // lj energy
            << " , "
            << "potentialNearEleEnergyAvr"   // add near energy
            << " , "                         //
            << "PotentialFarEleEnergyAvr"    // add far energy
            << " , "                         //
            << "Residual" 
            << ", "
            << "kBT"                         // add kBT
            << ", "
            << "KinteicEnergy"
            << ", "
            << "PotentialEnergy"
            << ", "
            << "NonBonEnergy"
            << ", "
            << "TotalEnergy"
            << ", "
            << "BondEnergy"
            << ", "
            << "AngleEnergy"
            << ", "
            << "CumulativeTime"
            << ", "
            << "SpecialEnergy"
            << std::endl;
    }
    try
    {
      _file << _executioner->CurrentStep() << " , " 
            << _time<< " , "
            << _lj_potential_energy_avr << " , "
            << _near_ele_potential_energy_avr << " , "
            << _far_ele_potential_energy_avr << " , " 
            << _residual << ", "
            << _tempT  << ", "
            << _kinteic_energy << ", "
            << _potential_energy << ", "
            << _potential_energy_avr << ", "
            << _total_energy << ", "
            << _bond_energy << ", "
            << _angle_energy << ", "
            << _cumulative_time + _timer.GetElapsedTime() << ", "
            << _spec_far_ele_potential_energy_avr 
            << std::endl;
    }
    catch (const std::exception& e)
    {
      _file.close();
      console::Error(e.what());
    }
  }
}

void MoleculesTableOutput::StatisticalStatus()
{
  if (_executioner->CurrentStep() >= 1)
  {
    _status[POTENTIAL_ENERGY].push_back(_potential_energy);
    _status[TOTAL_ENERGY].push_back(_total_energy);
    _status[KIN_ENERGY].push_back(_kinteic_energy);
    _status[TEMPERATURE].push_back(_temperature);
  }
}

void MoleculesTableOutput::PostData() 
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

void MoleculesTableOutput::PostExecute() 
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
