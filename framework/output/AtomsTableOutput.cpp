#include "Executioner.h"
#include "AtomsTableOutput.h"
#include <vtkm/cont/Algorithm.h>
#include "Application.h"
#include "ConsoleOutput.h"
#include "math/Math.h"
#include "output/worklet/OutPutWorklet.h"
#include "FieldName.h"
#include "UnitFactor.h"

const std::string HEADER_POTENTIAL_ENERGY_NAME = "POTENTIAL_ENERGY";
const std::string HEADER_RESIDUAL_NAME = "RESIDUAL";
const std::string HEADER_KINETIC_ENERGY_NAME = "KINETIC_ENERGY";
const std::string HEADER_TOTAL_ENERGY_NAME = "TOTAL_ENERGY";
const std::string HEADER_KBT_NAME = "temperature";
const std::string HEADER_PRESSURE_NAME = "pressure";
const std::string HEADER_CUMULATIVE_TIME_NAME = "CUMULATIVE_TIME";

//RegisterObject(AtomsTableOutput);

AtomsTableOutput::AtomsTableOutput(const Configuration& cfg)
  : ConsoleOutput(cfg)
  , _binary(Get<bool>("binary", false))
  , _file(_name + ".csv")
  , _interval(Get<int>("interval", 1))
  , _out_initial(Get<bool>("out_initial", false))
  , _output_file(Get<bool>("output_file"))
  , _compute(Get<bool>("compute"))
  , _cut_off(0.0)
  , _Vlength(0.0)
  , _volume(0.0)
  , _Kmax(0)
  , _alpha(0.0)
  , _rho(0.0)
  , _tempT_sum(0.0)
  , _tempT(0.0)
  , _pressure(0.0)
{

}

AtomsTableOutput::~AtomsTableOutput()
{
  _file.close();
}

void AtomsTableOutput::Init()
{
  if ( _compute)  //_output_screen &&
  {
    ConsoleOutput::Init();
    AddHeader(HEADER_POTENTIAL_ENERGY_NAME);
    //AddHeader(HEADER_RESIDUAL_NAME);
    AddHeader(HEADER_KBT_NAME);
    AddHeader(HEADER_PRESSURE_NAME);
    AddHeader(HEADER_KINETIC_ENERGY_NAME);
    AddHeader(HEADER_TOTAL_ENERGY_NAME);
    //AddHeader(HEADER_CUMULATIVE_TIME_NAME);

    _cut_off = _system.GetParameter<Real>(PARA_CUTOFF);
    _volume = _system.GetParameter<Real>(PARA_VOLUME);
    _Vlength = _system.GetParameter<Real>(PARA_VLENGTH);
    _Kmax = _system.GetParameter<IdComponent>(PARA_KMAX);
    _alpha = _system.GetParameter<Real>(PARA_ALPHA);
    _rho = _system.GetParameter<Real>(PARA_RHO);
  }
}

void AtomsTableOutput::Execute()
{
  if (_compute)
  {
    _tempT_sum = _system.GetParameter<Real>(PARA_TEMPT_SUM);
    _tempT = _system.GetParameter<Real>(PARA_TEMPT);
    _pressure = _system.GetParameter<Real>(PARA_PRESSURE);

    ComputePotentialEnergy();          ////////////////////         这里需要有判断
    //ComputeEAMPotentialEnergy();

    Residual();

    AddDataToTable();

    WriteToFile();
  }
  ConsoleOutput::Execute();
}

bool AtomsTableOutput::ShouldOutput()
{
  auto current_step = _executioner->CurrentStep();
  if (current_step == 0)
    return _out_initial;
  else if (current_step == _executioner->NumStep())
    return true;
  else
    return current_step % _interval == 0;
}

void AtomsTableOutput::ComputePotentialEnergy()
{ 
  auto atoms_id = _system.GetFieldAsArrayHandle<Id>(field::atom_id);
  auto charge = _system.GetFieldAsArrayHandle<Real>(field::charge);
  auto position = _system.GetFieldAsArrayHandle<Vec3f>(field::position);
  ArrayHandle<Real> lj_potential_energy;

  ContForceFunction force_function;

  ContTopology topology;
  SetTopology(topology);

  ContPointLocator locator;
  SetLocator(locator);
  OutPut::ComputePotentialEnergy(_cut_off, atoms_id,  locator, topology, force_function, lj_potential_energy);  

   //auto  range = _system.GetParameter<vtkm::Vec<vtkm::Range, 3>>(PARA_RANGE);
   //auto Vlength = range[0].Max - range[0].Min;
   //OutPut::ComputePotentialEnergyPBC(
    // _cut_off, Vlength, atoms_id, locator, topology, force_function, lj_potential_energy);  
  
  auto _lj_potential_energy_total = vtkm::cont::Algorithm::Reduce( lj_potential_energy, vtkm::TypeTraits<Real>::ZeroInitialization()); 
  _lj_potential_energy_avr = _lj_potential_energy_total / position.GetNumberOfValues();

  //TODO: turn to parameter
  auto N = position.GetNumberOfValues();

  auto unit_factor = _system.GetParameter<UnitFactor>(PARA_UNIT_FACTOR);
  // self potential energy
  ArrayHandle<Real> _self_energy;
  OutPut::ComputeSqCharge(charge, _self_energy);
  auto _self_potential_energy_total = -vtkm::Sqrt(_alpha / vtkm::Pi()) *
    vtkm::cont::Algorithm::Reduce(_self_energy, vtkm::TypeTraits<Real>::ZeroInitialization());
  _self_potential_energy_avr = _self_potential_energy_total / N;
  _self_potential_energy_avr = _self_potential_energy_avr * unit_factor._qqr2e;

  //_PotentialTimer.Start();
  ArrayHandle<Real> _near_ele_potential_energy;
  OutPut::ComputeNearElePotential(_cut_off, _alpha, atoms_id, locator, topology ,force_function, _near_ele_potential_energy);

  auto _near_ele_potential_energy_total = vtkm::cont::Algorithm::Reduce(
    _near_ele_potential_energy, vtkm::TypeTraits<Real>::ZeroInitialization());
  _near_ele_potential_energy_avr =
    _near_ele_potential_energy_total / position.GetNumberOfValues();
  _near_ele_potential_energy_avr = _near_ele_potential_energy_avr * unit_factor._qqr2e;

  _spec_far_ele_potential_energy_avr = 0.0;

  _near_ele_potential_energy_avr = _near_ele_potential_energy_avr - _spec_far_ele_potential_energy_avr;

  Real Volume = vtkm::Pow(_Vlength, 3);
  ArrayHandle<Real> Density_Real;
  ArrayHandle<Real> Density_Image;
  Real _far_ele_potential_energy_total = 0.0;
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
          _far_ele_potential_energy_total +=
            vtkm::Exp(-Range_K * Range_K / (4 * _alpha)) * Range_density2 /
            (Range_K * Range_K);
        }
      }
    }
  }

  _far_ele_potential_energy_total = _far_ele_potential_energy_total * (2 * vtkm::Pi() / Volume);
  _far_ele_potential_energy_avr = _far_ele_potential_energy_total / N;
  _far_ele_potential_energy_avr = _far_ele_potential_energy_avr * unit_factor._qqr2e;

  _far_ele_potential_energy_avr = _far_ele_potential_energy_avr + _self_potential_energy_avr;

  _potential_energy = _lj_potential_energy_avr + _far_ele_potential_energy_avr +
    _near_ele_potential_energy_avr;

  _kinteic_energy = _tempT_sum / N * 0.5;

  _total_energy = _potential_energy + _kinteic_energy;

  _temperature = _tempT;


}

void AtomsTableOutput::Residual()
{
  if (_executioner->CurrentStep() <= 0)
  {
    _residual = 0;
  }
  if (_executioner->CurrentStep() > 0)
  {
    _residual = _potential_energy - _potential_energy_avr_old;
  }
  _potential_energy_avr_old = _potential_energy;
}

void AtomsTableOutput::AddDataToTable()
{
  if (_output_screen)
  {
    AddData(HEADER_POTENTIAL_ENERGY_NAME, _potential_energy);
    //AddData(HEADER_RESIDUAL_NAME, _residual);
    AddData(HEADER_TOTAL_ENERGY_NAME, _total_energy);
    AddData(HEADER_KINETIC_ENERGY_NAME, _kinteic_energy);
    AddData(HEADER_KBT_NAME, _tempT);
    AddData(HEADER_PRESSURE_NAME, _pressure);
    //AddData(HEADER_CUMULATIVE_TIME_NAME, _cumulative_time + _timer.GetElapsedTime());
  }
}

void AtomsTableOutput::WriteToFile()
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
            << "temperature "                         // add kBT
            << ", "
            << "pressure"
            << ", "
            << "KinteicEnergy"
            << ", "
            << "PotentialEnergy"
            << ", "
            << "TotalEnergy"
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
            << _pressure << ", "
            << _kinteic_energy << ", "
            << _potential_energy << ", "
            << _total_energy << ", "
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

void AtomsTableOutput::ComputeEAMPotentialEnergy()
{
  auto atoms_id = _system.GetFieldAsArrayHandle<Id>(field::atom_id);
  auto N = atoms_id.GetNumberOfValues();
  //
  auto cut_off = _system.GetParameter<Real>(EAM_PARA_CUTOFF);
  auto Vlength = _system.GetParameter<Real>(PARA_VLENGTH);

  auto rhor_spline = _system.GetFieldAsArrayHandle<Vec7f>(field::rhor_spline);
  auto frho_spline = _system.GetFieldAsArrayHandle<Vec7f>(field::frho_spline);
  auto z2r_spline = _system.GetFieldAsArrayHandle<Vec7f>(field::z2r_spline);

  //
  ContForceFunction force_function;
  SetForceFunction(force_function);

  ContTopology topology;
  SetTopology(topology);

  ContPointLocator locator;
  SetLocator(locator);

  //1:compute EAM_rho;
  ArrayHandle<Real> EAM_rho;
  OutPut::EAM_rho(
    cut_off, Vlength, atoms_id, rhor_spline, locator, topology, force_function, EAM_rho);

  //2: embedding_energy_atom
  ArrayHandle<Real> embedding_energy_atom;
  OutPut::EAM_EmbeddingEnergy(
    atoms_id, EAM_rho, frho_spline, locator, topology, force_function, embedding_energy_atom);

  auto embedding_energy_total = vtkm::cont::Algorithm::Reduce(
    embedding_energy_atom, vtkm::TypeTraits<Real>::ZeroInitialization());

  //3: pair_energy_atom
  ArrayHandle<Real> pair_energy_atom;
  OutPut::EAM_PairEnergy(
    cut_off, Vlength, atoms_id, z2r_spline, locator, topology, force_function, pair_energy_atom);

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
