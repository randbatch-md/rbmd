#include "EAMSystem.h"
#include "Application.h"
#include <fstream>
#include <vtkm/Pair.h>
#include <vtkm/VectorAnalysis.h>
#include <vtkm/cont/Algorithm.h>
#include <vtkm/cont/ArrayHandleGroupVecVariable.h>
#include <vtkm/cont/ConvertNumComponentsToOffsets.h>
#include <vtkm/Math.h>
#include <cmath>  // erfc(x)
#include "math/Math.h"
#include "math/Utiles.h"
#include "RBEPSample.h"
#include "system/worklet/SystemWorklet.h"
#include "MeshFreeCondition.h"
#include "FieldName.h"
#include "system/worklet/MolecularWorklet.h"


//RegisterObject(EAMSystem);
#define MIN(A, B) ((A) < (B) ? (A) : (B))
#define MAX(A, B) ((A) > (B) ? (A) : (B))

EAMSystem::EAMSystem(const Configuration& cfg)
  : MDSystem(cfg) 
  , _RBE_P(Get<IdComponent>("rbeP"))  
  , _executioner((_app.GetExecutioner()))
  , _kbT(Get<IdComponent>("kbT"))
  , _Kmax(Get<IdComponent>("kmax"))
  , _alpha(Get<Real>("alpha"))
  , _farforce_type(Get<std::string>("farforce_type"))
  , _nearforce_type(Get<std::string>("nearforce_type"))
  , _temp_con_type(Get<std::string>("temp_con_type"))
  , _potential_file(Get<std::string>("potential_file"))
{
  SetParameter(PARA_RBE_P, _RBE_P);
  SetParameter(PARA_ALPHA, _alpha);
  SetParameter(PARA_KMAX, _Kmax);
  SetParameter(PARA_TEMPT_SUM, Real{ 0.0 });
  SetParameter(PARA_TEMPT, Real{ 0.0 });
}

void EAMSystem::Init()
{
  //ReadPotentialFile:
  ReadPotentialFile(_potential_file);
  InitStyle();

  MDSystem::Init();

  //init variable
  InitialCondition();

  ComputeForce(); // Presolve force

}

void EAMSystem::InitialCondition()
{
  MDSystem::InitialCondition();
  _rho = GetParameter<Real>(PARA_RHO);
  _nosehooverxi = 0.0;

  SetParameter(PARA_RDF_RHO, _rho);
  SetCharge();
  PreForce();
}

void EAMSystem::ReadPotentialFile(std::ifstream& input_file)
{
  // 检查文件是否打开成功
  if (!input_file.is_open())
  {
    std::cerr << "Unable to open the file." << std::endl;
  }

  // 跳过前两行
  for (int i = 0; i < 2; ++i)
  {
    input_file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  }

  // 开始读取第三行的值
  input_file >> file.nrho >> file.drho >> file.nr >> file.dr >> file.cut_off;

  //
  file.frho.resize(file.nrho + 1);
  file.zr.resize(file.nr + 1);
  file.rhor.resize(file.nrho + 1);

  // 读取并保存 frho 数组

  for (int i = 0; i < file.nrho; ++i)
  {
    input_file >> file.frho[i];
  }

  // 读取并保存 zr 数组

  for (int i = 0; i < file.nr; ++i)
  {
    input_file >> file.zr[i];
  }

  // 读取并保存 rhor 数组

  for (int i = 0; i < file.nrho; ++i)
  {
    input_file >> file.rhor[i];
  }

  // 关闭文件
  input_file.close();
}

void EAMSystem::AllocateEAM() {}

void EAMSystem::file2array()
{
  Id i, j, k, m, n;
  Real sixth = 1.0 / 6.0;
  // auto ntypes = _header._num_atoms_type;

  Real rmax;
  dr = drho = rmax = rhomax = 0.0;

  dr = MAX(dr, file.dr);
  drho = MAX(drho, file.drho);
  rmax = MAX(rmax, (file.nr - 1) * file.dr);
  rhomax = MAX(rhomax, (file.nrho - 1) * file.drho);

  // set nr,nrho from cutoff and spacings
  // 0.5 is for round-off in divide

  nr = static_cast<int>(rmax / dr + 0.5);
  nrho = static_cast<int>(rhomax / drho + 0.5);

  // ------------------------------------------------------------------
  // setup frho arrays
  // ------------------------------------------------------------------
  frho.resize(nrho + 1);

  Real r, p, cof1, cof2, cof3, cof4;
  for (m = 1; m <= nrho; m++)
  {
    r = (m - 1) * drho;
    p = r / file.drho + 1.0;
    k = static_cast<int>(p);
    k = MIN(k, file.nrho - 2);
    k = MAX(k, 2);
    p -= k;
    p = MIN(p, 2.0);
    cof1 = -sixth * p * (p - 1.0) * (p - 2.0);
    cof2 = 0.5 * (p * p - 1.0) * (p - 2.0);
    cof3 = -0.5 * p * (p + 1.0) * (p - 2.0);
    cof4 = sixth * p * (p * p - 1.0);
    frho[m] = cof1 * file.frho[k - 1] + cof2 * file.frho[k] + cof3 * file.frho[k + 1] +
      cof4 * file.frho[k + 2];
  }

  // ------------------------------------------------------------------
  // setup rhor arrays
  // ------------------------------------------------------------------
  rhor.resize(nrho + 1);
  for (m = 1; m <= nr; m++)
  {
    r = (m - 1) * dr;
    p = r / file.dr + 1.0;
    k = static_cast<int>(p);
    k = MIN(k, file.nr - 2);
    k = MAX(k, 2);
    p -= k;
    p = MIN(p, 2.0);
    auto cof1 = -sixth * p * (p - 1.0) * (p - 2.0);
    auto cof2 = 0.5 * (p * p - 1.0) * (p - 2.0);
    auto cof3 = -0.5 * p * (p + 1.0) * (p - 2.0);
    auto cof4 = sixth * p * (p * p - 1.0);
    rhor[m] = cof1 * file.rhor[k - 1] + cof2 * file.rhor[k] + cof3 * file.rhor[k + 1] +
      cof4 * file.rhor[k + 2];
  }

  // ------------------------------------------------------------------
  // setup z2r arrays
  // ------------------------------------------------------------------
  z2r.resize(nr + 1);

  double zri;
  for (m = 1; m <= nr; m++)
  {
    r = (m - 1) * dr;

    p = r / file.dr + 1.0;
    k = static_cast<int>(p);
    k = MIN(k, file.nr - 2);
    k = MAX(k, 2);
    p -= k;
    p = MIN(p, 2.0);
    cof1 = -sixth * p * (p - 1.0) * (p - 2.0);
    cof2 = 0.5 * (p * p - 1.0) * (p - 2.0);
    cof3 = -0.5 * p * (p + 1.0) * (p - 2.0);
    cof4 = sixth * p * (p * p - 1.0);
    zri = cof1 * file.zr[k - 1] + cof2 * file.zr[k] + cof3 * file.zr[k + 1] + cof4 * file.zr[k + 2];

    z2r[m] = 27.2 * 0.529 * zri * zri;
  }
}

void EAMSystem::interpolate(Id n, Real delta, std::vector<Real>& f, std::vector<Vec7f>& spline)
{
  for (int m = 1; m <= n; m++)
  {
    spline[m][6] = f[m];
  }

  spline[1][5] = spline[2][6] -
    spline[1][6]; //f'(x) = (f(x + h) - f(x)) / h    [5] 为一阶导数的系数， 能量表达式的系数
  spline[2][5] = 0.5 * (spline[3][6] - spline[1][6]);
  spline[n - 1][5] = 0.5 * (spline[n][6] - spline[n - 2][6]);
  spline[n][5] = spline[n][6] - spline[n - 1][6];

  for (int m = 3; m <= n - 2; m++)
  {
    spline[m][5] =
      ((spline[m - 2][6] - spline[m + 2][6]) + 8.0 * (spline[m + 1][6] - spline[m - 1][6])) /
      12.0; //使用更远的样本点以获得更准确的估计
  }

  for (int m = 1; m <= n - 1; m++)
  {
    spline[m][4] = 3.0 * (spline[m + 1][6] - spline[m][6]) - 2.0 * spline[m][5] -
      spline[m + 1][5]; //[4] 为二阶导数的系数
    spline[m][3] = spline[m][5] + spline[m + 1][5] -
      2.0 * (spline[m + 1][6] - spline[m][6]); // [3]为三阶导数的系数
  }

  spline[n][4] = 0.0;
  spline[n][3] = 0.0; //最后一个样本点处的二阶和三阶导数的系数为零,
    //为了使插值曲线在两端更平滑，可以将边界处的高阶导数系数设置为零。
    //这是因为样条插值通常在内部样本点上使用高阶多项式插值，而在边界处使用较低阶的多项式以确保平滑性。

  for (int m = 1; m <= n; m++)
  {
    spline[m][2] = spline[m][5] / delta;       //二次导数的系数。   力表达式的系数
    spline[m][1] = 2.0 * spline[m][4] / delta; //一次导数的系数
    spline[m][0] = 3.0 * spline[m][3] / delta; //零次导数（即函数值）的系数。
  }
}

void EAMSystem::array2spline()
{
  frho_spline.resize(nrho + 1);
  rhor_spline.resize(nrho + 1);
  z2r_spline.resize(nr + 1);

  interpolate(nrho, drho, frho, frho_spline);
  interpolate(nr, dr, rhor, rhor_spline);
  interpolate(nr, dr, z2r, z2r_spline);
}

void EAMSystem::SetEAM()
{
  SetParameter(EAM_PARA_CUTOFF, file.cut_off);
  SetParameter(EAM_PARA_RHOMAX, rhomax);
  SetParameter(EAM_PARA_NRHO, nrho);
  SetParameter(EAM_PARA_DRHO, drho);
  SetParameter(EAM_PARA_NR, nr);
  SetParameter(EAM_PARA_DR, dr);
  //

  AddField(field::rhor_spline, ArrayHandle<Vec7f>{});
  auto rhor_spline_get = GetFieldAsArrayHandle<Vec7f>(field::rhor_spline);
  vtkm::cont::ArrayCopy(vtkm::cont::make_ArrayHandle<Vec7f>(rhor_spline), rhor_spline_get);

  AddField(field::frho_spline, ArrayHandle<Vec7f>{});
  auto frho_spline_get = GetFieldAsArrayHandle<Vec7f>(field::frho_spline);
  vtkm::cont::ArrayCopy(vtkm::cont::make_ArrayHandle<Vec7f>(frho_spline), frho_spline_get);

  AddField(field::z2r_spline, ArrayHandle<Vec7f>{});
  auto z2r_spline_get = GetFieldAsArrayHandle<Vec7f>(field::z2r_spline);
  vtkm::cont::ArrayCopy(vtkm::cont::make_ArrayHandle<Vec7f>(z2r_spline), z2r_spline_get);
}

void EAMSystem::InitStyle()
{
  AllocateEAM();
  file2array();
  array2spline();
  SetEAM();
}

void EAMSystem::ComputeForce()
{
  //ComputeAllForce();
  //TempConTypeForce();
  if (_nearforce_type == "RBL")
  {
    ComputeRBLEAMForce(_force);
  }
  else if (_nearforce_type == "VERLETLIST")
  {
    ComputeVerletlistEAMForce(_force);
  }
  else if (_nearforce_type == "ORIGINAL")
  {
    ComputeOriginalEAMForce(_force);
  }
}

void EAMSystem::UpdateVelocity()
{
  try
  {
    SystemWorklet::UpdateVelocity(_dt, _unit_factor._fmt2v ,_force, _mass, _velocity);
    ComputeTempe();
    UpdateVelocityByTempConType();
  }
  catch (const std::exception& e)
  {
    std::cout << e.what() << std::endl;
  }
}

void EAMSystem::UpdatePosition()
{
  auto&& position_flag = GetFieldAsArrayHandle<Id3>(field::position_flag);
  SystemWorklet::UpdatePositionFlag(_dt, _velocity, _locator, _position, position_flag);
  _locator.SetPosition(_position);
  SetCenterTargetPositions();
}

void EAMSystem::SetCharge()
{
  _charge = GetFieldAsArrayHandle<Real>(field::charge);
  auto n = _position.GetNumberOfValues();
  _charge.AllocateAndFill(n, 0);
  _topology.SetCharge(_charge);
}

void EAMSystem::ComputeTempe()
{
  auto n = _position.GetNumberOfValues();
  ArrayHandle<Real> sq_velocity;
  sq_velocity.Allocate(n);

  SystemWorklet::ComputerKineticEnergy(_velocity, _mass, sq_velocity);
  Invoker{}(MolecularWorklet::UnitRescaleWorklet{ _unit_factor._mvv2e }, sq_velocity);

  _tempT_sum =
    vtkm::cont::Algorithm::Reduce(sq_velocity, vtkm::TypeTraits<Real>::ZeroInitialization());
  
  Real temperature_kB = _unit_factor._kB;
  _tempT = 0.5 * _tempT_sum / (3 * n * temperature_kB / 2.0);

  SetParameter(PARA_TEMPT_SUM, _tempT_sum);
  SetParameter(PARA_TEMPT, _tempT);
}

void EAMSystem::UpdateVelocityByTempConType()
{
  if (_temp_con_type == "NOSE_HOOVER")
  {
    //Nose Hoover
    //Maybe tauT = 20.0 * dt is a good choice for dt = 2e-3 from the test.
    //Because the temperature curve is the smoothest of all test simulations.
    //In fact, 5.0 and 10.0 are also optional.
    //As long as the coefficent is not too large, such as larger than 100 * dt.
    SystemWorklet::UpdateVelocityNoseHoover(_dt, _unit_factor._fmt2v, _nosehooverxi, _force, _mass, _velocity);
    Real tauT = 5000 * _dt;           //2000 is is ok for Cu, Ni ; 5000 is ok for Ag , Au, Pd ; 10000  is ok for Pt
    _nosehooverxi += 0.5 * _dt * (_tempT / _kbT - 1.0) / tauT;
  }
  else if (_temp_con_type == "TEMP_RESCALE")
  {
    Real coeff_rescale = vtkm::Sqrt(_kbT / _tempT);
    SystemWorklet::UpdateVelocityRescale(coeff_rescale, _velocity);
  }
  else if (_temp_con_type == "BERENDSEN")
  {
    //
    //Velocity Rescale: Berendsen
    //Maybe dt_divide_taut = 0.05 is a good choice for dt = 2e-3. 0.005, 0.01, 0.1 is optional.
    //The selection of dt_divide_taut determines the temperature equilibrium time.
    //
    Real dt_divide_taut = 0.5;   // 0.01 is ok for Cu ; 0.5 is ok for Ni ,Ag 
    Real coeff_Berendsen = vtkm::Sqrt(1.0 + dt_divide_taut * (_kbT / _tempT - 1.0));
    SystemWorklet::UpdateVelocityRescale(coeff_Berendsen, _velocity);
  }
}

void EAMSystem::PreSolve()
{
  _locator.SetPosition(_position);
  PreForce();
}
 
void EAMSystem::Solve() 
{
  // stage1:
  UpdateVelocity();
  end_of_step();

  // stage2:
  UpdatePosition();

  // stage3:
  ComputeForce();
  UpdateVelocity();
}

void EAMSystem::PostSolve()
{}

void EAMSystem::PreForce()
{
  _Vlength = GetParameter<Real>(PARA_VLENGTH);
  _dt = _executioner->Dt();
  // prepare for RBE force
  auto velocity_type = GetParameter<std::string>(gtest::velocity_type);
  auto random = (velocity_type != "TEST") ? true : false;
  RBEPSAMPLE rbe_presolve_psample = { _alpha, _Vlength, _RBE_P };
  rbe_presolve_psample._RBE_random = random;
  _psample = rbe_presolve_psample.Fetch_P_Sample(
    Real(0.0), (vtkm::Sqrt(_alpha / 2.0) * _Vlength / vtkm::Pi()));

  // prepare for Langevin dynamics
  RBEPSAMPLE sample_presolve_1d;
  sample_presolve_1d._RBE_random = random;
  Real gaussianx = sample_presolve_1d.FetchSample_1D(Real(0.0), Real(1.0));
  Real gaussiany = sample_presolve_1d.FetchSample_1D(Real(0.0), Real(1.0));
  Real gaussianz = sample_presolve_1d.FetchSample_1D(Real(0.0), Real(1.0));
  _gaussian = { gaussianx, gaussiany, gaussianz };
}

void EAMSystem::InitField()
{
  MDSystem::InitField();
  AddField(field::pts_type , ArrayHandle<Id>{});
  AddField(field::position_flag, ArrayHandle<Id3>{});
  AddField(field::center_position, ArrayHandle<Vec3f>{});
  AddField(field::target_position, ArrayHandle<Vec3f>{});
  AddField(field::epsilon, ArrayHandle<Real>{});
  AddField(field::sigma, ArrayHandle<Real>{});
}

void EAMSystem::SetCenterTargetPositions()
{
  //auto num_pos = _position.GetNumberOfValues();
  //vtkm::cont::ArrayHandle<vtkm::Vec3f> center_position_temp;
  //center_position_temp.Allocate(num_pos / 2);
  //auto&& write_prot_center = center_position_temp.WritePortal();
  //auto&& read_prot_center = _position.ReadPortal();
  //
  //vtkm::cont::ArrayHandle<vtkm::Vec3f> target_position_temp;
  //target_position_temp.Allocate(num_pos / 2);
  //auto&& write_prot_target = target_position_temp.WritePortal();
  //auto&& read_prot_target = _position.ReadPortal();
  //
  //for (int i = 0; i < num_pos / 2; i++)
  //{
  //  write_prot_center.Set(i, read_prot_center.Get(i));
  //  write_prot_target.Set(i, read_prot_target.Get(i + (num_pos / 2)));
  //}
  //
  //auto center_position = GetFieldAsArrayHandle<Vec3f>(field::center_position);
  //vtkm::cont::ArrayCopy(center_position_temp, center_position);
  //
  //auto target_position = GetFieldAsArrayHandle<Vec3f>(field::target_position);
  //vtkm::cont::ArrayCopy(target_position_temp, target_position);

  auto center_position = GetFieldAsArrayHandle<Vec3f>(field::center_position);
  vtkm::cont::ArrayCopy(_position, center_position);

  auto target_position = GetFieldAsArrayHandle<Vec3f>(field::target_position);
  vtkm::cont::ArrayCopy(_position, target_position);
}

void EAMSystem::TimeIntegration() {}

void EAMSystem::SetForceFunction()
{
  InitERF();
  auto rhomax = GetParameter<Real>(EAM_PARA_RHOMAX);
  auto nrho = GetParameter<Id>(EAM_PARA_NRHO);
  auto drho = GetParameter<Real>(EAM_PARA_DRHO);
  auto nr = GetParameter<Id>(EAM_PARA_NR);
  auto dr = GetParameter<Real>(EAM_PARA_DR);

  _force_function.SetEAMParameters(rhomax, nrho, drho, nr, dr);
}

void EAMSystem::SetTopology()
{
  auto pts_type = GetFieldAsArrayHandle<Id>(field::pts_type);
  auto molecule_id = GetFieldAsArrayHandle<Id>(field::molecule_id);
  auto epsilon = GetFieldAsArrayHandle<Real>(field::epsilon);
  auto sigma = GetFieldAsArrayHandle<Real>(field::sigma);
  _topology.SetAtomsType(pts_type);
  _topology.SetMolecularId(molecule_id);
  _topology.SetEpsAndSigma(epsilon, sigma);
}

void EAMSystem::ComputeVirial()
{
  auto cut_off = GetParameter<Real>(EAM_PARA_CUTOFF);
  auto Vlength = GetParameter<Real>(PARA_VLENGTH);

  auto rhor_spline = GetFieldAsArrayHandle<Vec7f>(field::rhor_spline);
  auto frho_spline = GetFieldAsArrayHandle<Vec7f>(field::frho_spline);
  auto z2r_spline = GetFieldAsArrayHandle<Vec7f>(field::z2r_spline);

  ArrayHandle<Real> EAM_rho;
  ArrayHandle<Real> fp;

  //1:compute _EAM_rho   = density at each atom
  SystemWorklet::EAM_rho(
    cut_off, Vlength, _atoms_id, rhor_spline, _locator, _topology, _force_function, EAM_rho);

  // 2:compute fp    = derivative of embedding energy at each atom
  SystemWorklet::EAM_fp(_atoms_id, EAM_rho, frho_spline, _locator, _topology, _force_function, fp);

  // 3:compute  virial  = EAM_virial
  SystemWorklet::EAM_virial(cut_off,
                           Vlength,
                           _atoms_id,
                           fp,
                           rhor_spline,
                           z2r_spline,
                           _locator,
                           _topology,
                           _force_function,
                            _virial_atom);


  //for (int i = 0; i <_virial_atom.GetNumberOfValues();++i)
  //{
  //  std::cout << "i=" << i << ",_virial_atom=" << _virial_atom.ReadPortal().Get(i)[0] << "," 
  //            << _virial_atom.ReadPortal().Get(i)[1] << "," << _virial_atom.ReadPortal().Get(i)[2] << "," 
  //            << _virial_atom.ReadPortal().Get(i)[3] << "," << _virial_atom.ReadPortal().Get(i)[4] << "," 
  //            << _virial_atom.ReadPortal().Get(i)[5]
  //      << std::endl;
  //}
}

void EAMSystem::Compute_Pressure_Scalar()
{

  ComputeVirial();

  virial = { 0, 0, 0, 0, 0, 0 };
  for (int i = 0; i < _virial_atom.GetNumberOfValues(); ++i)
  {
    // 获取当前原子的virial
    Vec6f vatom = _virial_atom.ReadPortal().Get(i);
    for (int j = 0; j < 6; ++j)
    {
      virial[j] += vatom[j];
    }
  }

  //
  //for (int i = 0; i < 6; ++i)
  //{
  //  std::cout << "total_virial[" << i << "] = " << virial[i] << std::endl;
  //}

  //compute pressure scalar
  auto volume = GetParameter<Real>(PARA_VOLUME);
  auto temperature = GetParameter<Real>(PARA_TEMPT);
  auto inv_volume = 1.0 / volume;
  //
  auto n = _position.GetNumberOfValues();
  auto extra_dof = 3; //dimension =3

  auto dof = 3 * n - extra_dof;

  pressure_scalar = (dof * _unit_factor.boltz * temperature + virial[0] + virial[1] + virial[2]) / 3.0 
                      * inv_volume * _unit_factor.nktv2p;

  std::cout << " scalar=" << pressure_scalar << std::endl;
}

void EAMSystem::Compute_Temp_Scalar() {}

void EAMSystem::Couple() 
{
  p_current[0] = p_current[1] = p_current[2] = pressure_scalar;

}

void EAMSystem::end_of_step()
{
  bulkmodulus = 10.0;
  p_start[0] = p_start[1] = p_start[2] = 1.0;
  p_stop[0] = p_stop[1] = p_stop[2] = 1.0;
  p_period[0] = p_period[1] = p_period[2] = 1.0;

  // compute new T,P

  Compute_Pressure_Scalar();
  Couple();

  auto currentstep = _app.GetExecutioner()->CurrentStep();
  auto beginstep = 0;
  auto endstep = _app.GetExecutioner()->NumStep();

  auto delta = currentstep - beginstep;
  std::cout << delta << std::endl;
  if (delta != 0.0)
  {
    delta = delta / (endstep - beginstep);
  }
  for (int i = 0; i < 3; i++)
  {
    p_target[i] = p_start[i] + delta * (p_stop[i] - p_start[i]);
    dilation[i] = pow(1.0 - _dt / p_period[i] * (p_target[i] - p_current[i]) / bulkmodulus, 1.0 / 3.0);
  }

  for (int i = 0; i < 3; i++)
  {
    std::cout << "i=" << i  << ",dilation=" << dilation[i] << std::endl;
  }
  // remap simulation box and atoms
  // redo KSpace coeffs since volume has changed

  remap();


}

void EAMSystem::x2lamda(Id n)
{
  //
  set_global_box();
  //
  Vec3f delta;
  auto position = GetFieldAsArrayHandle<Vec3f>(field::position);
  //
  vtkm::Vec<vtkm::Range, 3> range = GetParameter<vtkm::Vec<vtkm::Range, 3>>(PARA_RANGE);

  for (int i = 0; i < n; i++)
  {
    delta[0] = position.ReadPortal().Get(i)[0] - range[0].Min;
    delta[1] = position.ReadPortal().Get(i)[1] - range[1].Min;
    delta[2] = position.ReadPortal().Get(i)[2] - range[2].Min;

    position.WritePortal().Get(i)[0] = h_inv[0] * delta[0] + h_inv[5] * delta[1] + h_inv[4] * delta[2];
    position.WritePortal().Get(i)[1] = h_inv[1] * delta[1] + h_inv[3] * delta[2];
    position.WritePortal().Get(i)[2] = h_inv[2] * delta[2];
  }
  //for (int i = 0; i < n; i++)
  //{
  //  for (int j = 0; j < 3; j++)
  //  {
  //    std::cout << "i=" << i << ", lamda_position=" << position.ReadPortal().Get(i)[j] << std::endl;
  //  }
  //}
}

void EAMSystem::lamda2x(Id n)
{
  auto position = GetFieldAsArrayHandle<Vec3f>(field::position);
  vtkm::Vec<vtkm::Range, 3> range = GetParameter<vtkm::Vec<vtkm::Range, 3>>(PARA_RANGE);

  for (int i = 0; i < n; i++)
  {
    position.WritePortal().Get(i)[0] = h[0] * position.ReadPortal().Get(i)[0] + 
                                       h[5] * position.ReadPortal().Get(i)[1] +
                                       h[4] * position.ReadPortal().Get(i)[2] + range[0].Min;

    position.WritePortal().Get(i)[1] = h[1] * position.ReadPortal().Get(i)[1] +
                                       h[3] * position.ReadPortal().Get(i)[2] +  
                                       range[1].Min;

    position.WritePortal().Get(i)[2] = h[2] * position.ReadPortal().Get(i)[2] + range[2].Min;
  }
  //for (int i = 0; i < n; i++)
  //{
  //  for (int j = 0; j < 3; j++)
  //  {
  //    std::cout << "i=" << i << ", position=" << position.ReadPortal().Get(i)[j] << std::endl;
  //  }
  //}
}

void EAMSystem::set_global_box()
{
  vtkm::Vec<vtkm::Range, 3> range = GetParameter<vtkm::Vec<vtkm::Range, 3>>(PARA_RANGE);
  prd[0] = xprd = range[0].Max - range[0].Min;
  prd[1] = yprd = range[1].Max - range[1].Min;
  prd[2] = zprd = range[2].Max - range[2].Min;

  h[0] = xprd;
  h[1] = yprd;
  h[2] = zprd;
  //
  auto orthogonal  = 1;
  if (orthogonal)
  {
    h_inv[0] = 1.0 / h[0];
    h_inv[1] = 1.0 / h[1];
    h_inv[2] = 1.0 / h[2];
    h_inv[3] = 0;
    h_inv[4] = 0;
    h_inv[5] = 0;
  }


}

void EAMSystem::remap() 
{
  Real oldlo, oldhi, ctr;
  auto position = GetFieldAsArrayHandle<Vec3f>(field::position);
  auto n = position.GetNumberOfValues();

  vtkm::Vec<vtkm::Range, 3> range = GetParameter<vtkm::Vec<vtkm::Range, 3>>(PARA_RANGE);

  // convert pertinent atoms and rigid bodies to lamda coords
  x2lamda(n);

  // reset global and local box to new size/shape

  for (int i = 0; i < 3; i++)
  {
    oldlo = range[i].Min;
    oldhi = range[i].Max;
    ctr = 0.5 * (oldlo + oldhi);
    range[i].Min = (oldlo - ctr) * dilation[i] + ctr;
    range[i].Max = (oldhi - ctr) * dilation[i] + ctr;
  }

  set_global_box();

  // convert pertinent atoms and rigid bodies back to box coords

  lamda2x(n);

}

//void EAMSystem::ev_tall() 
//{
//  auto Vlength = GetParameter<Real>(PARA_VLENGTH);
//  SystemWorklet::ComputeFpair(_force, fpair);
//  for (int i = 0; i < fpair.GetNumberOfValues(); ++i)
//  {
//    std::cout << "i=" << i << ",fpair=" << fpair.ReadPortal().Get(i) << std::endl;
//  }
//  SystemWorklet::ComputeVirial0(Vlength, _atoms_id, fpair, _locator, _force_function, virial);
//  //for (int i = 0; i < virial.GetNumberOfValues();++i)
//  //{
//  //  std::cout << "i=" << i << ",virial="  << virial.ReadPortal().Get(0) << 
//  //               virial.ReadPortal().Get(1) <<
//  //               virial.ReadPortal().Get(2) <<
//  //               virial.ReadPortal().Get(3) << 
//  //               virial.ReadPortal().Get(4) << 
//  //                  virial.ReadPortal().Get(5)
//  //      << std::endl;
//
//  //}
//
//}

