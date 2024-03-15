#include "EAMAlloySystem.h"
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


//RegisterObject(EAMAlloySystem);
#define MIN(A, B) ((A) < (B) ? (A) : (B))
#define MAX(A, B) ((A) > (B) ? (A) : (B))

EAMAlloySystem::EAMAlloySystem(const Configuration& cfg)
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

void EAMAlloySystem::Init()
{
  //ReadPotentialFile:
  ReadPotentialFile(_potential_file);
  InitStyle();

  MDSystem::Init();

  //init variable
  InitialCondition();

  ComputeForce(); // Presolve force
}

void EAMAlloySystem::InitialCondition()
{
  MDSystem::InitialCondition();
  _rho = GetParameter<Real>(PARA_RHO);
  _nosehooverxi = 0.0;

  SetParameter(PARA_RDF_RHO, _rho);
  SetCharge();
  PreForce();
}

void EAMAlloySystem::ReadPotentialFile(std::ifstream& input_file)
{
  // 检查文件是否打开成功
  if (!input_file.is_open())
  {
    std::cerr << "Unable to open the file." << std::endl;
  }

  // 跳过前三行
  for (int i = 0; i < 3; ++i)
  {
    input_file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  }

  // 开始读取第四行的值
  std::string elementsA, elementsB;
  input_file >> alloy_file.nelements >> elementsA >> elementsB;

  // 开始读取第五行的值
  input_file >> alloy_file.nrho >> alloy_file.drho >> alloy_file.nr >> alloy_file.dr >>
    alloy_file.cut_off;

  //分配内存
  auto nelements = alloy_file.nelements;
  alloy_file.frho.resize(nelements, std::vector<Real>(alloy_file.nrho + 1));
  alloy_file.rhor.resize(nelements, std::vector<Real>(alloy_file.nr + 1));
  alloy_file.z2r.resize(nelements);
  for (int i = 0; i < nelements; ++i)
  {
    alloy_file.z2r[i].resize(nelements, std::vector<Real>(alloy_file.nr + 1));
  }

  // 第6行
  int number;
  float mass, lattice_parameter;
  std::string lattice_style;

  input_file >> number >> mass >> lattice_parameter >> lattice_style;
  // 读取并保存 Ni_frho 数组

  for (int i = 0; i < alloy_file.nrho; ++i)
  {
    input_file >> alloy_file.frho[0][i];
    //std::cout << "i= " << i << ", alloy_file.frho: " << alloy_file.frho[0][i] << std::endl;
  }

  // 读取并保存  Ni_rhor 数组

  for (int i = 0; i < alloy_file.nrho; ++i)
  {
    input_file >> alloy_file.rhor[0][i];
    //std::cout << "i= " << i << ", alloy_file.rhor: " << alloy_file.rhor[0][i] << std::endl;
  }

  input_file >> number >> mass >> lattice_parameter >> lattice_style;
  // 读取并保存 Cu_frho 数组

  for (int i = 0; i < alloy_file.nrho; ++i)
  {
    input_file >> alloy_file.frho[1][i];
    //std::cout << "i= " << i << ", alloy_file.frho: " << alloy_file.frho[1][i] << std::endl;
  }

  // 读取并保存  Cu_rhor 数组

  for (int i = 0; i < alloy_file.nrho; ++i)
  {
    input_file >> alloy_file.rhor[1][i];
    //std::cout << "i= " << i << ", alloy_file.rhor: " << alloy_file.rhor[1][i] << std::endl;
  }

  // 读取并保存 z2r 数组

  for (int i = 0; i < nelements; ++i)
  {
    for (int j = 0; j <= i; ++j)
    {
      for (int k = 0; k <= alloy_file.nr; ++k)
      {
        input_file >> alloy_file.z2r[i][j][k];
        // std::cout << "alloy_file.z2r[" << i << "][" << j << "][" << k  << "] = " << alloy_file.z2r[i][j][k] << std::endl;
      }
    }
  }

  // 关闭文件
  input_file.close();

  // 打印读取的值
  std::cout << "alloy_file.nrho: " << alloy_file.nrho << ", alloy_file.drho: " << alloy_file.drho
            << ", alloy_file.nr: " << alloy_file.nr << ", alloy_file.dr: " << alloy_file.dr
            << ", alloy_file.cut_off: " << alloy_file.cut_off << std::endl;
}

void EAMAlloySystem::AllocateEAM() 
{
  Id i, j, m, n;

  //auto ntypes = GetParameter<Id>(NTYPES);
  auto ntypes = 2;

  nrho = alloy_file.nrho;
  nr = alloy_file.nr;
  drho = alloy_file.drho;
  dr = alloy_file.dr;
  rhomax = (nrho - 1) * drho;

  // ------------------------------------------------------------------
  // setup frho arrays
  // ------------------------------------------------------------------
  auto nelements = alloy_file.nelements;
  nfrho = nelements + 1;
  alloy_frho.resize(nfrho, std::vector<Real>(nrho + 1));

  // copy each element's frho to global frho
  for (i = 0; i < nelements; i++)
  {
    for (m = 1; m <= nrho; m++)
    {
      alloy_frho[i][m] = alloy_file.frho[i][m];
      //std::cout << i << "," << m <<", "  << alloy_frho[i][m] << std::endl;
    }
  }

  // add extra frho of zeroes for non-EAM types to point to (pair hybrid)
  // this is necessary b/c fp is still computed for non-EAM atoms

  for (m = 1; m <= nrho; m++)
  {
    alloy_frho[nfrho - 1][m] = 0.0;
  }


  for (i = 1; i <= ntypes; i++)
  {
    if (map[i] >= 0)
      type2frho[i] = map[i];
    else
      type2frho[i] = nfrho - 1;
  }



  // ------------------------------------------------------------------
  // setup rhor arrays
  // ------------------------------------------------------------------
  nrhor = nelements;
  alloy_rhor.resize(nrhor, std::vector<Real>(nr + 1));

  // copy each element's rhor to global rhor
  for (i = 0; i < nelements; i++)
  {
    for (j = 1; j <= nr; j++)
    {
      alloy_rhor[i][j] = alloy_file.rhor[i][j];
      //std::cout << i << "," << j << ", " << alloy_rhor[i][j] << std::endl;
    }
  }

  // type2rhor[i][j] = which rhor array (0 to nrhor-1) each type pair maps to
  // for setfl files, I,J mapping only depends on I
  // OK if map = -1 (non-EAM atom in pair hybrid) b/c type2rhor not used

  for (i = 1; i <= ntypes; i++)
  {
    for (j = 1; j <= ntypes; j++)
    {
      type2rhor[i][j] = map[i];
    }
  }


  // ------------------------------------------------------------------
  // setup z2r arrays
  // ------------------------------------------------------------------

  // nz2r = N*(N+1)/2 where N = # of setfl elements

  nz2r = nelements * (nelements + 1) / 2;
  alloy_z2r.resize(nz2r, std::vector<Real>(nr + 1));

  // copy each element pair z2r to global z2r, only for I >= J
  n = 0;
  for (i = 0; i < nelements; i++)
  {
    for (j = 0; j <= i; j++)
    {
      for (m = 1; m <= nr; m++)
      {
        alloy_z2r[n][m] = alloy_file.z2r[i][j][m];
        //std::cout << n << "," << m << ", " << alloy_z2r[n][m] << std::endl;
      }
      n++;
    }
  }

  // type2z2r[i][j] = which z2r array (0 to nz2r-1) each type pair maps to
  // set of z2r arrays only fill lower triangular Nelement matrix
  // value = n = sum over rows of lower-triangular matrix until reach irow,icol
  // swap indices when irow < icol to stay lower triangular
  // if map = -1 (non-EAM atom in pair hybrid):
  //   type2z2r is not used by non-opt
  //   but set type2z2r to 0 since accessed by opt

  int irow, icol;
  for (i = 1; i <= ntypes; i++)
  {
    for (j = 1; j <= ntypes; j++)
    {
      irow = map[i];
      icol = map[j];
      if (irow == -1 || icol == -1)
      {
        type2z2r[i][j] = 0;
        continue;
      }
      if (irow < icol)
      {
        irow = map[j];
        icol = map[i];
      }
      n = 0;
      for (m = 0; m < irow; m++)
      {
        n += m + 1;
      }
      n += icol;
      type2z2r[i][j] = n;
    }
  }
}

void EAMAlloySystem::file2array()
{
  Id i, j, m, n;

  //auto ntypes = GetParameter<Id>(NTYPES);
  auto ntypes = 2;

  nrho = alloy_file.nrho;
  nr = alloy_file.nr;
  drho = alloy_file.drho;
  dr = alloy_file.dr;
  rhomax = (nrho - 1) * drho;

  // ------------------------------------------------------------------
  // setup frho arrays
  // ------------------------------------------------------------------
  auto nelements = alloy_file.nelements;
  nfrho = nelements + 1;
  alloy_frho.resize(nfrho, std::vector<Real>(nrho + 1));

  // copy each element's frho to global frho
  for (i = 0; i < nelements; i++)
  {
    for (m = 1; m <= nrho; m++)
    {
      alloy_frho[i][m] = alloy_file.frho[i][m];
      //std::cout << i << "," << m <<", "  << alloy_frho[i][m] << std::endl;
    }
  }

  // add extra frho of zeroes for non-EAM types to point to (pair hybrid)
  // this is necessary b/c fp is still computed for non-EAM atoms

  for (m = 1; m <= nrho; m++)
  {
    alloy_frho[nfrho - 1][m] = 0.0;
  }


  for (i = 1; i <= ntypes; i++)
  {
    if (map[i] >= 0)
      type2frho[i] = map[i];
    else
      type2frho[i] = nfrho - 1;
  }



  // ------------------------------------------------------------------
  // setup rhor arrays
  // ------------------------------------------------------------------
  nrhor = nelements;
  alloy_rhor.resize(nrhor, std::vector<Real>(nr + 1));

  // copy each element's rhor to global rhor
  for (i = 0; i < nelements; i++)
  {
    for (j = 1; j <= nr; j++)
    {
      alloy_rhor[i][j] = alloy_file.rhor[i][j];
      //std::cout << i << "," << j << ", " << alloy_rhor[i][j] << std::endl;
    }
  }

  // type2rhor[i][j] = which rhor array (0 to nrhor-1) each type pair maps to
  // for setfl files, I,J mapping only depends on I
  // OK if map = -1 (non-EAM atom in pair hybrid) b/c type2rhor not used

  for (i = 1; i <= ntypes; i++)
  {
    for (j = 1; j <= ntypes; j++)
    {
      type2rhor[i][j] = map[i];
    }
  }


  // ------------------------------------------------------------------
  // setup z2r arrays
  // ------------------------------------------------------------------

  // nz2r = N*(N+1)/2 where N = # of setfl elements

  nz2r = nelements * (nelements + 1) / 2;
  alloy_z2r.resize(nz2r, std::vector<Real>(nr + 1));

  // copy each element pair z2r to global z2r, only for I >= J
  n = 0;
  for (i = 0; i < nelements; i++)
  {
    for (j = 0; j <= i; j++)
    {
      for (m = 1; m <= nr; m++)
      {
        alloy_z2r[n][m] = alloy_file.z2r[i][j][m];
        //std::cout << n << "," << m << ", " << alloy_z2r[n][m] << std::endl;
      }
      n++;
    }
  }

  // type2z2r[i][j] = which z2r array (0 to nz2r-1) each type pair maps to
  // set of z2r arrays only fill lower triangular Nelement matrix
  // value = n = sum over rows of lower-triangular matrix until reach irow,icol
  // swap indices when irow < icol to stay lower triangular
  // if map = -1 (non-EAM atom in pair hybrid):
  //   type2z2r is not used by non-opt
  //   but set type2z2r to 0 since accessed by opt

  int irow, icol;
  for (i = 1; i <= ntypes; i++)
  {
    for (j = 1; j <= ntypes; j++)
    {
      irow = map[i];
      icol = map[j];
      if (irow == -1 || icol == -1)
      {
        type2z2r[i][j] = 0;
        continue;
      }
      if (irow < icol)
      {
        irow = map[j];
        icol = map[i];
      }
      n = 0;
      for (m = 0; m < irow; m++)
      {
        n += m + 1;
      }
      n += icol;
      type2z2r[i][j] = n;
    }
  }
}

void EAMAlloySystem::interpolate(Id n, Real delta, std::vector<Real>& f, std::vector<Vec7f>& spline)
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

void EAMAlloySystem::array2spline()
{
  //分配内存

  alloy_frho_spline.resize(nfrho, std::vector<Vec7f>(nrho + 1));
  alloy_rhor_spline.resize(nrhor, std::vector<Vec7f>(nr + 1));
  alloy_z2r_spline.resize(nz2r, std::vector<Vec7f>(nr + 1));


  for (int i = 0; i < nfrho; i++)
  {
    interpolate(nrho, drho, alloy_frho[i], alloy_frho_spline[i]);
  }
  for (int i = 0; i < nrhor; i++)
  {
    interpolate(nrho, dr, alloy_rhor[i], alloy_rhor_spline[i]);
  }
  for (int i = 0; i < nz2r; i++)
  {
    interpolate(nrho, dr, alloy_z2r[i], alloy_z2r_spline[i]);
  }
  //
  //for (int i = 0; i < alloy_frho_spline.size(); ++i)
  //{
  //  for (int j = 0; j < alloy_frho_spline[i].size(); ++j)
  //  {
  //    std::cout << "j=" << j << ", "<<
  //              alloy_frho_spline[i][j][6] << ", " << alloy_frho_spline[i][j][5] << ", "
  //              << alloy_frho_spline[i][j][4] << ", " << alloy_frho_spline[i][j][3] << ", "
  //              << alloy_frho_spline[i][j][2] << ", " << alloy_frho_spline[i][j][1] << ", "
  //              << alloy_frho_spline[i][j][0] << ", " << std::endl;
  //  }
  //}


  //for (int i = 0; i < alloy_rhor_spline.size(); ++i)
  //{
  //  for (int j = 0; j < alloy_rhor_spline[i].size(); ++j)
  //  {

  //    for (int k = 0; k < 7; ++k)
  //    {
  //      std::cout << alloy_rhor_spline[i][j][6] << ", " <<
  //                alloy_rhor_spline[i][j][5] << ", " <<
  //                alloy_rhor_spline[i][j][4] << ", " <<
  //                alloy_rhor_spline[i][j][3] << ", "<<
  //                alloy_rhor_spline[i][j][2] << ", " <<
  //                alloy_rhor_spline[i][j][1] << ", " <<
  //                alloy_rhor_spline[i][j][0] << ", "  << std::endl;
  //    }
  //  }
  //}

  //for (int i = 0; i < alloy_z2r_spline.size(); ++i)
  //{
  //  for (int j = 0; j < alloy_z2r_spline[i].size(); ++j)
  //  {

  //    for (int k = 0; k < 7; ++k)
  //    {
  //      std::cout << alloy_z2r_spline[i][j][6] << ", " << alloy_z2r_spline[i][j][5] << ", "
  //                << alloy_z2r_spline[i][j][4] << ", " << alloy_z2r_spline[i][j][3] << ", "
  //                << alloy_z2r_spline[i][j][2] << ", " << alloy_z2r_spline[i][j][1] << ", "
  //                << alloy_z2r_spline[i][j][0] << ", " << std::endl;
  //    }
  //  }
  //}
}

void EAMAlloySystem::SetEAM()
{
  SetParameter(EAM_PARA_CUTOFF, alloy_file.cut_off);
  SetParameter(EAM_PARA_RHOMAX, rhomax);
  SetParameter(EAM_PARA_NRHO, nrho);
  SetParameter(EAM_PARA_DRHO, drho);
  SetParameter(EAM_PARA_NR, nr);
  SetParameter(EAM_PARA_DR, dr);
  //


  std::vector<Vec7f> flattened_alloy_frho_spline;
  std::vector<vtkm::Id> frho_countVector;
  for (const std::vector<Vec7f>& innerVector : alloy_frho_spline)
  {
    //扁平化
    flattened_alloy_frho_spline.insert(
      flattened_alloy_frho_spline.end(), innerVector.begin(), innerVector.end());
    //计数数组
    frho_countVector.push_back(innerVector.size());
  }


  //for (int i = 0; i < flattened_vector.size(); ++i)
  //{
  //
  //  std::cout << " i=" << i << ","
  //   << flattened_vector[i][0] << ","
  //  <<flattened_vector[i][1] << ","
  //  << flattened_vector[i][2] << ","
  //  << flattened_vector[i][3] << ","
  //  << flattened_vector[i][4] << ","
  //  << flattened_vector[i][5] << ","
  //  << flattened_vector[i][6]
  //      << std::endl;
  //}

  //sourceArray
  vtkm::cont::ArrayHandle<Vec7f> source_alloy_frho_spline =
    vtkm::cont::make_ArrayHandle<Vec7f>(flattened_alloy_frho_spline);
  //for (int i = 0; i < sourceArray.GetNumberOfValues(); ++i)
  //{
  //    std::cout << " i=" << i << "," <<
  //    sourceArray.ReadPortal().Get(i)[0] << "," <<
  //    sourceArray.ReadPortal().Get(i)[1] << "," <<
  //    sourceArray.ReadPortal().Get(i)[2] << "," <<
  //    sourceArray.ReadPortal().Get(i)[3] << "," <<
  //    sourceArray.ReadPortal().Get(i)[4] << "," <<
  //    sourceArray.ReadPortal().Get(i)[5] << "," <<
  //    sourceArray.ReadPortal().Get(i)[6]
  //            << std::endl;
  //}


  //offsetArray
  vtkm::cont::ArrayHandle<vtkm::Id> frho_countArray =
    vtkm::cont::make_ArrayHandle<vtkm::Id>(frho_countVector);

  vtkm::Id sourceArraySize;
  vtkm::cont::ArrayHandle<vtkm::Id> frho_offsetsArray =
    vtkm::cont::ConvertNumComponentsToOffsets(frho_countArray, sourceArraySize);

  //for (int i = 0; i < frho_ffsetsArray.GetNumberOfValues(); ++i)
  //{
  //  std::cout << " i=" << i << "," << frho_ffsetsArray.ReadPortal().Get(i)
  //          << std::endl;
  //}

  std::vector<Vec7f> flattened_alloy_rhor_spline;
  std::vector<vtkm::Id> rhor_countVector;
  for (const std::vector<Vec7f>& innerVector : alloy_rhor_spline)
  {
    //扁平化
    flattened_alloy_rhor_spline.insert(
      flattened_alloy_rhor_spline.end(), innerVector.begin(), innerVector.end());
    //计数数组
    rhor_countVector.push_back(innerVector.size());
  }
  vtkm::cont::ArrayHandle<Vec7f> source_alloy_rhor_spline =
    vtkm::cont::make_ArrayHandle<Vec7f>(flattened_alloy_rhor_spline);


  //for (int i = 0; i < source_alloy_rhor_spline.GetNumberOfValues(); ++i)
  //{
  //    std::cout << " i=" << i << "," <<
  //    source_alloy_rhor_spline.ReadPortal().Get(i)[0] << "," <<
  //    source_alloy_rhor_spline.ReadPortal().Get(i)[1] << "," <<
  //    source_alloy_rhor_spline.ReadPortal().Get(i)[2] << "," <<
  //    source_alloy_rhor_spline.ReadPortal().Get(i)[3] << "," <<
  //    source_alloy_rhor_spline.ReadPortal().Get(i)[4] << "," <<
  //    source_alloy_rhor_spline.ReadPortal().Get(i)[5] << "," <<
  //    source_alloy_rhor_spline.ReadPortal().Get(i)[6]
  //            << std::endl;
  //}

  vtkm::cont::ArrayHandle<vtkm::Id> rhor_countArray =
    vtkm::cont::make_ArrayHandle<vtkm::Id>(rhor_countVector);

  vtkm::Id rhor_sourceArraySize;
  vtkm::cont::ArrayHandle<vtkm::Id> rhor_offsetsArray =
    vtkm::cont::ConvertNumComponentsToOffsets(rhor_countArray, rhor_sourceArraySize);

  //
  std::vector<Vec7f> flattened_alloy_z2r_spline;
  std::vector<vtkm::Id> z2r_countVector;
  for (const std::vector<Vec7f>& innerVector : alloy_z2r_spline)
  {
    //扁平化
    flattened_alloy_z2r_spline.insert(
      flattened_alloy_z2r_spline.end(), innerVector.begin(), innerVector.end());
    //计数数组
    z2r_countVector.push_back(innerVector.size());
  }
  vtkm::cont::ArrayHandle<Vec7f> source_alloy_z2r_spline =
    vtkm::cont::make_ArrayHandle<Vec7f>(flattened_alloy_z2r_spline);

  //for (int i = 0; i < source_alloy_z2r_spline.GetNumberOfValues(); ++i)
  //{
  //    std::cout << " i=" << i << "," <<
  //    source_alloy_z2r_spline.ReadPortal().Get(i)[0] << "," <<
  //    source_alloy_z2r_spline.ReadPortal().Get(i)[1] << "," <<
  //    source_alloy_z2r_spline.ReadPortal().Get(i)[2] << "," <<
  //    source_alloy_z2r_spline.ReadPortal().Get(i)[3] << "," <<
  //    source_alloy_z2r_spline.ReadPortal().Get(i)[4] << "," <<
  //    source_alloy_z2r_spline.ReadPortal().Get(i)[5] << "," <<
  //    source_alloy_z2r_spline.ReadPortal().Get(i)[6]
  //            << std::endl;
  //}



  //
  vtkm::cont::ArrayHandle<vtkm::Id> z2r_countArray =
    vtkm::cont::make_ArrayHandle<vtkm::Id>(z2r_countVector);

  vtkm::Id z2r_sourceArraySize;
  vtkm::cont::ArrayHandle<vtkm::Id> z2r_offsetsArray =
    vtkm::cont::ConvertNumComponentsToOffsets(z2r_countArray, z2r_sourceArraySize);
}

void EAMAlloySystem::InitStyle()
{
  AllocateEAM();
  file2array();
  array2spline();
  SetEAM();
}

void EAMAlloySystem::ComputeForce()
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

void EAMAlloySystem::UpdateVelocity()
{
  try
  {
    SystemWorklet::UpdateVelocity(_dt, _unit_factor._fmt2v, _force, _mass, _velocity);
    ComputeTempe();
    UpdateVelocityByTempConType();
  }
  catch (const std::exception& e)
  {
    std::cout << e.what() << std::endl;
  }
}

void EAMAlloySystem::UpdatePosition()
{
  auto&& position_flag = GetFieldAsArrayHandle<Id3>(field::position_flag);
  SystemWorklet::UpdatePositionFlag(_dt, _velocity, _locator, _position, position_flag);
  _locator.SetPosition(_position);
  SetCenterTargetPositions();
}

void EAMAlloySystem::SetCharge()
{
  _charge = GetFieldAsArrayHandle<Real>(field::charge);
  auto n = _position.GetNumberOfValues();
  _charge.AllocateAndFill(n, 0);
  _topology.SetCharge(_charge);
}

void EAMAlloySystem::ComputeTempe()
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

void EAMAlloySystem::UpdateVelocityByTempConType()
{
  if (_temp_con_type == "NOSE_HOOVER")
  {
    //Nose Hoover
    //Maybe tauT = 20.0 * dt is a good choice for dt = 2e-3 from the test.
    //Because the temperature curve is the smoothest of all test simulations.
    //In fact, 5.0 and 10.0 are also optional.
    //As long as the coefficent is not too large, such as larger than 100 * dt.
    SystemWorklet::UpdateVelocityNoseHoover(
      _dt, _unit_factor._fmt2v, _nosehooverxi, _force, _mass, _velocity);
    Real tauT = 2000 * _dt;
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
    Real dt_divide_taut = 0.01;
    Real coeff_Berendsen = vtkm::Sqrt(1.0 + dt_divide_taut * (_kbT / _tempT - 1.0));
    SystemWorklet::UpdateVelocityRescale(coeff_Berendsen, _velocity);
  }
}

void EAMAlloySystem::PreSolve()
{
  _locator.SetPosition(_position);
  PreForce();
}

void EAMAlloySystem::Solve()
{
  // stage1:
  UpdateVelocity();

  // stage2:
  UpdatePosition();

  // stage3:
  ComputeForce();
  UpdateVelocity();
}

void EAMAlloySystem::PostSolve() {}

void EAMAlloySystem::PreForce()
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

void EAMAlloySystem::InitField()
{
  MDSystem::InitField();
  AddField(field::pts_type, ArrayHandle<Id>{});
  AddField(field::position_flag, ArrayHandle<Id3>{});
  AddField(field::center_position, ArrayHandle<Vec3f>{});
  AddField(field::target_position, ArrayHandle<Vec3f>{});
  AddField(field::epsilon, ArrayHandle<Real>{});
  AddField(field::sigma, ArrayHandle<Real>{});
}

void EAMAlloySystem::SetCenterTargetPositions()
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

void EAMAlloySystem::TimeIntegration() {}

void EAMAlloySystem::SetForceFunction()
{
  InitERF();
  auto rhomax = GetParameter<Real>(EAM_PARA_RHOMAX);
  auto nrho = GetParameter<Id>(EAM_PARA_NRHO);
  auto drho = GetParameter<Real>(EAM_PARA_DRHO);
  auto nr = GetParameter<Id>(EAM_PARA_NR);
  auto dr = GetParameter<Real>(EAM_PARA_DR);

  _force_function.SetEAMParameters(rhomax, nrho, drho, nr, dr);
}

void EAMAlloySystem::SetTopology()
{
  auto pts_type = GetFieldAsArrayHandle<Id>(field::pts_type);
  auto molecule_id = GetFieldAsArrayHandle<Id>(field::molecule_id);
  auto epsilon = GetFieldAsArrayHandle<Real>(field::epsilon);
  auto sigma = GetFieldAsArrayHandle<Real>(field::sigma);
  _topology.SetAtomsType(pts_type);
  _topology.SetMolecularId(molecule_id);
  _topology.SetEpsAndSigma(epsilon, sigma);
}