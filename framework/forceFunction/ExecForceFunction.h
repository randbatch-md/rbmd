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
#include <vtkm/VectorAnalysis.h>
#include <vtkm/cont/ArrayHandle.h>
#include <vtkm/cont/CoordinateSystem.h>
#include <vtkm/Deprecated.h>
#include <vtkm/Pair.h>
#include <vtkm/Math.h>
#include "Types.h"
#include "math/Math.h"
#define MIN(A, B) ((A) < (B) ? (A) : (B))
#define MAX(A, B) ((A) > (B) ? (A) : (B))

class ExecForceFunction
{
public:
  using IdPortalTypeId = typename vtkm::cont::ArrayHandle<vtkm::Id>::ReadPortalType;
  using IdPortalTypeReal = typename vtkm::cont::ArrayHandle<Real>::ReadPortalType;
  using IdPortalTypeVec2 = typename vtkm::cont::ArrayHandle<vtkm::Vec2i>::ReadPortalType;
  using Matrix3x3 = vtkm::Vec<vtkm::Vec3f, 3>; // 定义 3x3 矩阵

  ExecForceFunction(const Real cut_off,
                    const Real& alpha,
                    const Real& volume,
                    const Real& vlength,
                    const Vec3f& box,
                    const IdComponent& kmax,
                    const IdComponent& rbeP,
                    const Real& rhomax,
                    const Id& nrho,
                    const Real& drho,
                    const Id& nr,
                    const Real& dr)
    : _cut_Off(cut_off)
    , _alpha(alpha)
    , _volume(volume)
    , _Vlength(vlength)
    , _box(box)
    , Kmax(kmax)
    , RBEP(rbeP)
    , _rhomax(rhomax)
    , _nrho(nrho)
    , _drho(drho)
    , _nr(nr)
    , _dr(dr)
  { 
      _sum_gauss = Compute_S();
  }
 
  VTKM_EXEC void angmom_to_omega(const Vec3f& m,
      const Vec3f& ex,
      const Vec3f& ey,
      const Vec3f& ez,
      const Vec3f& idiag,
      Vec3f& w) const
  {
      Vec3f wbody;

      // 计算主体坐标系下的角速度 wbody
      wbody[0] = (idiag[0] == 0.0) ? 0.0 : (vtkm::Dot(m, ex) / idiag[0]);
      wbody[1] = (idiag[1] == 0.0) ? 0.0 : (vtkm::Dot(m, ey) / idiag[1]);
      wbody[2] = (idiag[2] == 0.0) ? 0.0 : (vtkm::Dot(m, ez) / idiag[2]);

      // 将主体坐标系下的角速度转换为空间坐标系下的角速度
      w[0] = wbody[0] * ex[0] + wbody[1] * ey[0] + wbody[2] * ez[0];
      w[1] = wbody[0] * ex[1] + wbody[1] * ey[1] + wbody[2] * ez[1];
      w[2] = wbody[0] * ex[2] + wbody[1] * ey[2] + wbody[2] * ez[2];
  }

  VTKM_EXEC  Id jacobi3(Matrix3x3& mat, Vec3f& eval, Matrix3x3& evec) const
  {
      // 1. 复制输入矩阵以便后续修改
      Matrix3x3 mat_cpy = mat;

      // 2. 初始化特征值、特征向量
      eval = Vec3f(0.0f, 0.0f, 0.0f); // 存储特征值
      evec = Matrix3x3();                   // 存储特征向量

      // 3. 调用 Jacobi 迭代法进行对角化
      // max_num_sweeps=50;
      Id ierror = Diagonalize(mat_cpy, eval, evec, true, 50);

      if (ierror != 0) {
          return ierror;  // 如果计算失败，则返回错误代码
      }

      // 4. 转置特征向量矩阵
      for (vtkm::Id i = 0; i < 3; i++) {
          for (vtkm::Id j = i + 1; j < 3; j++) {
              std::swap(evec[i][j], evec[j][i]);
          }
      }

      return 0;  // 成功完成
  }

  // Jacobi 对角化函数实现
  VTKM_EXEC  Id Diagonalize(const Matrix3x3& mat, Vec3f& eval, Matrix3x3& evec, bool calc_evec, int max_num_sweeps) const
  {
      // 初始化：将输入矩阵拷贝到本地矩阵
      Matrix3x3 M = mat;

      // 初始化特征值向量和特征向量矩阵
      eval = Vec3f(0.0, 0.0, 0.0);
      evec = Matrix3x3();  // 初始化为单位矩阵
      if (calc_evec)
      {
          for (vtkm::Id i = 0; i < 3; i++)
          {
              for (vtkm::Id j = 0; j < 3; j++) {
                  evec[i][j] = (i == j) ? 1.0f : 0.0f;
              }
          }
      }

      // 初始化最大值索引
      Id3 max_idx_row = { MaxEntryRow(M, 0), MaxEntryRow(M, 1), MaxEntryRow(M, 2) };

      // 迭代过程
      Id n_iters;
      Id max_num_iters = max_num_sweeps * 3 * (3 - 1) / 2;
      for (n_iters = 0; n_iters < max_num_iters; ++n_iters) {
          Id i, j;
          MaxEntry(M, i, j);  // 找到最大值的位置

          // 检查是否接近零
          if ((M[i][i] + M[i][j] == M[i][i]) && (M[j][j] + M[i][j] == M[j][j])) {
              M[i][j] = 0.0;
              max_idx_row[i] = MaxEntryRow(M, i);
          }

          // 如果 M[i][j] 已经为 0，则停止迭代
          if (M[i][j] == 0.0) {
              break;
          }

          // 否则，计算旋转矩阵并应用旋转
          vtkm::FloatDefault c, s;
          CalcRot(M, i, j, c, s);  // 计算旋转矩阵参数
          ApplyRot(M, i, j, c, s); // 应用旋转到矩阵 M

          // 如果需要计算特征向量，则应用旋转到特征向量矩阵
          if (calc_evec) {
              ApplyRotLeft(evec, i, j, c, s);
          }
      }

      // 提取对角线元素作为特征值
      for (Id i = 0; i < 3; ++i) {
          eval[i] = M[i][i];
      }

      // 对特征值和特征向量进行排序
      SortRows(eval, evec);

      return (n_iters == max_num_iters) ? 1 : 0; // 1 表示达到最大迭代次数，0 表示成功
  }

  // 计算 M 的第 row 行中非对角元素最大值的索引
  VTKM_EXEC  Id MaxEntryRow(const Matrix3x3& M, Id row) const
  {
      vtkm::FloatDefault max_val = 0.0;
      Id max_idx = row;
      for (Id j = 0; j < 3; ++j) {
          if (j != row && vtkm::Abs(M[row][j]) > max_val) {
              max_val = vtkm::Abs(M[row][j]);
              max_idx = j;
          }
      }
      return max_idx;
  }

  // 找到矩阵中最大的非对角元素
  VTKM_EXEC void MaxEntry(const Matrix3x3& M, Id& i, Id& j) const
  {
      vtkm::FloatDefault max_val = 0.0;
      for (Id row = 0; row < 3; ++row) {
          for (Id col = row + 1; col < 3; ++col) {
              vtkm::FloatDefault val = vtkm::Abs(M[row][col]);
              if (val > max_val) {
                  max_val = val;
                  i = row;
                  j = col;
              }
          }
      }
  }

  // 计算旋转矩阵参数 c 和 s
  VTKM_EXEC void CalcRot(Matrix3x3& M, Id i, Id j, vtkm::FloatDefault& c, vtkm::FloatDefault& s) const
  {
      vtkm::FloatDefault theta = (M[j][j] - M[i][i]) / (2.0 * M[i][j]);
      vtkm::FloatDefault t;

      if (theta >= 0.0) {
          t = 1.0 / (theta + vtkm::Sqrt(1.0 + theta * theta));
      }
      else {
          t = -1.0 / (-theta + vtkm::Sqrt(1.0 + theta * theta));
      }

      c = 1.0 / vtkm::Sqrt(1.0 + t * t);
      s = t * c;
  }

  // 应用旋转到矩阵 M
  VTKM_EXEC  void ApplyRot(Matrix3x3& M, Id i, Id j, vtkm::FloatDefault c, vtkm::FloatDefault s) const
  {
      for (Id k = 0; k < 3; ++k) {
          vtkm::FloatDefault Mik = M[i][k];
          vtkm::FloatDefault Mjk = M[j][k];
          M[i][k] = c * Mik - s * Mjk;
          M[j][k] = s * Mik + c * Mjk;
      }
  }

  // 应用旋转到特征向量矩阵 evec
  VTKM_EXEC  void ApplyRotLeft(Matrix3x3& evec, Id i, Id j, vtkm::FloatDefault c, vtkm::FloatDefault s) const
  {
      for (Id k = 0; k < 3; ++k) {
          vtkm::FloatDefault eik = evec[i][k];
          vtkm::FloatDefault ejk = evec[j][k];
          evec[i][k] = c * eik - s * ejk;
          evec[j][k] = s * eik + c * ejk;
      }
  }

  // 对特征值和特征向量排序
  VTKM_EXEC void SortRows(Vec3f& eval, Matrix3x3& evec) const
  {
      for (Id i = 0; i < 3; ++i) {
          for (Id j = i + 1; j < 3; ++j) {
              if (eval[j] > eval[i]) {
                  std::swap(eval[i], eval[j]);
                  for (Id k = 0; k < 3; ++k) {
                      std::swap(evec[i][k], evec[j][k]);
                  }
              }
          }
      }
  }

  VTKM_EXEC Vec6f ComputeLJVirial(const Vec3f& r_ij,
      const Real& eps_i,
      const Real& eps_j,
      const Real& sigma_i,
      const Real& sigma_j,
      const Real& cut_off) const
  {
      Vec6f LJVirial = { 0, 0, 0, 0, 0, 0 };

      const Real small_value = 0.0001;
      Vec3f force{ 0, 0, 0 };
      const Real dis_2 = r_ij[0] * r_ij[0] + r_ij[1] * r_ij[1] + r_ij[2] * r_ij[2];
      const Real cut_off_2 = cut_off * cut_off;

      if (dis_2 < cut_off_2 && dis_2 > small_value)
      {
          Real sigma_ij = (sigma_i + sigma_j) / 2;

          Real sigmaij_6 = sigma_ij * sigma_ij * sigma_ij * sigma_ij * sigma_ij * sigma_ij;
          Real dis_6 = dis_2 * dis_2 * dis_2;
          Real sigmaij_dis_6 = sigmaij_6 / dis_6;
          Real eps_ij = vtkm::Sqrt(eps_i * eps_j);
          auto f = 24 * eps_ij * ((2 * sigmaij_dis_6 - 1) * sigmaij_dis_6) / dis_2;
          auto LJPair = 0.5 * f * r_ij;

          //compute LJVirial
          LJVirial[0] = r_ij[0] * LJPair[0]; //xx
          LJVirial[1] = r_ij[1] * LJPair[1]; //yy
          LJVirial[2] = r_ij[2] * LJPair[2]; //zz
          LJVirial[3] = r_ij[0] * LJPair[1]; //xy
          LJVirial[4] = r_ij[0] * LJPair[2]; //xz
          LJVirial[5] = r_ij[1] * LJPair[2]; //yz
      }
      return LJVirial;
  }

  VTKM_EXEC Vec6f ComputePairVirial(const Vec3f& r_ij, const Vec3f& fpair) const
  {
      Vec6f LJVirial = { 0, 0, 0, 0, 0, 0 };
      auto LJPair = -0.5 * fpair;

      //compute LJVirial
      LJVirial[0] = r_ij[0] * LJPair[0]; //xx
      LJVirial[1] = r_ij[1] * LJPair[1]; //yy
      LJVirial[2] = r_ij[2] * LJPair[2]; //zz
      LJVirial[3] = r_ij[0] * LJPair[1]; //xy
      LJVirial[4] = r_ij[0] * LJPair[2]; //xz
      LJVirial[5] = r_ij[1] * LJPair[2]; //yz
      return LJVirial;
  }

  VTKM_EXEC Vec6f ComputePairVirial_fix(const Vec3f& r_ij, const Vec3f& fpair) const
  {
      Vec6f LJVirial = { 0, 0, 0, 0, 0, 0 };
      auto LJPair = 0.5 * fpair;

      //compute LJVirial
      LJVirial[0] = r_ij[0] * LJPair[0]; //xx
      LJVirial[1] = r_ij[1] * LJPair[1]; //yy
      LJVirial[2] = r_ij[2] * LJPair[2]; //zz
      LJVirial[3] = r_ij[0] * LJPair[1]; //xy
      LJVirial[4] = r_ij[0] * LJPair[2]; //xz
      LJVirial[5] = r_ij[1] * LJPair[2]; //yz
      return LJVirial;
  }

  VTKM_EXEC Vec6f ComputeCoulVirial(const Vec3f& r_ij,
      const Real& charge_pi,
      const Real& charge_pj,
      const Real& cut_off) const
  {
      Vec6f CoulVirial{ 0,0,0,0,0,0, };
      Vec3f CoulForce{ 0, 0, 0 };

      const Real small_value = 0.0001;
      Real dis = vtkm::Magnitude(r_ij);
      const Real dis_2 = r_ij[0] * r_ij[0] + r_ij[1] * r_ij[1] + r_ij[2] * r_ij[2];
      const Real cut_off_2 = cut_off * cut_off;

      if (dis_2 < cut_off_2 && dis_2 > small_value)
      {
          auto f = charge_pi * charge_pj * Gnear(_alpha, dis) / dis;
          CoulForce = -0.5 * f * r_ij;

          //compute CoulVirial
          CoulVirial[0] = r_ij[0] * CoulForce[0]; //xx
          CoulVirial[1] = r_ij[1] * CoulForce[1]; //yy
          CoulVirial[2] = r_ij[2] * CoulForce[2]; //zz
          CoulVirial[3] = r_ij[0] * CoulForce[1]; //xy
          CoulVirial[4] = r_ij[0] * CoulForce[2]; //xz
          CoulVirial[5] = r_ij[1] * CoulForce[2]; //yz
      }
      return CoulVirial;
  }

  VTKM_EXEC Vec6f ComputeLongVirial(const Vec3f& M,
      const Vec3f& r_i,
      const Real& charge_p_i,
      const Vec2f& rhok_ri)
  {
      Vec6f LongVirial{ 0, 0, 0, 0, 0, 0 };

      Vec3f K{ 0, 0, 0 };
      Real LongForce = 0;

      K[0] = 2 * vtkm::Pi() * M[0] / _box[0];
      K[1] = 2 * vtkm::Pi() * M[1] / _box[1];
      K[2] = 2 * vtkm::Pi() * M[2] / _box[2];

      Real range_K_2 = K[0] * K[0] + K[1] * K[1] + K[2] * K[2];
      auto factor_a = 4 * vtkm::Pi() * charge_p_i;
      auto factor_b = vtkm::Exp(-range_K_2 / (4 * _alpha));
      auto factor_c = vtkm::Cos(vtkm::Dot(K, r_i)) * rhok_ri[1];
      auto factor_d = vtkm::Sin(vtkm::Dot(K, r_i)) * rhok_ri[0];

      auto f = factor_a / (_volume * range_K_2) * factor_b * (factor_c - factor_d);
      LongForce = 0.5 * f;

      //compute LongVirial
      LongVirial[0] = (r_i[0] * K[0] + r_i[0] * K[0]) * LongForce; //_xx
      LongVirial[1] = (r_i[1] * K[1] + r_i[1] * K[1]) * LongForce; // yy
      LongVirial[2] = (r_i[2] * K[2] + r_i[2] * K[2]) * LongForce; // zz
      LongVirial[3] = (r_i[0] * K[1] + r_i[1] * K[0]) * LongForce; // xy
      LongVirial[4] = (r_i[0] * K[2] + r_i[2] * K[0]) * LongForce; // xz
      LongVirial[5] = (r_i[1] * K[2] + r_i[2] * K[1]) * LongForce; // yz

      return LongVirial;
  }

  VTKM_EXEC Vec6f ComputeKspaceVirial(const Vec3f& K,
      const Vec3f& r_i,
      const Real& KspaceForce)
  {
      Vec6f KspaceVirial{ 0, 0, 0, 0, 0, 0 };

      //compute LongVirial
      KspaceVirial[0] =  (r_i[0] * K[0] + r_i[0] * K[0]) * KspaceForce; //_xx
      KspaceVirial[1] =  (r_i[1] * K[1] + r_i[1] * K[1]) * KspaceForce; // yy
      KspaceVirial[2] =  (r_i[2] * K[2] + r_i[2] * K[2]) * KspaceForce; // zz
      KspaceVirial[3] =  (r_i[0] * K[1] + r_i[1] * K[0]) * KspaceForce; // xy
      KspaceVirial[4] =  (r_i[0] * K[2] + r_i[2] * K[0]) * KspaceForce; // xz
      KspaceVirial[5] =  (r_i[1] * K[2] + r_i[2] * K[1]) * KspaceForce; // yz

      return KspaceVirial;
  }

  VTKM_EXEC vtkm::Vec3f ComputeLJForce(const Vec3f& r_ij,
                                       const Real& eps_i,
                                       const Real& eps_j,
                                       const Real& sigma_i,
                                       const Real& sigma_j,
                                       const Real& cut_off) const
  {
    const Real small_value = 0.0001;
    Vec3f force{ 0, 0, 0 };
    const Real dis_2 = r_ij[0] * r_ij[0] + r_ij[1] * r_ij[1] + r_ij[2] * r_ij[2];
    const Real cut_off_2 = cut_off * cut_off;

    //if (dis_2 < cut_off_2 && dis_2 > small_value)
    if (dis_2 < cut_off_2)
    {
      Real sigma_ij = (sigma_i + sigma_j) / 2;

      Real sigmaij_6 = sigma_ij * sigma_ij * sigma_ij * sigma_ij * sigma_ij * sigma_ij;
      Real dis_6 = dis_2 * dis_2 * dis_2;
      Real sigmaij_dis_6 = sigmaij_6 / dis_6;
      force = -24 * vtkm::Sqrt(eps_i * eps_j) * ((2 * sigmaij_dis_6 - 1) * sigmaij_dis_6) / dis_2 * r_ij;
    }
    return force;
  }

  VTKM_EXEC vtkm::Vec3f ComputeLJForceRcs(const Vec3f& r_ij,
                                          const Real& eps_i,
                                          const Real& eps_j,
                                          const Real& sigma_i,
                                          const Real& sigma_j,
                                          const Real& cut_off,
                                          const Real& rs) const
  {
    Vec3f force{ 0, 0, 0 };
    const Real dis_2 = r_ij[0] * r_ij[0] + r_ij[1] * r_ij[1] + r_ij[2] * r_ij[2];

    const Real cut_off_2 = cut_off * cut_off;
    const Real rs_2 = rs * rs;

    if (dis_2 < cut_off_2 && dis_2 > rs_2)
    {
      Real sigma_ij = (sigma_i + sigma_j) / 2;

      Real sigmaij_6 = sigma_ij * sigma_ij * sigma_ij * sigma_ij * sigma_ij * sigma_ij;
      Real dis_6 = dis_2 * dis_2 * dis_2;
      Real sigmaij_dis_6 = sigmaij_6 / dis_6;

      force = -24 * vtkm::Sqrt(eps_i * eps_j) * ((2 * sigmaij_dis_6 - 1) * sigmaij_dis_6) / dis_2 * r_ij;
    }
    return force;
  }

  VTKM_EXEC Vec3f ComputeNearEnergyForce(const Vec3f& p_i,
                                         const Vec3f& p_j,
                                         const Real& charge_pi,
                                         const Real& charge_pj) const
  {
    const Real small_value = 0.01;
    Vec3f force{ 0, 0, 0 };
    auto r_ij = p_j - p_i;
    Real dis = vtkm::Magnitude(r_ij);

    if (dis < _cut_Off && dis > small_value)
    {
      force = -charge_pi * charge_pj * Gnear(_alpha, dis) * r_ij / dis;
    }
    return force;
  }

  VTKM_EXEC Vec3f ComputeNearEnergyForce1(const Vec3f& r_ij,
      const Real& charge_pi,
      const Real& charge_pj,
      const Real& cut_off) const
  {
      const Real small_value = 0.0001;
      Vec3f force{ 0, 0, 0 };

      Real dis = vtkm::Magnitude(r_ij);
      const Real dis_2 = r_ij[0] * r_ij[0] + r_ij[1] * r_ij[1] + r_ij[2] * r_ij[2];
      const Real cut_off_2 = cut_off * cut_off;

      if (dis_2 < cut_off_2 && dis_2 > small_value)
      {
          force = -charge_pi * charge_pj * Gnear(_alpha, dis) * r_ij / dis;
      }
      return force;
  }

  VTKM_EXEC void  ComputeNearEnergyForce_fix(const Vec3f& r_ij,
      const Real& charge_pi,
      const Real& charge_pj,
      const Real& cut_off,
      Vec3f& force_coul_factor,
      Vec3f& force_coul) const
  {
      const Real small_value = 0.0001;
      const Real dis_2 = r_ij[0] * r_ij[0] + r_ij[1] * r_ij[1] + r_ij[2] * r_ij[2];
      const Real dis = vtkm::Sqrt(dis_2);
      const Real dis_3 = dis * dis * dis;
      const Real cut_off_2 = cut_off * cut_off;

      if (dis_2 < cut_off_2 && dis_2 > small_value)
      {
          force_coul_factor = -charge_pi * charge_pj * r_ij / dis_3; //-
          force_coul = -charge_pi * charge_pj * Gnear(_alpha, dis) * r_ij / dis;
      }
  }

  VTKM_EXEC Vec3f ComputeNearEnergyForceRcs(const Vec3f& p_i,
                                         const Vec3f& p_j,
                                         const Real& charge_pi,
                                         const Real& charge_pj,
                                         const Real& rs) const
  {
    const Real small_value = 0.01;
    Vec3f force{ 0, 0, 0 };
    auto r_ij = p_j - p_i;
    Real dis = vtkm::Magnitude(r_ij);

    if (dis < _cut_Off && dis > rs)
    {
      force = -charge_pi * charge_pj * Gnear(_alpha, dis) * r_ij / dis;
    }
    return force;
  }

  VTKM_EXEC Vec3f ComputeNearEnergyForceERF(const Vec3f& r_ij,
                                         const Real& charge_pi,
                                         const Real& charge_pj,
                                         const Real& table_pij) const
  {
    const Real small_value = 0.0001;
    Vec3f force{ 0, 0, 0 };
    Real dis = vtkm::Magnitude(r_ij);
    Real dis_2 = vtkm::MagnitudeSquared(r_ij);
    auto rc_2 = _cut_Off * _cut_Off;
    if (dis_2 < rc_2 && dis_2 > small_value)
    {
      force = -charge_pi * charge_pj * table_pij * r_ij / dis;
    }
    return force;
  }

  VTKM_EXEC void  ComputeNearEnergyForceRsERF_fix(
      const Vec3f& r_ij,
      const Real& charge_pi,
      const Real& charge_pj,
      const Real& table_pij,
      const Real& rs,
      Vec3f& force_coul_factor,
      Vec3f& force_coul) const
  {
      const Real small_value = 0.0001;
      const Real dis_2 = r_ij[0] * r_ij[0] + r_ij[1] * r_ij[1] + r_ij[2] * r_ij[2];
      const Real dis = vtkm::Sqrt(dis_2);
      const Real dis_3 = dis * dis * dis;
      auto rs_2 = rs * rs;
      //if (dis_2 < rs_2 && dis_2 > small_value)
      if (dis_2 < rs_2)
      {
          force_coul_factor = -charge_pi * charge_pj * r_ij / dis_3;//+
          force_coul = -charge_pi * charge_pj * table_pij * r_ij / dis;
      }
  }

  VTKM_EXEC Vec3f ComputeNearEnergyForceRcsERF(const Vec3f& r_ij,
                                          const Real& charge_pi,
                                          const Real& charge_pj,
                                          const Real& table_pij,
                                          const Real& rs) const
  {
    const Real small_value = 0.0001;
    Vec3f force{ 0, 0, 0 };

    Real dis = vtkm::Magnitude(r_ij);
    Real dis_2 = vtkm::MagnitudeSquared(r_ij);
    auto rc_2 = _cut_Off * _cut_Off;
    auto rs_2 = rs * rs;
    if (dis_2 < rc_2 && dis_2 > rs_2)
    {
      force = -charge_pi * charge_pj * table_pij * r_ij / dis;
    }
    return force;
  }


  VTKM_EXEC void  ComputeNearEnergyForceRcsERF_fix(
      const Vec3f& r_ij,
      const Real& charge_pi,
      const Real& charge_pj,
      const Real& table_pij,
      const Real& rc,
      const Real& rs,
      Vec3f& force_coul_factor,
      Vec3f& force_coul) const
  {
      const Real small_value = 0.0001;
      const Real dis_2 = r_ij[0] * r_ij[0] + r_ij[1] * r_ij[1] + r_ij[2] * r_ij[2];
      const Real dis = vtkm::Sqrt(dis_2);
      const Real dis_3 = dis * dis * dis;

      auto rc_2 = rc * rc;
      auto rs_2 = rs * rs;
      if (dis_2 < rc_2 && dis_2 > rs_2)
      {
          force_coul_factor = -charge_pi * charge_pj * r_ij / dis_3;//+
          force_coul = -charge_pi * charge_pj * table_pij * r_ij / dis;
      }
  }

  VTKM_EXEC Vec3f ComputeNearEnergyForceERF_box(const Vec3f& r_ij,
                                                const Real& charge_pi,
                                                const Real& charge_pj,
                                                const Real& table_pij) const
  {
    const Real small_value = 0.0001;
    Vec3f force{ 0, 0, 0 };
    Real dis = vtkm::Magnitude(r_ij);
    Real dis_2 = vtkm::MagnitudeSquared(r_ij);
    auto rc_2 = _cut_Off * _cut_Off;
    if (dis_2 < rc_2 && dis_2 > small_value)
    {
      force = -charge_pi * charge_pj * table_pij * r_ij / dis;
    }
    return force;
  }

  VTKM_EXEC Real Gnear(const Real& _nearalpha, const Real& dis) const
  {
    Real erfcx = vtkm::Sqrt(_nearalpha) * dis;
    Real expx = -_nearalpha * dis * dis;
    Real Gnearvalue = (1.0 - vtkm::ERF(erfcx)) / (dis * dis) +
      2 * vtkm::Sqrt(_nearalpha) * vtkm::Exp(expx) / (vtkm::Sqrt(vtkm::Pif()) * dis);
    return Gnearvalue;
  }

  VTKM_EXEC Vec3f ComputeFarEle(const Vec3f& M,
                                const Vec3f& r_i,
                                const Real& charge_p_i,
                                const Vec2f& rhok_ri)
  {
    Vec3f force = { 0, 0, 0 };
    Vec3f K{ 0, 0, 0 };
    K[0] = 2 * vtkm::Pi() * M[0] / _box[0];
    K[1] = 2 * vtkm::Pi() * M[1] / _box[1];
    K[2] = 2 * vtkm::Pi() * M[2] / _box[2];
    //Vec3f K = 2 * vtkm::Pi() * M / _Vlength; // TODO: Lx
    Real range_K_2 = K[0] * K[0] + K[1] * K[1] + K[2] * K[2];
    auto factor_a = -4 * vtkm::Pi() * charge_p_i * K;
    auto factor_b = vtkm::Exp(-range_K_2 / (4 * _alpha));
    auto factor_c = vtkm::Cos(vtkm::Dot(K, r_i)) * rhok_ri[1];
    auto factor_d = vtkm::Sin(vtkm::Dot(K, r_i)) * rhok_ri[0];

    force = factor_a / (_volume * range_K_2) * factor_b * (factor_c - factor_d);  
    return force;
  }

  VTKM_EXEC Real Computekspace_fix(const Vec3f& K,
      const Vec3f& r_i,
      const Real& charge_p_i,
      const Vec2f& rhok_ri)
  {
      Real force = 0;
      Real range_K_2 = K[0] * K[0] + K[1] * K[1] + K[2] * K[2];

      auto factor_a = -4 * vtkm::Pi() * charge_p_i;
      auto factor_b = vtkm::Exp(-range_K_2 / (4 * _alpha));
      auto factor_c = vtkm::Cos(vtkm::Dot(K, r_i)) * rhok_ri[1];
      auto factor_d = vtkm::Sin(vtkm::Dot(K, r_i)) * rhok_ri[0];

      force = factor_a / (_volume * range_K_2) * factor_b * (factor_c - factor_d);
      return force;
  }

  VTKM_EXEC Vec3f ComputeRBEForceSum(Vec3f& kl,
                                     const Vec3f& current_pts,
                                     const Real& current_charge,
                                     const Vec2f& rhok_ri)
  {
    Vec3f force = { 0, 0, 0 };
    kl[0] = 2.0 * vtkm::Pi() * kl[0] / _box[0];
    kl[1] = 2.0 * vtkm::Pi() * kl[1] / _box[1];
    kl[2] = 2.0 * vtkm::Pi() * kl[2] / _box[2];
    Real range_kl_2 = kl[0] * kl[0] + kl[1] * kl[1] + kl[2] * kl[2];

    auto factor_a = -4 * vtkm::Pi() * current_charge * kl;
    auto factor_b = vtkm::Cos(vtkm::Dot(kl, current_pts)) * rhok_ri[1];
    auto factor_c = vtkm::Sin(vtkm::Dot(kl, current_pts)) * rhok_ri[0];

    force = (factor_a / (_volume * range_kl_2)) * (factor_b - factor_c);
    return force;
  }

  VTKM_EXEC Real ComputeRBEForceSum_fix(Vec3f& kl,
      const Vec3f& current_pts,
      const Real& current_charge,
      const Vec2f& rhok_ri)
  {
      Real force = 0;
      Real range_kl_2 = kl[0] * kl[0] + kl[1] * kl[1] + kl[2] * kl[2];

      auto factor_a = -4 * vtkm::Pi() * current_charge;
      auto factor_b = vtkm::Cos(vtkm::Dot(kl, current_pts)) * rhok_ri[1];
      auto factor_c = vtkm::Sin(vtkm::Dot(kl, current_pts)) * rhok_ri[0];

      force = (factor_a / (_volume * range_kl_2)) * (factor_b - factor_c);
      return force;
  }

  VTKM_EXEC Vec3f ComputeRBEForce(const Id& p, Vec3f& force)
  {
    return force * _sum_gauss / p;
  }

  VTKM_EXEC Vec6f ComputeRBEVirial_fix(const Id& p, Vec6f& virial)
  {
      return virial * _sum_gauss / p;
  }

  VTKM_EXEC Real Compute_S() const
  {
    const Vec3f& factor = Compute_H();
    Real factor_3 = factor[0] * factor[1] * factor[2];
    Real S = factor_3 - 1;
    return S;
  }

  VTKM_EXEC Vec3f Compute_H() const
  {
    Vec3f H = { 0, 0, 0 };
    for (Id i = 0; i < 3; ++i)
    {
      const Real factor = -(_alpha * _box[i] * _box[i]);
      for (int m = -10; m <= 10; m++)
      {
        Real expx = m * m * factor;
        H[i] += vtkm::Exp(expx);
      }
      H[i] *= vtkm::Sqrt(-(factor) / vtkm::Pi());
    }

    return H;
  }

  VTKM_EXEC Real ComputePotentialEn(const Vec3f& r_ij,
                                    const Real& eps_ij,
                                    const Real& sigma_ij,
                                    const Real& cut_off) const
  {
    const Real small_value = 0.0001;
    Real cut_off_2 = cut_off * cut_off;
    Real ComptePE_ij = 0;
    Real dis_2 = r_ij[0] * r_ij[0] + r_ij[1] * r_ij[1] + r_ij[2] * r_ij[2];

    if (dis_2 < cut_off_2 && dis_2 > small_value)
    {
      Real sigmaij_6 = sigma_ij * sigma_ij * sigma_ij * sigma_ij * sigma_ij * sigma_ij;
      Real dis_6 = dis_2 * dis_2 * dis_2;
      ComptePE_ij = 4 * eps_ij * (sigmaij_6 / dis_6 - 1) * (sigmaij_6 / dis_6);
    }
    return 0.5 * ComptePE_ij;
  }

  VTKM_EXEC Real ComputePotentialEn0(const Vec3f& r_ij,
      const Real& eps_i,
      const Real& eps_j,
      const Real& sigma_i,
      const Real& sigma_j,
      const Real& cut_off) const
  {
      const Real small_value = 0.0001;
      Real cut_off_2 = cut_off * cut_off;
      Real ComptePE_ij = 0;
      Real dis_2 = r_ij[0] * r_ij[0] + r_ij[1] * r_ij[1] + r_ij[2] * r_ij[2];

      if (dis_2 < cut_off_2 && dis_2 > small_value)
      {
          Real sigma_ij = (sigma_i + sigma_j) / 2;
          Real sigmaij_6 = sigma_ij * sigma_ij * sigma_ij * sigma_ij * sigma_ij * sigma_ij;
          Real dis_6 = dis_2 * dis_2 * dis_2;

          Real eps_ij = vtkm::Sqrt(eps_i * eps_j);

          ComptePE_ij = 4 * eps_ij * (sigmaij_6 / dis_6 - 1) * (sigmaij_6 / dis_6);
      }
      return 0.5 * ComptePE_ij;
  }

  VTKM_EXEC Real ComputeNearEleEnergy(const Vec3f& r_ij,
                                    const Real& charge_i,
                                    const Real& charge_j,
                                    const Real& cut_off,
                                    const Real& nearalpha) const
  {
    const Real small_value = 0.01;
    Real ComputePE_ij = 0;
    Real dis = vtkm::Magnitude(r_ij);
    if (dis < cut_off && dis > small_value)
    {
      ComputePE_ij = charge_i * charge_j * (1.0 - vtkm::ERF(vtkm::Sqrt(nearalpha) * dis)) / dis;
    }
    return 0.5 * ComputePE_ij;
  }

  VTKM_EXEC Real  ComputeNearEleEnergy_fix(const Vec3f& r_ij,
      const Real& charge_i,
      const Real& charge_j,
      const Real& cut_off,
      const Real& nearalpha) const
  {
      const Real small_value = 0.0001;
      Real energy_coul = 0;

      Real cut_off_2 = cut_off * cut_off;
      Real dis_2 = r_ij[0] * r_ij[0] + r_ij[1] * r_ij[1] + r_ij[2] * r_ij[2];
      Real dis = vtkm::Sqrt(dis_2);
      if (dis_2 < cut_off_2 && dis_2 > small_value)
      {
          energy_coul = 0.5 * charge_i * charge_j * (1.0 - vtkm::ERF(vtkm::Sqrt(nearalpha) * dis)) / dis;
      }
      return energy_coul;
  }

  VTKM_EXEC Real  ComputeCoulFactor_fix(const Vec3f& r_ij,
      const Real& charge_i,
      const Real& charge_j,
      const Real& cut_off) const
  {
      const Real small_value = 0.0001;
      Real energy_coul_factor = 0;

      Real cut_off_2 = cut_off * cut_off;
      Real dis_2 = r_ij[0] * r_ij[0] + r_ij[1] * r_ij[1] + r_ij[2] * r_ij[2];
      Real dis = vtkm::Sqrt(dis_2);
      if (dis_2 < cut_off_2 && dis_2 > small_value)
      {
          energy_coul_factor = 0.5 * charge_i * charge_j / dis;
      }
      return energy_coul_factor;
  }

    //compute EAMforce;
  template<typename rhor_splineType>
  VTKM_EXEC Real ComputeEAMrho(const Real& rc,
                               const Vec3f& r_ij,
                               const rhor_splineType& rhor_spline) const
  {
    Real EAM_rho = 0;
    auto rdr = 1 / _dr;

    auto rsq = r_ij[0] * r_ij[0] + r_ij[1] * r_ij[1] + r_ij[2] * r_ij[2];
    auto cutsq = rc * rc;

    if (rsq < cutsq && rsq > 0.01)
    {
      auto p = vtkm::Sqrt(rsq) * rdr + 1.0;
      auto m = static_cast<Id>(p);
      // Id m = p;
      m = MIN(m, _nr - 1);
      p -= m;
      p = MIN(p, 1.0);
      auto coeff = rhor_spline.Get(m);
      EAM_rho = ((coeff[3] * p + coeff[4]) * p + coeff[5]) * p + coeff[6];
    }
    return EAM_rho;
  }

  template<typename rhor_splineType>
  VTKM_EXEC Real ComputeEAMrhoRs(const Real& rc,
                                 const Real& rs,
                                 const Vec3f& r_ij,
                                 const rhor_splineType& rhor_spline) const
  {
    Real EAM_rho = 0;
    auto rdr = 1 / _dr;

    auto rsq = r_ij[0] * r_ij[0] + r_ij[1] * r_ij[1] + r_ij[2] * r_ij[2];
    auto cutsq = rc * rc;
    auto rs_2 = rs * rs;
    if (rsq < cutsq && rsq > rs_2)
    {
      auto p = vtkm::Sqrt(rsq) * rdr + 1.0;
      auto m = static_cast<Id>(p);
      // Id m = p;
      m = MIN(m, _nr - 1);
      p -= m;
      p = MIN(p, 1.0);
      auto coeff = rhor_spline.Get(m);
      EAM_rho = ((coeff[3] * p + coeff[4]) * p + coeff[5]) * p + coeff[6];
    }
    return EAM_rho;
  }

  template<typename rhoType, typename frho_splineType>
  VTKM_EXEC Real ComputeEAMfp(const Id& atoms_id,
                              const rhoType& EAM_rho,
                              const frho_splineType& frho_spline) const
  {
    Real fp = 0;
    auto rdrho = 1 / _drho;
    //auto p = EAM_rho.Get(atoms_id) * rdrho + 1.0;
    auto p = EAM_rho * rdrho + 1.0;
    auto m = static_cast<int>(p);
    //Id m = p;
    m = MAX(1, MIN(m, _nrho - 1));
    p -= m;
    p = MIN(p, 1.0);
    auto coeff = frho_spline.Get(m);
    fp = (coeff[0] * p + coeff[1]) * p + coeff[2];
    return fp;
  }

  template<typename rhoType, typename frho_splineType>
  VTKM_EXEC Real ComputeEAMfpOriginal(const Id& atoms_id,
                              const rhoType& EAM_rho,
                              const frho_splineType& frho_spline) const
  {
    Real fp = 0;
    auto rdrho = 1 / _drho;
    auto p = EAM_rho.Get(atoms_id) * rdrho + 1.0;
    //auto p = EAM_rho * rdrho + 1.0;
    auto m = static_cast<int>(p);
    //Id m = p;
    m = MAX(1, MIN(m, _nrho - 1));
    p -= m;
    p = MIN(p, 1.0);
    auto coeff = frho_spline.Get(m);
    fp = (coeff[0] * p + coeff[1]) * p + coeff[2];
    return fp;
  }


  template<typename fpType, typename rhor_splineType, typename z2r_splineType>
  VTKM_EXEC Vec3f ComputeEAMforce(const Real& eam_cut_off,
                                  const Id& atoms_id,
                                  const Id& pts_id_j,
                                  const Vec3f& r_ij,
                                  const fpType fp,
                                  const rhor_splineType& rhor_spline,
                                  const z2r_splineType& z2r_spline) const
  {
    Vec3f eam_force = { 0, 0, 0 };

    auto rdr = 1 / _dr;
    auto rsq = r_ij[0] * r_ij[0] + r_ij[1] * r_ij[1] + r_ij[2] * r_ij[2];
    auto cutsq = eam_cut_off * eam_cut_off;

    if (rsq < cutsq && rsq > 0.01)
    {
      auto r = vtkm::Sqrt(rsq);
      auto p = r * rdr + 1.0;
      auto m = static_cast<int>(p);
      //Id m = p;
      m = MIN(m, _nr - 1);
      p -= m;
      p = MIN(p, 1.0);

      // rhoip = derivative of (density at atom j due to atom i)
      // rhojp = derivative of (density at atom i due to atom j)
      // phi = pair potential energy
      // phip = phi'
      // z2 = phi * r
      // z2p = (phi * r)' = (phi' r) + phi
      // psip needs both fp[i] and fp[j] terms since r_ij appears in two
      //   terms of embed eng: Fi(sum rho_ij) and Fj(sum rho_ji)
      //   hence embed' = Fi(sum rho_ij) rhojp + Fj(sum rho_ji) rhoip
      // scale factor can be applied by thermodynamic integration

      auto coeffi = rhor_spline.Get(m);
      auto rhoip = (coeffi[0] * p + coeffi[1]) * p + coeffi[2];

      auto coeffj = rhor_spline.Get(m);
      auto rhojp = (coeffj[0] * p + coeffj[1]) * p + coeffj[2];

      auto coeff = z2r_spline.Get(m);
      auto z2p = (coeff[0] * p + coeff[1]) * p + coeff[2];
      auto z2 = ((coeff[3] * p + coeff[4]) * p + coeff[5]) * p + coeff[6];

      auto recip = 1.0 / r;
      auto phi = z2 * recip;                 //pair potential energy
      auto phip = z2p * recip - phi * recip; //pair force
      auto psip = fp.Get(atoms_id) * rhojp + fp.Get(pts_id_j) * rhoip + phip;
      auto fpair = -psip * recip;
      //std::cout << fpair << "," << r_ij[0] << "," << r_ij[1] << ","  << r_ij[2] << std::endl;

      //compute f
      eam_force[0] = r_ij[0] * fpair;
      eam_force[1] = r_ij[1] * fpair;
      eam_force[2] = r_ij[2] * fpair;
    }
    return eam_force;
  }

  template<typename fpType, typename rhor_splineType, typename z2r_splineType>
  VTKM_EXEC Vec3f ComputeEAMforceRBL(const Real& eam_cut_off,
                                     const Real& rs,
                                     const Id& atoms_id,
                                     const Id& pts_id_j,
                                     const Vec3f& r_ij,
                                     const fpType fp,
                                     const rhor_splineType& rhor_spline,
                                     const z2r_splineType& z2r_spline) const
  {
    Vec3f eam_force = { 0, 0, 0 };

    auto rdr = 1 / _dr;
    auto rsq = r_ij[0] * r_ij[0] + r_ij[1] * r_ij[1] + r_ij[2] * r_ij[2];
    auto cutsq = eam_cut_off * eam_cut_off;
    auto rs_2 = rs * rs;
    if (rsq < cutsq && rsq > rs_2)
    {
      auto r = vtkm::Sqrt(rsq);
      auto p = r * rdr + 1.0;
      auto m = static_cast<int>(p);
      //Id m = p;
      m = MIN(m, _nr - 1);
      p -= m;
      p = MIN(p, 1.0);

      auto coeffi = rhor_spline.Get(m);
      auto rhoip = (coeffi[0] * p + coeffi[1]) * p + coeffi[2];

      auto coeffj = rhor_spline.Get(m);
      auto rhojp = (coeffj[0] * p + coeffj[1]) * p + coeffj[2];

      auto coeff = z2r_spline.Get(m);
      auto z2p = (coeff[0] * p + coeff[1]) * p + coeff[2];
      auto z2 = ((coeff[3] * p + coeff[4]) * p + coeff[5]) * p + coeff[6];

      auto recip = 1.0 / r;
      auto phi = z2 * recip;                 //pair potential energy
      auto phip = z2p * recip - phi * recip; //pair force
      auto psip = fp.Get(atoms_id) * rhojp + fp.Get(pts_id_j) * rhoip + phip;
      auto fpair = -psip * recip;
      //std::cout << fpair << "," << r_ij[0] << "," << r_ij[1] << ","  << r_ij[2] << std::endl;

      //compute f
      eam_force[0] = r_ij[0] * fpair;
      eam_force[1] = r_ij[1] * fpair;
      eam_force[2] = r_ij[2] * fpair;
    }
    return eam_force;
  }


  //compute EAM Energy = EmbeddingEnergy + PairEnergy;
  template<typename rhor_splineType>
  VTKM_EXEC Real ComputeEAMrhoOUT(const Real& eam_cut_off,
                                  const Vec3f& r_ij,
                                  const rhor_splineType& rhor_spline) const
  {
    Real EAM_rho = 0;
    auto rdr = 1 / _dr;

    auto rsq = r_ij[0] * r_ij[0] + r_ij[1] * r_ij[1] + r_ij[2] * r_ij[2];
    auto cutsq = eam_cut_off * eam_cut_off;

    if (rsq < cutsq && rsq > 0.01)
    {
      auto p = vtkm::Sqrt(rsq) * rdr + 1.0;
      auto m = static_cast<Id>(p);
      //Id m = p;
      m = MIN(m, _nr - 1);
      p -= m;
      p = MIN(p, 1.0);
      auto coeff = rhor_spline.Get(m);
      EAM_rho = ((coeff[3] * p + coeff[4]) * p + coeff[5]) * p + coeff[6];
    }
    return EAM_rho;
  }

  template<typename rhoType, typename frho_splineType>
  VTKM_EXEC Real ComputeEmbeddingEnergy(const Id& atoms_id,
                                        const rhoType& EAM_rho,
                                        const frho_splineType& frho_spline) const
  {
    Real phi = 0;

    auto rdrho = 1 / _drho;
    auto p = EAM_rho.Get(atoms_id) * rdrho + 1.0;
    auto m = static_cast<int>(p);
    //Id m = p;
    m = MAX(1, MIN(m, _nrho - 1));
    p -= m;
    p = MIN(p, 1.0);
    auto coeff = frho_spline.Get(m);
    auto fp = (coeff[0] * p + coeff[1]) * p + coeff[2];
    phi = ((coeff[3] * p + coeff[4]) * p + coeff[5]) * p + coeff[6];
    if (EAM_rho.Get(atoms_id) > _rhomax)
    {
      phi += fp * (EAM_rho.Get(atoms_id) - _rhomax);
    }
    return phi;
  }

  template<typename z2r_splineType>
  VTKM_EXEC Real ComputePairEnergy(const Real& eam_cut_off,
                                   const Vec3f& r_ij,
                                   const z2r_splineType& z2r_spline) const
  {
    Real phi = 0;
    auto rdr = 1 / _dr;
    auto rsq = r_ij[0] * r_ij[0] + r_ij[1] * r_ij[1] + r_ij[2] * r_ij[2];
    auto cutsq = eam_cut_off * eam_cut_off;

    if (rsq < cutsq && rsq > 0.01)
    {
      auto r = vtkm::Sqrt(rsq);
      auto p = r * rdr + 1.0;
      auto m = static_cast<int>(p);
      //Id m = p;
      m = MIN(m, _nr - 1);
      p -= m;
      p = MIN(p, 1.0);

      auto coeff = z2r_spline.Get(m);
      auto z2p = (coeff[0] * p + coeff[1]) * p + coeff[2];
      auto z2 = ((coeff[3] * p + coeff[4]) * p + coeff[5]) * p + coeff[6];

      auto recip = 1.0 / r;
      phi = z2 * recip; //pair potential energy
    }
    return 0.5 * phi;
  }

  private:
  Real _cut_Off;
  Real _alpha;
  Real _volume;
  Real _Vlength;
  Vec3f _box;
  Real _sum_gauss;
  IdComponent Kmax;
  IdComponent RBEP;

  Real _rhomax;
  Id _nr;
  Id _nrho;
  Real _dr;
  Real _drho;

};
