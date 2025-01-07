//==================================================================================
//  RBMD 2.2.0 is developed for random batch molecular dynamics calculation.
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
#include "Types.h"
#include <random>
#include <numeric>

template<typename T>
T RandomValue(const Real& Min, const Real& Max){

  std::mt19937 gen(std::random_device{}());
  std::uniform_real_distribution<T> dis(Min, Max);

  return dis(gen);
};

namespace MathRandom
{
/**
* @brief：生成给定范围内不重复的随机整型值数组
* @min：随机数最小值
* @max：随机数最大值
* @num：生成随机数的个数
* @out：随机数数组
*/
static std::vector<Id> NoRepeatIntArray(const Id& min,
                                                    const Id& max,
                                                    const Id& num)
{
  // 设置随机数生成器
  std::random_device rd;
  std::mt19937 gen(rd());

  // 生成均匀分布的整数数组
  const int numRandomNumbers = num;
  std::vector<Id> allNumbers(max - min + 1); // 包括0到num_array的所有整数
  std::iota(allNumbers.begin(), allNumbers.end(), 0);

  // 使用 Fisher-Yates 洗牌算法打乱数组
  std::shuffle(allNumbers.begin(), allNumbers.end(), gen);

  // 从洗牌后的数组中取前numRandomNumbers个整数
  std::vector<Id> randomNumbers(allNumbers.begin(), allNumbers.begin() + numRandomNumbers);

  //从小到达排序
  std::sort(randomNumbers.begin(), randomNumbers.end());

  //生成arrayhandle
  return randomNumbers;

  // 打印生成的随机数
  //std::cout << "Random integers uniformly distributed between 0 and 1000 (non-repeating):"
  //          << std::endl;
  //for (int number : randomNumbers)
  //{
  //  std::cout << number << " ";
  //}
  //std::cout << std::endl;
}

}