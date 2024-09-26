//==================================================================================
//  RBMD 1.0 is developed for random batch molecular dynamics calculation.
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
#include "math/Math.h"
#include <vtkm/Deprecated.h>
#include <vtkm/Math.h>
#include <vtkm/Pair.h>
#include <vtkm/VectorAnalysis.h>
#include <vtkm/cont/ArrayHandle.h>
#include <vtkm/cont/CoordinateSystem.h>
#include <vtkm/cont/Timer.h>
//float exp=3;man=8
const int MANTISSA_BITS = 8;
const int MASK = 67076096; // 2^26-2^15

//float exp=3,man=9;
//const IdComponent MANTISSA_BITS = 9;
//const IdComponent MASK = 67092480; // 2^26-2^14

//float exp=3;man=10
//const int MANTISSA_BITS = 10;
//const int MASK = 67100672; // 2^26-2^13

//float exp=3;man=12
//const int MANTISSA_BITS = 12;
//const int MASK = 67106816; // 2^26-2^11

//float exp=3;man=13
//const int MANTISSA_BITS = 13;
//const int MASK = 67107840; // 2^26-2^10

//float exp=3;man=14
//const int MANTISSA_BITS = 14;
//const int MASK = 67108352; // 2^26-2^9

//float exp=3;man=16
//const int MANTISSA_BITS = 16;
//const int MASK = 67108736; // 2^26-2^7

class ExecStaticTable
{
public:
  using IdPortalTypeId = typename vtkm::cont::ArrayHandle<IdComponent>::ReadPortalType;
  using IdPortalTypeReal = typename vtkm::cont::ArrayHandle<Real>::ReadPortalType;

  ExecStaticTable(const IdPortalTypeId table_index_1,
                    const IdPortalTypeReal table_rij_1,
                    const IdPortalTypeReal table_drij_1,
                    const IdPortalTypeReal table_function_rij_1,
                    const IdPortalTypeReal table_dfunction_rij_1,
                    const IdPortalTypeId table_index_2,
                    const IdPortalTypeReal table_rij_2,
                    const IdPortalTypeReal table_drij_2,
                    const IdPortalTypeReal table_function_rij_2,
                    const IdPortalTypeReal table_dfunction_rij_2
                    )
    : _table_index_1(table_index_1)
    , _table_rij_1(table_rij_1)
    , _table_drij_1(table_drij_1)
    , _table_function_rij_1(table_function_rij_1)
    , _table_dfunction_rij_1(table_dfunction_rij_1)
    , _table_index_2(table_index_2)
    , _table_rij_2(table_rij_2)
    , _table_drij_2(table_drij_2)
    , _table_function_rij_2(table_function_rij_2)
    , _table_dfunction_rij_2(table_dfunction_rij_2)
  { 
  }


  VTKM_EXEC Real TableGnearValue(const Real& dis, IdComponent& index_rij) const
  {
    Real fraction, table;
    //vtkm::cont::Timer index_time;
    //index_time.Start();
    //IdComponent index_rij = Extract(dis);
    //Real time4index = index_time.GetElapsedTime();
    //printf("time for index is = %f\n", time4index);

    //vtkm::cont::Timer compute_time;
    //compute_time.Start();
    if (dis < 2.0)
    {
      index_rij = index_rij - _table_index_1.Get(0);
      table = _table_function_rij_1.Get(index_rij);
      //fraction = (dis - _table_rij_1.Get(index_rij)) * _table_drij_1.Get(index_rij);
      //table = _table_function_rij_1.Get(index_rij) + fraction * _table_dfunction_rij_1.Get(index_rij);
    }
    else
    {
      index_rij = index_rij - _table_index_2.Get(0);
      table = _table_function_rij_2.Get(index_rij);
      //fraction = (dis - _table_rij_2.Get(index_rij)) * _table_drij_2.Get(index_rij);
      //table = _table_function_rij_2.Get(index_rij) + fraction * _table_dfunction_rij_2.Get(index_rij);
    }
    //Real time4compute = compute_time.GetElapsedTime();
    //printf("time for compute is = %f\n", time4compute);

    return table;
  }

  VTKM_EXEC IdComponent Extract(const Real& t) const
  {
    IdComponent* int_p;
    IdComponent index;
    int_p = (IdComponent*)&t;
    index = (*int_p & MASK);
    index = index >> (23 - MANTISSA_BITS);
    return index;
  }

private:
  IdPortalTypeId _table_index_1;
  IdPortalTypeReal _table_rij_1;
  IdPortalTypeReal _table_drij_1;
  IdPortalTypeReal _table_function_rij_1;
  IdPortalTypeReal _table_dfunction_rij_1;
  IdPortalTypeId _table_index_2;
  IdPortalTypeReal _table_rij_2;
  IdPortalTypeReal _table_drij_2;
  IdPortalTypeReal _table_function_rij_2;
  IdPortalTypeReal _table_dfunction_rij_2;  
};
