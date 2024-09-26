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
//  Contact Email : [your - email@example.com]
//==================================================================================

#pragma once
#include <vtkm/VectorAnalysis.h>
#include <vtkm/cont/ArrayHandle.h>
#include "Types.h"
#include <vtkm/cont/ConvertNumComponentsToOffsets.h>
#include <vtkm/cont/ArrayHandleGroupVecVariable.h>
#include <vtkm/cont/ArrayHandleGroupVec.h>

class ExecTopology
{
public:
  using IdPortalTypeId = typename vtkm::cont::ArrayHandle<vtkm::Id>::ReadPortalType;
  using IdPortalTypeReal = typename vtkm::cont::ArrayHandle<Real>::ReadPortalType;
  using GroupIdPortalTypeId = typename vtkm::cont::ArrayHandleGroupVecVariable<vtkm::cont::ArrayHandle<vtkm::Id>,
                                                                               vtkm::cont::ArrayHandle<vtkm::Id>>::ReadPortalType;

  ExecTopology(const GroupIdPortalTypeId groupVecArray,
               const IdPortalTypeId molecularId,
               const IdPortalTypeId atoms_type,
               const IdPortalTypeReal charge,
               const IdPortalTypeReal epsilon,
               const IdPortalTypeReal sigma)
               //,
               //const GroupIdPortalTypeId groupVecNeighbour) 
     : _group_vec_portal(groupVecArray)
     , _molecular_id_portal(molecularId)
     , _atoms_type_portal(atoms_type)
     , _charge_portal(charge)
     , _eps_portal(epsilon)
     , _sig_portal(sigma)
     //, _group_vec_neighbour_portal(groupVecNeighbour)
  {
  }

  template<typename T >
  VTKM_EXEC vtkm::VecFromPortal<T> GetGroupVecArray( vtkm::Id&  atoms_id) const
  { 
		return this->_group_vec_portal.Get(atoms_id);
  }

  VTKM_EXEC vtkm::Id GetMolecularId(const Id& atoms_id) const
  {
        return this->_molecular_id_portal.Get(atoms_id);
  }

  VTKM_EXEC vtkm::Id GetAtomsType(const Id& atoms_id) const
  {
        return this->_atoms_type_portal.Get(atoms_id);
  }
  VTKM_EXEC Real GetCharge(const Id& atoms_id) const { return this->_charge_portal.Get(atoms_id); }

  VTKM_EXEC Real GetEpsilon(const Id& atoms_type) const{return this->_eps_portal.Get(atoms_type);}

  VTKM_EXEC Real GetSigma(const Id& atoms_type) const { return this->_sig_portal.Get(atoms_type);}

  // template<typename T >
  // VTKM_EXEC vtkm::VecFromPortal<T> GetGroupVecNeighbourArray( vtkm::Id&  atoms_id) const
  // { 
	// 	return this->_group_vec_neighbour_portal.Get(atoms_id);
  // }

private:

  GroupIdPortalTypeId _group_vec_portal;
  IdPortalTypeId _molecular_id_portal;
  IdPortalTypeId _atoms_type_portal;
  IdPortalTypeReal _charge_portal;
  IdPortalTypeReal _eps_portal;
  IdPortalTypeReal _sig_portal;

  //GroupIdPortalTypeId _group_vec_neighbour_portal;
};