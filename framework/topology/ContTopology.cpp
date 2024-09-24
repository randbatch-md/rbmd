#include "topology/ContTopology.h"
#include <vtkm/cont/ArrayHandle.h>
#include  "ContTopology.h"

VTKM_CONT ExecTopology ContTopology::PrepareForExecution(vtkm::cont::DeviceAdapterId device,
                                                         vtkm::cont::Token& token) const
{
  return ExecTopology(this->_group_vec_array.PrepareForInput(device, token), 
                      this->_molecular_id.PrepareForInput(device, token),
                      this->_atoms_type.PrepareForInput(device, token),
                      this->_charge.PrepareForInput(device, token),
                      this->_epsilon.PrepareForInput(device, token),
                      this->_sigma.PrepareForInput(device, token) );
                      // this->_group_vec_neighbour.PrepareForInput(device, token) 
                      
}