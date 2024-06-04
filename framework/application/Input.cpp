#include "Input.h"

Input::Input(int argc, char** argv)
  : Application(argc,argv)
{
  //_cfg = std::make_shared<Configuration>();
}

//void Input::ExecuteCommand(const Object& obj) 
//{
//  obj.Get<vtkm::IdComponent>("abd");
//}
void Input::ExecuteCommandom() {}

void Input::InitConfigurationInBuild() 
{
  auto _dims = _obj->GetVectorValue<int>("dims");
  auto _x_range= _obj->GetVectorValue<Real>("x_range");
  auto _y_range= _obj->GetVectorValue<Real>("y_range");
  auto _z_range= _obj->GetVectorValue<Real>("z_range");

}

void Input::InitConfigurationReadData() {}
