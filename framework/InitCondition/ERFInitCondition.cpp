//#include "ERFInitCondition.h"
//#include <iostream>
//#include <vector>
//#include <fstream>
//#include <string>
//#include <algorithm>
//#include <vtkm/cont/ArrayCopy.h>
//#include "FieldName.h"
//
////RegisterObject(ERFInitCondition);
//ERFInitCondition::ERFInitCondition(const Configuration& cfg)
//  : InitCondition(cfg) 
//  , _use_erf(Get<bool>("use_erf"))
//  , _file_gnear_geq2("./res/Gnear_table_geq2_2.txt")
//  , _file_gnear_leq2("./res/Gnear_table_leq2_2.txt")
//{
//  _para.SetParameter(PARA_USE_ERF, _use_erf);
//}
//
//void ERFInitCondition::Execute()
//{
//  if (_use_erf == true)
//  {
//    try
//    {
//      ReadGnearGeq2File();
//      ReadGnearLeq2File();
//    }
//    catch (const std::exception& e)
//    {
//      _file_gnear_geq2.close();
//      _file_gnear_leq2.close();
//      std::cout << e.what() << std::endl;
//    }
//  }
//}
//
//void ERFInitCondition::ReadGnearGeq2File() 
//{
//    std::string line;
//    IdComponent index;
//    Real dis, dis_dif, gnear, gnear_der;
//    std::ofstream _file1("index_geq2.txt");
//    std::ofstream _file2("dis_geq2.txt");
//    std::ofstream _file3("dis_dif_geq2.txt");
//    std::ofstream _file4("gnear_geq2.txt");
//    std::ofstream _file5("gnear_der_geq2.txt");
//    while (getline(_file_gnear_geq2, line))
//    {
//      if (!line.empty() && line != "\r")
//      {
//        std::istringstream iss(line);
//        iss >> index >> dis >> dis_dif >> gnear >> gnear_der;
//        _file1 << index << ",";
//        _file2 << dis << ",";
//        _file3 << dis_dif << ",";
//        _file4 << gnear << ",";
//        _file5 << gnear_der << ",";
//        index_geq2.push_back(index);
//        dis_geq2.push_back(dis);
//        dis_dif_geq2.push_back(dis_dif);
//        gnear_geq2.push_back(gnear);
//        gnear_der_geq2.push_back(gnear_der);
//      }
//    }
//
//    auto erf_index_geq2 = _para.GetFieldAsArrayHandle<IdComponent>(field::erf_index_geq2);
//    vtkm::cont::ArrayCopy(vtkm::cont::make_ArrayHandle(index_geq2), erf_index_geq2);
//
//    auto erf_dis_geq2 = _para.GetFieldAsArrayHandle<Real>(field::erf_dis_geq2);
//    vtkm::cont::ArrayCopy(vtkm::cont::make_ArrayHandle(dis_geq2), erf_dis_geq2);
//
//    auto erf_dis_dif_geq2 = _para.GetFieldAsArrayHandle<Real>(field::erf_dis_dif_geq2);
//    vtkm::cont::ArrayCopy(vtkm::cont::make_ArrayHandle(dis_dif_geq2), erf_dis_dif_geq2);
//
//    auto erf_gnear_geq2 = _para.GetFieldAsArrayHandle<Real>(field::erf_gnear_geq2);
//    vtkm::cont::ArrayCopy(vtkm::cont::make_ArrayHandle(gnear_geq2), erf_gnear_geq2);
//
//    auto erf_gnear_der_geq2 = _para.GetFieldAsArrayHandle<Real>(field::erf_gnear_der_geq2);
//    vtkm::cont::ArrayCopy(vtkm::cont::make_ArrayHandle(gnear_der_geq2), erf_gnear_der_geq2);
//}
//
//void ERFInitCondition::ReadGnearLeq2File()
//{
//  std::string line;
//  IdComponent index;
//  Real dis, dis_dif, gnear, gnear_der;
//  std::ofstream _file1("index_leq2.txt");
//  std::ofstream _file2("dis_leq2.txt");
//  std::ofstream _file3("dis_dif_leq2.txt");
//  std::ofstream _file4("gnear_leq2.txt");
//  std::ofstream _file5("gnear_der_leq2.txt");
//  while (getline(_file_gnear_leq2, line))
//  {
//    if (!line.empty() && line != "\r")
//    {
//      std::istringstream iss(line);
//      iss >> index >> dis >> dis_dif >> gnear >> gnear_der;
//      _file1 << index << ",";
//      _file2 << dis << ",";
//      _file3 << dis_dif << ",";
//      _file4 << gnear << ",";
//      _file5 << gnear_der << ",";
//      index_leq2.push_back(index);
//      dis_leq2.push_back(dis);
//      dis_dif_leq2.push_back(dis_dif);
//      gnear_leq2.push_back(gnear);
//      gnear_der_leq2.push_back(gnear_der);
//    }
//  }
//
//  auto erf_index_leq2 = _para.GetFieldAsArrayHandle<IdComponent>(field::erf_index_leq2);
//  vtkm::cont::ArrayCopy(vtkm::cont::make_ArrayHandle(index_leq2), erf_index_leq2);
//
//  auto erf_dis_leq2 = _para.GetFieldAsArrayHandle<Real>(field::erf_dis_leq2);
//  vtkm::cont::ArrayCopy(vtkm::cont::make_ArrayHandle(dis_leq2), erf_dis_leq2);
//
//  auto erf_dis_dif_leq2 = _para.GetFieldAsArrayHandle<Real>(field::erf_dis_dif_leq2);
//  vtkm::cont::ArrayCopy(vtkm::cont::make_ArrayHandle(dis_dif_leq2), erf_dis_dif_leq2);
//
//  auto erf_gnear_leq2 = _para.GetFieldAsArrayHandle<Real>(field::erf_gnear_leq2);
//  vtkm::cont::ArrayCopy(vtkm::cont::make_ArrayHandle(gnear_leq2), erf_gnear_leq2);
//
//  auto erf_gnear_der_leq2 = _para.GetFieldAsArrayHandle<Real>(field::erf_gnear_der_leq2);
//  vtkm::cont::ArrayCopy(vtkm::cont::make_ArrayHandle(gnear_der_leq2), erf_gnear_der_leq2);
//}
