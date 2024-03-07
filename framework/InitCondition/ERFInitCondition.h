//#pragma once
//#include "InitCondition.h"
//#include <vector>
//
//class ERFInitCondition : public InitCondition
//{
//public:
//  ERFInitCondition(const Configuration& cfg);
//  virtual ~ERFInitCondition() = default;
//
//  void Execute() override;
//
//protected:
//  void UpdateField() override{};
//
//private:
//  void ReadGnearGeq2File();
//  void ReadGnearLeq2File();
// 
//private:
//  bool _use_erf;
//  std::ifstream _file_gnear_geq2;
//  std::ifstream _file_gnear_leq2;
//  std::vector<IdComponent> index_geq2;
//  std::vector<IdComponent> index_leq2;
//  std::vector<Real> dis_geq2;
//  std::vector<Real> dis_leq2;
//  std::vector<Real> dis_dif_geq2;
//  std::vector<Real> dis_dif_leq2;
//  std::vector<Real> gnear_geq2;
//  std::vector<Real> gnear_leq2;
//  std::vector<Real> gnear_der_geq2;
//  std::vector<Real> gnear_der_leq2;
//};
