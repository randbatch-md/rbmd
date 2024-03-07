#include "TestMadelung.h"
#include "MDApplication.h"
#include "System.h"
#include <gtest/gtest.h>
#include "FieldName.h"
#include <fmt/format.h>
#include <fmt/args.h>
#include <fstream>

  TestMadelung::TestMadelung()
  : _madelung_result(0.0)
{
}

  void TestMadelung::MadelungFormater(const std::string& type, const int dims)
  {
    std::string madelung_parameter = R"(
    # 配置参数
    [Parallel]
      type = {_type}
      num_threads = 8
    []
    
    [InitCondition]
       type = MadelungInitCondition
       dims = '{_dims} {_dims} {_dims}'
       x_range = '0.0 {_dims}'
       x_range = '0.0 {_dims}'
       y_range = '0.0 {_dims}'
       z_range = '0.0 {_dims}'
    []
    
    [System]
      type = MadelungSystem
    []

    [Executioner]
       dt = 1e-05
       num_steps = 1
    []
    )";

    fmt::dynamic_format_arg_store<fmt::format_context> args;
    args.push_back(fmt::arg("_type", type));
    args.push_back(fmt::arg("_dims", dims));

    std::ofstream outFile;
    outFile.open("madelung_test.i"); 
    outFile << fmt::vformat(madelung_parameter, args);
  }

  void TestMadelung::ComputeMadelung()
  {
    int argc = 3;
    char* argv[3] = {
      " ",
      "-i",
      "madelung_test.i",
    };
  
    std::shared_ptr<Application> app = std::make_shared<MDApplication>(argc, argv);
    app->Run();
  
    _madelung_result = (app->GetSystem())->GetParameter<Real>(gtest::madelung);
  }

TEST_F(TestMadelung, MadeLung)
{
  //madelung test
  std::cout << "当前测试文件为：TestMadelung" << std::endl;
  MadelungFormater("serial");
  ComputeMadelung();
  EXPECT_NEAR(_madelung_result, -1.759, 0.001);
  
  MadelungFormater("cuda");
  ComputeMadelung();
  EXPECT_NEAR(_madelung_result, -1.759, 0.001);
}