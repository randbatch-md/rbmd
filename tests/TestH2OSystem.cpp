#include "TestH2OSystem.h"
#include "MDApplication.h"
#include "System.h"
#include <gtest/gtest.h>
#include "FieldName.h"
#include <fmt/format.h>
#include <fmt/args.h>
#include <fstream>

TestH2OSystem::TestH2OSystem()
  : _kbt(0.0)
  , _ave_kbt(0.0)
  , _ave_pe(0.0)
  , _ave_ke(0.0)
{
}

void TestH2OSystem::H2OSystemFormater(const std::string& type,
                                      const int force_type,
                                      const double dt,
                                      const int num_steps)
{
    std::string H20System_parameter = R"(
    # 配置参数
    [Parallel]
      type = {_type}
      num_threads = 8
    []

    [Locator]
      cut_off = 10
    []

    [InitCondition]
       type = ModelFileInitCondition
       file = 'H2O_v.data'
       velocity_type = false
    []

    [System]
      type = H2OSystem
      alpha = 0.1
      rbeP = 100
      Kmax = 4
      kbT = 298
      force_type = {_force_type} #0:RBE, 1:Ewald
      temp_con_type = 3 #0:NOSE_HOOVER, 1:LANGEVIN, 2:TEMP_RESCALE, 3:BERENDSEN
      unit = REAL   #0:REAL, 1:LJ
      use_shake = false
      flag_steps = 0
    []

    [Executioner]
       dt = {_dt}
       num_steps = {_num_steps}
    []

    [Outputs]
    [./TableOutput]
      type = ThermoOutput
      interval = 1
      binary = false
      out_initial = false
      output_screen=false
      output_file=false
      compute=true
    [../]

    [./FileOutput]
      type = RDFOutput
      interval = 10
      binary = false
      out_initial = false
      output_file=false
      compute=false
      radius = 10 
      min_radius = 0
      dr = 0.01
      statistics_rdf_steps = 1000
      center_type = 0;
      target_type = 1;
    [../]

    [./MSDOutput]
      type = MSDOutput
      interval = 10
      binary = false
      out_initial = false
      output_file=false
      compute=false
      start_step = 500
      end_step = 600
    [../]

    [./Progress]
      type = Progress
    [../]
    []
    )";

    fmt::dynamic_format_arg_store<fmt::format_context> args;
    args.push_back(fmt::arg("_type", type));
    args.push_back(fmt::arg("_force_type", force_type));
    args.push_back(fmt::arg("_dt", dt));
    args.push_back(fmt::arg("_num_steps", num_steps));

    std::ofstream outFile;
    outFile.open("LJFluid_real.i");
    outFile << fmt::vformat(H20System_parameter, args);
}

void TestH2OSystem::ComputeH2OSystem()
{
    int argc = 3;
    char* argv[3] = {
    " ",
    "-i",
    "LJFluid_real.i",
    };

    std::shared_ptr<Application> app = std::make_shared<MDApplication>(argc, argv);
    app->Run();

    _kbt = (app->GetParameter())->GetParameter<Real>(PARA_TEMPT);
    _ave_kbt = (app->GetParameter())->GetParameter<Real>(gtest::ave_temp_t);
    _ave_pe = (app->GetParameter())->GetParameter<Real>(gtest::ave_potential_energy);
    _ave_ke = (app->GetParameter())->GetParameter<Real>(gtest::ave_kin_energy);
}

TEST_F(TestH2OSystem, thermo_out)
{
  // H2OSystem Ewald Test
  std::cout << "当前测试文件为：TestH2OSystem Ewald" << std::endl;
  H2OSystemFormater("serial", 1, 1, 50);
  ComputeH2OSystem();
  EXPECT_NEAR(_kbt, 328.073, 0.001);     //EXPECT_EQ(_kbt, 328.073);
  EXPECT_NEAR(_ave_kbt, 291.644, 0.001); //EXPECT_EQ(_ave_kbt, 291.644);
  EXPECT_NEAR(_ave_pe, -2.30006, 0.001); //EXPECT_EQ(_ave_pe, -2.30006);
  EXPECT_NEAR(_ave_ke, 0.869335, 0.001); //EXPECT_EQ(_ave_ke, 0.869335);
  
  H2OSystemFormater("cuda", 1, 1, 50);
  ComputeH2OSystem();
  EXPECT_NEAR(_kbt, 328.075, 0.001);     //EXPECT_EQ(_kbt, 328.075);
  EXPECT_NEAR(_ave_kbt, 291.644, 0.001); //EXPECT_EQ(_ave_kbt, 291.644);
  EXPECT_NEAR(_ave_pe, -2.30006, 0.001); //EXPECT_EQ(_ave_pe, -2.30008);
  EXPECT_NEAR(_ave_ke, 0.869335, 0.001); //EXPECT_EQ(_ave_ke, 0.869335);

  // H2OSystem RBE test
  std::cout << "当前测试文件为：TestH2OSystem RBE" << std::endl;
  H2OSystemFormater("serial", 0, 1, 50);
  ComputeH2OSystem();
  EXPECT_NEAR(_kbt, 338.694855, 0.001);     //EXPECT_EQ(_kbt, 328.073);
  EXPECT_NEAR(_ave_kbt, 294.75485, 0.001); //EXPECT_EQ(_ave_kbt, 291.644);
  EXPECT_NEAR(_ave_pe, -2.2284889, 0.001);  //EXPECT_EQ(_ave_pe, -2.30006);
  EXPECT_NEAR(_ave_ke, 0.878608, 0.001);   //EXPECT_EQ(_ave_ke, 0.869335);
  
  H2OSystemFormater("cuda", 0, 1, 50);
  ComputeH2OSystem();
  EXPECT_NEAR(_kbt, 338.694885, 0.001);      //EXPECT_EQ(_kbt, 328.075);
  EXPECT_NEAR(_ave_kbt, 294.75476, 0.001); //EXPECT_EQ(_ave_kbt, 291.644);
  EXPECT_NEAR(_ave_pe, -2.228508, 0.001);    //EXPECT_EQ(_ave_pe, -2.30008);
  EXPECT_NEAR(_ave_ke, 0.8786079, 0.001);  //EXPECT_EQ(_ave_ke, 0.869335);
}