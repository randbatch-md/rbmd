#include "FieldName.h"
#include "MDApplication.h"
#include "System.h"
#include "TestNaClSystem.h"
#include <fmt/args.h>
#include <fmt/format.h>
#include <fstream>
#include <gtest/gtest.h>

TestNaClSystem::TestNaClSystem()
  : _kbt(0.0)
  , _ave_kbt(0.0)
  , _ave_pe(0.0)
  , _ave_ke(0.0)
{
}

void TestNaClSystem::NaClSystemFormater(const std::string& type,
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
       file = 'NaClsystem.data'
       velocity_type = false
    []

    [System]
      type = NaClSystem
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
      type = MoleculesTableOutput
      interval = 1
      binary = false
      out_initial = false
      output_screen=false
      output_file=false
      compute=true
    [../]

    [./FileOutput]
      type = MoleculesFileOutput
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

void TestNaClSystem::ComputeNaClSystem()
{
  int argc = 3;
  char* argv[3] = {
    " ",
    "-i",
    "LJFluid_real.i",
  };

  std::shared_ptr<Application> app = std::make_shared<MDApplication>(argc, argv);
  app->Run();

  _kbt = (app->GetSystem())->GetParameter<Real>(PARA_TEMPT);
  _ave_kbt = (app->GetSystem())->GetParameter<Real>(gtest::ave_temp_t);
  _ave_pe = (app->GetSystem())->GetParameter<Real>(gtest::ave_potential_energy);
  _ave_ke = (app->GetSystem())->GetParameter<Real>(gtest::ave_kin_energy);
}

TEST_F(TestNaClSystem, thermo_out)
{
  // NaClSystem Ewald Test
  std::cout << "当前测试文件为：TestNaClSystem Ewald" << std::endl;
  NaClSystemFormater("serial", 1, 0.2, 10);
  ComputeNaClSystem();
  EXPECT_NEAR(_kbt, 333.34304809570312, 0.001); //EXPECT_EQ(_kbt, 328.073);
  EXPECT_NEAR(_ave_kbt, 226.84980773925781, 0.001); //EXPECT_EQ(_ave_kbt, 291.644);
  EXPECT_NEAR(_ave_pe, -0.14420869946479797, 0.001); //EXPECT_EQ(_ave_pe, -2.30006);
  EXPECT_NEAR(_ave_ke, 0.67619627714157104, 0.001);  //EXPECT_EQ(_ave_ke, 0.869335);
  
  NaClSystemFormater("cuda", 1, 0.2, 10);
  ComputeNaClSystem();
  EXPECT_NEAR(_kbt, 333.34298706054688, 0.001); //EXPECT_EQ(_kbt, 328.075);
  EXPECT_NEAR(_ave_kbt, 226.84988403320312, 0.001); //EXPECT_EQ(_ave_kbt, 291.644);
  EXPECT_NEAR(_ave_pe, -0.14516164362430573, 0.001); //EXPECT_EQ(_ave_pe, -2.30008);
  EXPECT_NEAR(_ave_ke, 0.67619645595550537, 0.001);  //EXPECT_EQ(_ave_ke, 0.869335);

  // NaClSystem RBE test(RBE算法会发散)
  std::cout << "当前测试文件为：TestNaClSystem RBE" << std::endl;
  NaClSystemFormater("serial", 0, 0.2, 10);
  ComputeNaClSystem();
  EXPECT_NEAR(_kbt, 312.53936767578125, 0.001); //EXPECT_EQ(_kbt, 328.073);
  EXPECT_NEAR(_ave_kbt, 257.19378662109375, 0.001); //EXPECT_EQ(_ave_kbt, 291.644);
  EXPECT_NEAR(_ave_pe, 0.05721588060259819, 0.001); //EXPECT_EQ(_ave_pe, -2.30006);
  EXPECT_NEAR(_ave_ke, 0.76664584875106812, 0.001); //EXPECT_EQ(_ave_ke, 0.869335);

  NaClSystemFormater("cuda", 0, 0.2, 10);
  ComputeNaClSystem();
  EXPECT_NEAR(_kbt, 312.53903198242188, 0.001); //EXPECT_EQ(_kbt, 328.075);
  EXPECT_NEAR(_ave_kbt, 257.19363403320312, 0.001); //EXPECT_EQ(_ave_kbt, 291.644);
  EXPECT_NEAR(_ave_pe, 0.056220255792140961, 0.001); //EXPECT_EQ(_ave_pe, -2.30008);
  EXPECT_NEAR(_ave_ke, 0.76664537191390991, 0.001);  //EXPECT_EQ(_ave_ke, 0.869335);
}