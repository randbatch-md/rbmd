﻿//==================================================================================
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

#include "FieldName.h"
#include "MDApplication.h"
#include "System.h"
#include "TestLennardJonesSystem.h"
#include <fmt/args.h>
#include <fmt/format.h>
#include <fstream>
#include <gtest/gtest.h>

TestLennardJonesSystem::TestLennardJonesSystem()
  : _kbt(0.0)
{
}

void TestLennardJonesSystem::LennardJonesSystemFormater(const std::string& type,
                                               const int dims,
                                               const double dt,
                                               const int num_steps)
{
    std::string LJSystem_parameter = R"(
    # 配置参数
    [Parallel]
      type = {_type}
      num_threads = 8
    []

    [Locator]
      cut_off = 5
    []

    [InitCondition]
       type = LJInitCondition
       max_steps = 0
       dt = 1e-05
       dims = '{_dims} {_dims} {_dims}'
       x_range = '0.0 {_dims}'
       x_range = '0.0 {_dims}'
       y_range = '0.0 {_dims}'
       z_range = '0.0 {_dims}'
       velocity_type = false
    []

    [System]
      type = LennardJonesSystem
      alpha = 0.1
      rbeP = 100
      Kmax = 2
      kbT = 1.0
      force_type = 1 #0:RBE, 1:Ewald
      temp_con_type = 3 #0:NOSE_HOOVER, 1:LANGEVIN, 2:TEMP_RESCALE, 3:BERENDSEN
      unit = LJ   #0:REAL, 1:LJ
    []

    [Executioner]
       type = Transient
       dt = {_dt}
       num_steps = {_num_steps}
    []

    [Outputs]
    [./TableOutput]
      type = AtomsTableOutput
      interval = 1
      binary = false
      out_initial = false
      output_screen=false
      output_file=false
      compute=true
    [../]

    [./FileOutput]
      type = AtomsFileOutput
      interval = 10
      binary = false
      out_initial = false
      output_file = false
      compute = false
      radius = 5
      min_radius = 0.0
      dr = 0.01
      statistics_rdf_steps = 1000
      center_type = 1;
      target_type = 1;
    [../]

    [./MSDOutput]
      type = MSDOutput
      interval = 1
      binary = false
      out_initial = false
      output_file=false
      compute=false
      start_step = 500
      end_step = 600
    [../]
    []
    )";

    fmt::dynamic_format_arg_store<fmt::format_context> args;
    args.push_back(fmt::arg("_type", type));
    args.push_back(fmt::arg("_dims", dims));
    args.push_back(fmt::arg("_dt", dt));
    args.push_back(fmt::arg("_num_steps", num_steps));

    std::ofstream outFile;
    outFile.open("LJFluid.i");
    outFile << fmt::vformat(LJSystem_parameter, args);
}

void TestLennardJonesSystem::ComputeLennardJonesSystem()
{
  int argc = 3;
  char* argv[3] = {
    " ",
    "-i",
    "LJFluid.i",
  };

  std::shared_ptr<Application> app = std::make_shared<MDApplication>(argc, argv);
  app->Run();
  _kbt = (app->GetParameter())->GetParameter<Real>(PARA_TEMPT);
}

TEST_F(TestLennardJonesSystem, thermo_out)
{
  //LJSystem test
  std::cout << "当前测试文件为：TestLennardJonesSystem" << std::endl;
  LennardJonesSystemFormater("serial", 10, 0.002, 50);
  ComputeLennardJonesSystem();
  EXPECT_NEAR(_kbt, 0.387303, 0.001); //EXPECT_EQ(_kbt, 1.46822);

  LennardJonesSystemFormater("cuda", 10, 0.002, 50); 
  ComputeLennardJonesSystem();
  EXPECT_NEAR(_kbt, 0.3873029, 0.001); //EXPECT_EQ(_kbt, 0.82059);

  LennardJonesSystemFormater("cuda", 10, 0.005, 50); 
  ComputeLennardJonesSystem();
  EXPECT_NEAR(_kbt, 0.3873047, 0.001); //EXPECT_EQ(_kbt, 1.37518);
}