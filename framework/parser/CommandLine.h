//==================================================================================
//  RBMD 2.1.0 is developed for random batch molecular dynamics calculation.
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

#include "Logging.h"
#include "cxxopts.hpp"
#include <memory>
#include <string>
#include <vector>

class CommandLine
{
public:
  CommandLine(int argc, char* argv[])
  {
    for (size_t i = 0; i < argc; i++)
    {
      _args.push_back(argv[i]);
    }
    auto opts = std::make_shared<cxxopts::Options>("opts");
    opts->add_options()("j", "json file", cxxopts::value<std::string>());
    _co = opts->parse(argc, argv);

    std::string ifile;
    if (_co.count("j"))
    {
      ifile = _co["j"].as<std::string>();
    }

    if (this->Argc() == 1 || _co.count("j") == 0 || ifile.empty())
    {
      this->PrintUsage();
      std::exit(0);
    }
  }

  ~CommandLine() {}

  void PrintUsage()
  {
    auto usage = fmt::format("Usage: {} -j <input.json>", _args[0]);
    std::cout << usage;
  }

  //GetPot& Get() { return *_cl; }
  int Argc() { return _args.size(); }

  template<typename T>
  T GetValue(const std::string& name)
  {
    try
    {
      return _co[name].as<T>();
    }
    catch (const std::exception& e)
    {
      console::Error(e.what());
    }
  }

private:
  cxxopts::ParseResult _co;
  std::vector<std::string> _args;
};
