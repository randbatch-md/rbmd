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
