#pragma once

#include "Output.h"
#include "model/FormatTable.h"
#include "System.h"
class ConsoleOutput : public Output
{
public:
  ConsoleOutput(const Configuration& cfg);
  virtual ~ConsoleOutput(){};

  void Init() override;
  void Execute() override;

protected:
  bool AddHeader(const std::string& name);

  template<typename DataType>
  bool AddData(const std::string& headerName, const DataType& data)
  {
    return _format_table->AddData(headerName, std::to_string(data));
  }

private:
  void EssientialOutput();

protected:
  std::unique_ptr<FormatTable> _format_table;
  std::string _name;
  Timer& _timer;
  Real _time;
  Real _cumulative_time;
  bool _output_screen;
};
