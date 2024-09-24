#pragma once
#include "TableBase.h"
#include <fstream>
#include <iostream>
class FormatTable : public TableBase
{
public:
  FormatTable();
  virtual ~FormatTable();

  bool AddHeader(const std::string& name) override;
  bool AddData(const std::string& name, const std::string& data) override;
  void Print() override;
  bool IsEmpty();
  void LogFile();

protected:
  void PrintHeader();
  void PrintRowData();

protected:
  const unsigned int _header_spacing; // use to print header
  const unsigned int _column_width;
  long _row_count;
  std::ofstream _log_file;
};