#pragma once

#include <unordered_map>
#include <string>
class TableBase
{
public:
  TableBase() = default;
  virtual ~TableBase() = default;

  virtual bool AddHeader(const std::string& name) = 0;
  virtual bool AddData(const std::string& name, const std::string& data) = 0;
  virtual void Print() = 0;

protected:
  std::unordered_map<std::string, std::string> _row_data;
};