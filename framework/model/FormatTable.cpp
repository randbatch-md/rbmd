#include "FormatTable.h"
#include "Logging.h"
#include <fmt/printf.h>
#include <fmt/color.h>
#include <fmt/args.h>

FormatTable::FormatTable()
  : _header_spacing(16)
  , _column_width(16)
  , _row_count(0)
  , _log_file("LogFile.csv")
{
}

FormatTable::~FormatTable()
{
  _log_file.close();
}

bool FormatTable::AddHeader(const std::string& name)
{
  auto header = _row_data.find(name);
  if (_row_data.end() != header)
  {
    return false;
  }

  _row_data.insert(std::make_pair(name, ""));

  return true;
}

bool FormatTable::AddData(const std::string& name, const std::string& data)
{
  auto column = _row_data.find(name);
  if (_row_data.end() == column)
  {
    return false;
  }

  _row_data[name] = data;

  return true;
}

void FormatTable::Print()
{
  if (_row_count % _header_spacing == 0)
  {
    PrintHeader();
  }

  PrintRowData();
}

void FormatTable::PrintHeader() 
{
  auto begin = _row_data.begin();
  auto end = _row_data.end();
  auto color = fmt::fg(fmt::color::yellow);
  for (auto iter = begin; iter != end; ++iter)
  {
    if (iter == begin)
      fmt::print(color, "┌{0:-^{1}}", "", _column_width);
    else if (std::next(iter) == end)
      fmt::print(color, "{0:-^{1}}┐\n", "", _column_width);
    else
      fmt::print(color, "{0:-^{1}}", "", _column_width);
  }

  for (auto iter = begin; iter != end; ++iter)
  {
    if (iter == begin)
      fmt::print(color, "│{0: ^{1}}", iter->first, _column_width);
    else if (std::next(iter) == end)
      fmt::print(color, "{0: ^{1}}│\n", iter->first, _column_width);
    else
      fmt::print(color, "{0: ^{1}}", iter->first, _column_width);
  }

  for (auto iter = begin; iter != end; ++iter)
  {
    if (iter == begin)
      fmt::print(color, "└{0:-^{1}}", "", _column_width);
    else if (std::next(iter) == end)
      fmt::print(color, "{0:-^{1}}┘\n", "", _column_width);
    else
      fmt::print(color, "{0:-^{1}}", "", _column_width);
  }
}

void FormatTable::PrintRowData() 
{
  auto begin = _row_data.begin();
  auto end = _row_data.end();
  auto color = fmt::fg(fmt::color::white);

  for (auto iter = begin; iter != end; ++iter)
  {
    if (iter == begin)
      fmt::print(color, "│{0: ^{1}}", iter->second, _column_width);
    else if (std::next(iter) == end)
      fmt::print(color, "{0: ^{1}}│\n", iter->second, _column_width);
    else
      fmt::print(color, "{0: ^{1}}", iter->second, _column_width);
  }

  for (auto iter = begin; iter != end; ++iter)
  {
    if (iter == begin)
      fmt::print(color, "└{0:-^{1}}", "", _column_width);
    else if (std::next(iter) == end)
      fmt::print(color, "{0:-^{1}}┘\n", "", _column_width);
    else
      fmt::print(color, "{0:-^{1}}", "", _column_width);
  }

  ++_row_count;
}

bool FormatTable::IsEmpty()
{
  return _row_data.empty();
}

void FormatTable::LogFile() 
{
  try
  {
    if ((_row_count - 1) % _header_spacing == 0)
    {
      for (const auto& value : _row_data)
      {
        _log_file << value.first << "   ";
      }
      _log_file << std::endl;
    }
    for (const auto& value : _row_data)
    {
      _log_file << value.second << "   ";
    }
    _log_file << std::endl;
  }
  catch (const std::exception& e)
  {
    _log_file.close();
    console::Error(e.what());
  }
}
