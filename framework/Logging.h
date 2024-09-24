#pragma once

#include <fmt/color.h>
#include <sstream>
#include <vtkm/cont/Logging.h>

inline void ToStringStream(std::ostringstream& oss) {}

template<typename T, typename... Args>
void ToStringStream(std::ostringstream& oss, T&& value, Args&&... args)
{
  oss << value << " ";
  ToStringStream(oss, std::forward<Args>(args)...);
}

namespace console
{
template<typename... Args>
[[noreturn]] void Error(Args&&... args)
{
  std::ostringstream oss;
  ToStringStream(oss, std::forward<Args>(args)...);
  fmt::print(fmt::fg(fmt::color::red), "[{0}] {1}\n", "Error", oss.str());
  std::exit(0);
}
template<typename... Args>
void Message(Args&&... args)
{
  std::ostringstream oss;
  ToStringStream(oss, std::forward<Args>(args)...);
  fmt::print(fmt::fg(fmt::color::white), "[{0}] {1}\n", "Msg", oss.str());
}
template<typename... Args>
void Info(Args&&... args)
{
  std::ostringstream oss;
  ToStringStream(oss, std::forward<Args>(args)...);
  fmt::print(fmt::fg(fmt::color::white), "[{0}] {1}\n", "Info", oss.str());
}
template<typename... Args>
void Success(Args&&... args)
{
  std::ostringstream oss;
  ToStringStream(oss, std::forward<Args>(args)...);
  fmt::print(fmt::fg(fmt::color::green), "[{0}] {1}\n", "Success", oss.str());
}

template<typename... Args>
void Warning(Args&&... args)
{
  std::ostringstream oss;
  ToStringStream(oss, std::forward<Args>(args)...);
  fmt::print(fmt::fg(fmt::color::yellow), "[{0}] {1}\n", "Warning", oss.str());
}
}