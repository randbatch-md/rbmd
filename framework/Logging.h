//==================================================================================
//  RBMD 1.0 is developed for random batch molecular dynamics calculation.
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
//  Contact Email : [your - email@example.com]
//==================================================================================

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