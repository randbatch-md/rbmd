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

#include <fstream>
#include "Logging.h"
#include "json/value.h"
#include "json/reader.h"

class JsonParser
{
public:
  JsonParser(const std::string& file)
  {
    if (!IsJsonFile(file))
    {
      console::Error(file, "is not json file!");
      return;
    }

    ParseJsonFile(file);
  }

  ~JsonParser() = default;

public:
  auto& GetJsonNode(const std::string& key)
  {
    if (!_json_node[key].isObject())
    {
      console::Warning("Can not find key: ", key);
    }

    return _json_node[key];
  }

  bool HasNode(const std::string& key) 
  {
      return _json_node.isMember(key);
  }

  std::string GetFileStr() 
  {
    return _buffer.str();
  }

private:

  bool IsJsonFile(const std::string& file) 
  { 
	  auto length = file.length();
      return (length >= 5 && file.substr(length - 5) == ".json");
  }

  void ParseJsonFile(const std::string& file) 
  {
    try
    {
      std::ifstream filestream(file);
      if (!filestream.is_open())
      {
        console::Error("failed to open file: ", file);
        return;
      }
      _buffer << filestream.rdbuf();
      filestream.clear();  // 清除 EOF 标志
      filestream.seekg(0); // 移动到文件开始位置
      Json::CharReaderBuilder readerBuilder;
      Json::parseFromStream(readerBuilder, filestream, &_json_node, nullptr);

      //关闭文件
      filestream.close();

      //此处可以删除注释，用于以后做支持注释的扩展
    }
    catch (const std::exception&)
    {
      console::Error("Error parsing JSON");
    }
  }

private:
  Json::Value _json_node;
  std::ostringstream _buffer;
};