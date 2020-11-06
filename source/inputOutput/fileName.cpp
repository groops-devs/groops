/***********************************************/
/**
* @file fileName.cpp
*
* @brief File names
*
* @author Torsten Mayer-Guerr
* @date 2008-07-28
*
*/
/***********************************************/

#include "base/importStd.h"
#include "base/string.h"
#include "parser/stringParser.h"
#include "fileName.h"

/*************************************************/

FileName FileName::append(const FileName& fileName) const
{
  if(empty())
    return fileName;
  std::string separator;
  if((name_.back() != '/') && (name_.back() != '\\'))
#ifdef _WIN32
    separator = '\\'; // = std::filesystem::path::preferred_separator;
#else
    separator = '/';
#endif
  return FileName(name_+separator+fileName.str());
}

/*************************************************/

FileName FileName::fullExtension() const
{
  auto pos = name_.rfind('.');
  if((pos == std::string::npos) || (pos+1 == name_.size()))
    return FileName();
  std::string ext = String::upperCase(name_.substr(pos+1));
  if((pos > 0) && (ext == "GZ" || ext == "Z"))
  {
    auto posNew = name_.rfind('.', pos-1);
    if(posNew != std::string::npos)
      pos = posNew;
  }
  return FileName(name_.substr(pos+1));
}

/*************************************************/

FileName FileName::typeExtension() const
{
  std::string ext = fullExtension();
  return FileName(ext.substr(0, ext.find('.')));
}

/*************************************************/

FileName FileName::packExtension() const
{
  std::string ext = fullExtension();
  auto pos = ext.find('.');
  if(pos == std::string::npos)
    return FileName();
  return FileName(ext.substr(pos+1));
}

/*************************************************/

FileName FileName::stripFullExtension() const
{
  auto pos = name_.rfind("."+fullExtension().str());
  if(pos == std::string::npos)
    return *this;
  return FileName(name_.substr(0, pos));
}

/*************************************************/

FileName FileName::replaceFullExtension(const std::string &text) const
{
  return FileName(stripFullExtension().str() + ((text.empty() || String::startsWith(text, ".")) ? "" : ".") + text);
}

/*************************************************/

FileName FileName::directory() const
{
  auto pos = name_.find_last_of("/\\");
  if(pos == std::string::npos)
    return FileName();
  return FileName(name_.substr(0, pos+1));
}

/*************************************************/

FileName FileName::stripDirectory() const
{
  std::string::size_type pos = name_.find_last_of("/\\");
  if(pos == std::string::npos)
    return *this;
  return FileName(name_.substr(pos+1));
}

/*************************************************/

FileName FileName::appendBaseName(const std::string &text) const
{
  auto pos = name_.rfind("."+fullExtension().str());
  if(pos == std::string::npos)
    return FileName(name_+text);
  std::string tmp = name_;
  return FileName(tmp.insert(pos, text));
}

/*************************************************/

FileName FileName::operator()(const VariableList &varList) const
{
  try
  {
    return FileName(StringParser::parse(name_, varList));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/*************************************************/
