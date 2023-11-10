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

FileName::FileName() : resolved(TRUE) {}

/*************************************************/

FileName::FileName(const std::string &name)
  : nameUnparsed(name), nameParsed(name), resolved(TRUE) {}

/*************************************************/

FileName::FileName(const std::string &name, const VariableList &varList)
  : nameUnparsed(name), resolved(FALSE), varList(varList) {}

/*************************************************/

void FileName::resolve() const
{
  try
  {
    if(resolved)
      return;
    nameParsed = StringParser::parse(nameUnparsed, varList);
    resolved = TRUE;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/*************************************************/

FileName FileName::append(const FileName &fileName) const
{
  resolve();
  if(empty())
    return fileName;
  std::string separator;
  if((nameParsed.back() != '/') && (nameParsed.back() != '\\'))
#ifdef _WIN32
    separator = '\\'; // = std::filesystem::path::preferred_separator;
#else
    separator = '/';
#endif
  return FileName(nameParsed+separator+fileName.str());
}

/*************************************************/

FileName FileName::fullExtension() const
{
  resolve();
  auto pos = nameParsed.rfind('.');
  if((pos == std::string::npos) || (pos+1 == nameParsed.size()))
    return FileName();
  std::string ext = String::upperCase(nameParsed.substr(pos+1));
  if((pos > 0) && (ext == "GZ" || ext == "Z"))
  {
    auto posNew = nameParsed.rfind('.', pos-1);
    if(posNew != std::string::npos)
      pos = posNew;
  }
  return FileName(nameParsed.substr(pos+1));
}

/*************************************************/

FileName FileName::typeExtension() const
{
  resolve();
  std::string ext = fullExtension();
  return FileName(ext.substr(0, ext.find('.')));
}

/*************************************************/

FileName FileName::packExtension() const
{
  resolve();
  std::string ext = fullExtension();
  auto pos = ext.find('.');
  if(pos == std::string::npos)
    return FileName();
  return FileName(ext.substr(pos+1));
}

/*************************************************/

FileName FileName::stripFullExtension() const
{
  resolve();
  auto pos = nameParsed.rfind("."+fullExtension().str());
  if(pos == std::string::npos)
    return *this;
  return FileName(nameParsed.substr(0, pos));
}

/*************************************************/

FileName FileName::replaceFullExtension(const std::string &text) const
{
  resolve();
  return FileName(stripFullExtension().str() + ((text.empty() || String::startsWith(text, ".")) ? "" : ".") + text);
}

/*************************************************/

FileName FileName::directory() const
{
  resolve();
  auto pos = nameParsed.find_last_of("/\\");
  if(pos == std::string::npos)
    return FileName();
  return FileName(nameParsed.substr(0, pos+1));
}

/*************************************************/

FileName FileName::stripDirectory() const
{
  resolve();
  std::string::size_type pos = nameParsed.find_last_of("/\\");
  if(pos == std::string::npos)
    return *this;
  return FileName(nameParsed.substr(pos+1));
}

/*************************************************/

FileName FileName::appendBaseName(const std::string &text) const
{
  resolve();
  auto pos = nameParsed.rfind("."+fullExtension().str());
  if(pos == std::string::npos)
    return FileName(nameParsed+text);
  std::string tmp = nameParsed;
  return FileName(tmp.insert(pos, text));
}

/*************************************************/

FileName FileName::operator()(const VariableList &varList) const
{
  try
  {
    VariableList varList2 = this->varList;
    varList2 += varList;
    return FileName(StringParser::parse(nameUnparsed, varList2));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/*************************************************/
