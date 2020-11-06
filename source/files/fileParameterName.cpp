/***********************************************/
/**
* @file fileParameterName.cpp
*
* @brief Read/write ParameterName.
*
* @author Torsten Mayer-Guerr
* @date 2017-07-31
*
*/
/***********************************************/

#define DOCSTRING_FILEFORMAT_ParameterName

#include "base/import.h"
#include "base/parameterName.h"
#include "inputOutput/fileArchive.h"
#include "files/fileFormatRegister.h"
#include "files/fileParameterName.h"

GROOPS_REGISTER_FILEFORMAT(ParameterName, FILE_PARAMETERNAME_TYPE)

/***********************************************/

template<> void save(OutArchive &ar, const ParameterName &x)
{
  save(ar, x.str());
}

/***********************************************/

template<> void load(InArchive  &ar, ParameterName &x)
{
  std::string str;
  load(ar, str);
  x = ParameterName::fromStr(str);
}

/***********************************************/

void writeFileParameterName(const FileName &fileName, const std::vector<ParameterName> &x)
{
  try
  {
    OutFileArchive file(fileName, FILE_PARAMETERNAME_TYPE);
    file.comment("object:type:temporal:interval");
    file.comment("=============================");
    file<<nameValue("parameterName", x);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void readFileParameterName(const FileName &fileName, std::vector<ParameterName> &x)
{
  try
  {
    InFileArchive file(fileName, FILE_PARAMETERNAME_TYPE);
    file>>nameValue("parameterName", x);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
