/***********************************************/
/**
* @file fileDoodsonHarmonic.cpp
*
* @brief Read/write harmonics coded with doodson frequencies.
*
* @author Torsten Mayer-Guerr
* @date 2010-06-10
*
*/
/***********************************************/

#define DOCSTRING_FILEFORMAT_DoodsonHarmonic

#include "base/import.h"
#include "base/doodson.h"
#include "inputOutput/fileArchive.h"
#include "files/fileFormatRegister.h"
#include "files/fileDoodsonHarmonic.h"

GROOPS_REGISTER_FILEFORMAT(DoodsonHarmonic, FILE_DOODSONHARMONIC_TYPE)

/***********************************************/

template<> void save(OutArchive &ar, const DoodsonHarmonic &x)
{
  UInt count = x.doodson.size();
  ar<<nameValue("GM", x.GM);
  ar<<nameValue("R",  x.R);
  ar<<nameValue("constituentsCount", count);
  for(UInt i=0; i<count; i++)
  {
    ar<<beginGroup("constituent");
    ar<<nameValue("doodson", x.doodson.at(i));
    ar<<nameValue("cnmCos",  x.cnmCos.at(i));
    ar<<nameValue("snmCos",  x.snmCos.at(i));
    ar<<nameValue("cnmSin",  x.cnmSin.at(i));
    ar<<nameValue("snmSin",  x.snmSin.at(i));
    ar<<endGroup("constituent");
  }
}

/***********************************************/

template<> void load(InArchive &ar, DoodsonHarmonic &x)
{
  UInt count;

  ar>>nameValue("GM", x.GM);
  ar>>nameValue("R",  x.R);
  ar>>nameValue("constituentsCount", count);

  x.doodson.resize(count);
  x.cnmCos.resize(count);
  x.cnmSin.resize(count);
  x.snmCos.resize(count);
  x.snmSin.resize(count);

  for(UInt i=0; i<count; i++)
  {
    ar>>beginGroup("constituent");
    ar>>nameValue("doodson", x.doodson.at(i));
    ar>>nameValue("cnmCos",  x.cnmCos.at(i));
    ar>>nameValue("snmCos",  x.snmCos.at(i));
    ar>>nameValue("cnmSin",  x.cnmSin.at(i));
    ar>>nameValue("snmSin",  x.snmSin.at(i));
    ar>>endGroup("constituent");
  }
}

/***********************************************/

void writeFileDoodsonHarmonic(const FileName &fileName, const DoodsonHarmonic &x)
{
  try
  {
    OutFileArchive file(fileName, FILE_DOODSONHARMONIC_TYPE);
    file<<nameValue("doodsonHarmonic", x);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void readFileDoodsonHarmonic(const FileName &fileName, DoodsonHarmonic &x)
{
  try
  {
    InFileArchive file(fileName, FILE_DOODSONHARMONIC_TYPE);
    file>>nameValue("doodsonHarmonic", x);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
