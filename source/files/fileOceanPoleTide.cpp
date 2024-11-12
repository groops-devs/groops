/***********************************************/
/**
* @file fileOceanPoleTide.cpp
*
* @brief Read/write OceanPoleTide.
*
* @author Torsten Mayer-Guerr
* @date 2017-12-09
*
*/
/***********************************************/

#define DOCSTRING_FILEFORMAT_OceanPoleTide

#include "base/import.h"
#include "base/sphericalHarmonics.h"
#include "inputOutput/fileArchive.h"
#include "files/fileFormatRegister.h"
#include "files/fileOceanPoleTide.h"

GROOPS_REGISTER_FILEFORMAT(OceanPoleTide, FILE_OCEANPOLETIDE_TYPE)

/***********************************************/

void writeFileOceanPoleTide(const FileName &fileName, const SphericalHarmonics &harmReal, const SphericalHarmonics &harmImag)
{
  try
  {
    OutFileArchive file(fileName, FILE_OCEANPOLETIDE_TYPE, FILE_OCEANPOLETIDE_VERSION);
    file<<nameValue("harmonicsReal",      harmReal);
    file<<nameValue("harmonicsImaginary", harmImag);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void readFileOceanPoleTide(const FileName &fileName, SphericalHarmonics &harmReal, SphericalHarmonics &harmImag)
{
  try
  {
    InFileArchive file(fileName, FILE_OCEANPOLETIDE_TYPE, FILE_OCEANPOLETIDE_VERSION);
    file>>nameValue("harmonicsReal",      harmReal);
    file>>nameValue("harmonicsImaginary", harmImag);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
