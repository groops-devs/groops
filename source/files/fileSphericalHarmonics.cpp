/***********************************************/
/**
* @file fileSphericalHarmonics.cpp
*
* @brief Read/write SphericalHarmonics.
*
* @author Torsten Mayer-Guerr
* @date 2005-01-14
*
*/
/***********************************************/

#define DOCSTRING_FILEFORMAT_PotentialCoefficients

#include "base/import.h"
#include "base/sphericalHarmonics.h"
#include "inputOutput/fileArchive.h"
#include "files/fileFormatRegister.h"
#include "files/fileSphericalHarmonics.h"

GROOPS_REGISTER_FILEFORMAT(PotentialCoefficients, FILE_POTENTIALCOEFFICIENTS_TYPE)

/***********************************************/

void writeFileSphericalHarmonics(const FileName &fileName, const SphericalHarmonics &x)
{
  try
  {
    OutFileArchive file(fileName, FILE_POTENTIALCOEFFICIENTS_TYPE, FILE_POTENTIALCOEFFICIENTS_VERSION);
    file<<nameValue("potentialCoefficients", x);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void readFileSphericalHarmonics(const FileName &fileName, SphericalHarmonics &x)
{
  try
  {
    InFileArchive file(fileName, FILE_POTENTIALCOEFFICIENTS_TYPE, FILE_POTENTIALCOEFFICIENTS_VERSION);
    file>>nameValue("potentialCoefficients", x);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
