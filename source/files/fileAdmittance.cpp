/***********************************************/
/**
* @file fileAdmittance.cpp
*
* @brief Read/write admittance.
*
* @author Torsten Mayer-Guerr
* @date 2013-02-08
*
*/
/***********************************************/

#define DOCSTRING_FILEFORMAT_Admittance

#include "base/import.h"
#include "base/doodson.h"
#include "inputOutput/fileArchive.h"
#include "files/fileFormatRegister.h"
#include "files/fileAdmittance.h"

GROOPS_REGISTER_FILEFORMAT(Admittance, FILE_ADMITTANCE_TYPE)

/***********************************************/

void writeFileAdmittance(const FileName &fileName, const Admittance &admittance)
{
  try
  {
    OutFileArchive file(fileName, FILE_ADMITTANCE_TYPE, FILE_ADMITTANCE_VERSION);
    file<<nameValue("major",      admittance.doodsonMajor);
    file<<nameValue("minor",      admittance.doodsonMinor);
    file<<nameValue("admittance", admittance.admittance);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void readFileAdmittance(const FileName &fileName, Admittance &admittance)
{
  try
  {
    InFileArchive file(fileName, FILE_ADMITTANCE_TYPE, FILE_ADMITTANCE_VERSION);
    file>>nameValue("major",      admittance.doodsonMajor);
    file>>nameValue("minor",      admittance.doodsonMinor);
    file>>nameValue("admittance", admittance.admittance);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
