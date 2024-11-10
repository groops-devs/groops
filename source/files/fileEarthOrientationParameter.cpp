/***********************************************/
/**
* @file fileEarthOrientationParameter.cpp
*
* @brief Read/write EarthOrientationParameter.
*
* @author Torsten Mayer-Guerr
* @date 2017-12-09
*
*/
/***********************************************/

#define DOCSTRING_FILEFORMAT_EarthOrientationParameter

#include "base/import.h"
#include "inputOutput/fileArchive.h"
#include "files/fileFormatRegister.h"
#include "files/fileEarthOrientationParameter.h"

GROOPS_REGISTER_FILEFORMAT(EarthOrientationParameter, FILE_EARTHORIENTATIONPARAMETER_TYPE)

/***********************************************/

void writeFileEarthOrientationParameter(const FileName &fileName, const const_MatrixSlice &EOP)
{
  try
  {
    OutFileArchive file(fileName, FILE_EARTHORIENTATIONPARAMETER_TYPE, FILE_EARTHORIENTATIONPARAMETER_VERSION);
    file<<nameValue("count", EOP.rows());
    file.comment("UTC [MJD]                 xp [arcsec]               yp [arcsec]               deltUT [sec]              LOD [sec]                 dX [arcsec]               dY [arcsec]             ");
    file.comment("====================================================================================================================================================================================");
    for(UInt i=0; i<EOP.rows(); i++)
    {
      file<<beginGroup("eop");
      file<<nameValue("mjd",     EOP(i,0));
      file<<nameValue("xp",      EOP(i,1));
      file<<nameValue("yp",      EOP(i,2));
      file<<nameValue("deltaUT", EOP(i,3));
      file<<nameValue("lod",     EOP(i,4));
      file<<nameValue("dX",      EOP(i,5));
      file<<nameValue("dY",      EOP(i,6));
      file<<endGroup("eop");
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void readFileEarthOrientationParameter(const FileName &fileName, Matrix &EOP)
{
  try
  {
    InFileArchive file(fileName, FILE_EARTHORIENTATIONPARAMETER_TYPE, FILE_EARTHORIENTATIONPARAMETER_VERSION);
    UInt count;
    file>>nameValue("count", count);
    EOP = Matrix(count, 7);
    for(UInt i=0; i<count; i++)
    {
      file>>beginGroup("eop");
      file>>nameValue("mjd",     EOP(i,0));
      file>>nameValue("xp",      EOP(i,1));
      file>>nameValue("yp",      EOP(i,2));
      file>>nameValue("deltaUT", EOP(i,3));
      file>>nameValue("lod",     EOP(i,4));
      file>>nameValue("dX",      EOP(i,5));
      file>>nameValue("dY",      EOP(i,6));
      file>>endGroup("eop");
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
