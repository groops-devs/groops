/***********************************************/
/**
* @file fileDoodsonEarthOrientationParameter.cpp
*
* @brief Read/write high frequent EOPs.
*
* @author Torsten Mayer-Guerr
* @date 2019-05-15
*
*/
/***********************************************/

#define DOCSTRING_FILEFORMAT_DoodsonEarthOrientationParameter

#include "base/import.h"
#include "base/doodson.h"
#include "inputOutput/fileArchive.h"
#include "files/fileFormatRegister.h"
#include "files/fileDoodsonEarthOrientationParameter.h"

GROOPS_REGISTER_FILEFORMAT(DoodsonEarthOrientationParameter, FILE_DOODSONEARTHORIENTATIONPARAMETER_TYPE)

/***********************************************/

template<> void save(OutArchive &ar, const DoodsonEop &x)
{
  ar<<nameValue("constituentsCount", x.doodson.size());
  ar.comment("dood.   xpCos [arcsec]            xpSin [arcsec]            ypCos [arcsec]            ypSin [arcsec]            ut1Cos [sec]              ut1Sin [sec]              lodCos [sec]              lodSin [sec]             name");
  ar.comment("===========================================================================================================================================================================================================================");
  for(UInt i=0; i<x.doodson.size(); i++)
  {
    ar<<beginGroup("constituent");
    ar<<nameValue("doodson", x.doodson.at(i));
    ar<<nameValue("xpCos",   x.coeff(i, 0));
    ar<<nameValue("xpSin",   x.coeff(i, 1));
    ar<<nameValue("ypCos",   x.coeff(i, 2));
    ar<<nameValue("ypSin",   x.coeff(i, 3));
    ar<<nameValue("ut1Cos",  x.coeff(i, 4));
    ar<<nameValue("ut1Sin",  x.coeff(i, 5));
    ar<<nameValue("lodCos",  x.coeff(i, 6));
    ar<<nameValue("lodSin",  x.coeff(i, 7));
    ar<<nameValue("name",    (x.doodson.at(i).name() != x.doodson.at(i).code()) ? x.doodson.at(i).name() : "");
    ar<<endGroup("constituent");
  }
}

/***********************************************/

template<> void load(InArchive &ar, DoodsonEop &x)
{
  UInt constituentsCount;
  ar>>nameValue("constituentsCount", constituentsCount);
  x.doodson.resize(constituentsCount);
  x.coeff = Matrix(constituentsCount, 8);
  for(UInt i=0; i<constituentsCount; i++)
  {
    std::string name;
    ar>>beginGroup("constituent");
    ar>>nameValue("doodson", x.doodson.at(i));
    ar>>nameValue("xpCos",   x.coeff(i, 0));
    ar>>nameValue("xpSin",   x.coeff(i, 1));
    ar>>nameValue("ypCos",   x.coeff(i, 2));
    ar>>nameValue("ypSin",   x.coeff(i, 3));
    ar>>nameValue("ut1Cos",  x.coeff(i, 4));
    ar>>nameValue("ut1Sin",  x.coeff(i, 5));
    ar>>nameValue("lodCos",  x.coeff(i, 6));
    ar>>nameValue("lodSin",  x.coeff(i, 7));
    ar>>nameValue("name",    name);
    ar>>endGroup("constituent");
  }
}

/***********************************************/

void writeFileDoodsonEarthOrientationParameter(const FileName &fileName, const DoodsonEop &x)
{
  try
  {
    OutFileArchive file(fileName, FILE_DOODSONEARTHORIENTATIONPARAMETER_TYPE);
    file<<nameValue("doodsonEarthOrientationParameter", x);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void readFileDoodsonEarthOrientationParameter(const FileName &fileName, DoodsonEop &x)
{
  try
  {
    InFileArchive file(fileName, FILE_DOODSONEARTHORIENTATIONPARAMETER_TYPE);
    file>>nameValue("doodsonEarthOrientationParameter", x);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
