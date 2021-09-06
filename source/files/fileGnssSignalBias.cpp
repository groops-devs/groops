/***********************************************/
/**
* @file fileGnssSignalBias.cpp
*
* @brief Code/Phase biases.
*
* @author Torsten Mayer-Guerr
* @date 2013-08-11
*
*/
/***********************************************/

#define DOCSTRING_FILEFORMAT_GnssSignalBias

#include "base/import.h"
#include "inputOutput/fileArchive.h"
#include "files/fileFormatRegister.h"
#include "files/fileGnssSignalBias.h"

GROOPS_REGISTER_FILEFORMAT(GnssSignalBias, FILE_GNSSSIGNALBIAS_TYPE)

/***********************************************/

Vector GnssSignalBias::compute(const std::vector<GnssType> &types) const
{
  try
  {
    Vector b(types.size());
    UInt idx;
    for(UInt idType=0; idType<types.size(); idType++)
      if(types.at(idType).isInList(this->types, idx))
        b(idType) = biases.at(idx);
    return b;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/

template<> void save(OutArchive &ar, const GnssSignalBias &x)
{
  try
  {
    ar<<nameValue("count", x.types.size());
    ar.comment("type   bias [m]");
    ar.comment("===============================");
    for(UInt i=0; i<x.types.size(); i++)
    {
      ar.endLine();
      ar<<nameValue("type", x.types.at(i));
      ar<<nameValue("bias", x.biases.at(i));
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

template<> void load(InArchive  &ar, GnssSignalBias &x)
{
  try
  {
    UInt count;
    ar>>nameValue("count", count);
    x.types.resize(count);
    x.biases.resize(count);
    for(UInt i=0; i<count; i++)
    {
      ar>>nameValue("type", x.types.at(i));
      ar>>nameValue("bias", x.biases.at(i));
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/

void writeFileGnssSignalBias(const FileName &fileName, const GnssSignalBias &x)
{
  try
  {
    OutFileArchive file(fileName, FILE_GNSSSIGNALBIAS_TYPE);
    file<<nameValue("signalBias", x);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void readFileGnssSignalBias(const FileName &fileName, GnssSignalBias &x)
{
  try
  {
    InFileArchive file(fileName, FILE_GNSSSIGNALBIAS_TYPE);
    file>>nameValue("signalBias", x);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
