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

Vector GnssSignalBias::compute(const std::vector<GnssType> &type) const
{
  try
  {
    Vector b(type.size());
    for(UInt idType=0; idType<type.size(); idType++)
    {
      const UInt i = GnssType::index(this->type, type.at(idType));
      if(i != NULLINDEX)
        b(idType) += bias.at(i);
    }
    return b;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssSignalBias::applyTecBias(Double tecBias)
{
  try
  {
    for(UInt idType=0; idType<type.size(); idType++)
    {
      GnssType t = type.at(idType);
      if(t.frequencyNumber() == 9999)
        t.setFrequencyNumber(0);
      if(type.at(idType) == GnssType::PHASE)
        bias.at(idType) += -Ionosphere::Ap/pow(t.frequency(), 2) * tecBias;
      if(type.at(idType) == GnssType::RANGE)
        bias.at(idType) +=  Ionosphere::Ap/pow(t.frequency(), 2) * tecBias;
    }
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
    ar<<nameValue("count", x.type.size());
    ar.comment("type   bias [m]");
    ar.comment("===============================");
    for(UInt i=0; i<x.type.size(); i++)
    {
      ar.endLine();
      ar<<nameValue("type", x.type.at(i));
      ar<<nameValue("bias", x.bias.at(i));
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
    x.type.resize(count);
    x.bias.resize(count);
    for(UInt i=0; i<count; i++)
    {
      ar>>nameValue("type", x.type.at(i));
      ar>>nameValue("bias", x.bias.at(i));
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
