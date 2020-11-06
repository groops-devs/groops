/***********************************************/
/**
* @file fileGnssAntennaDefinition.cpp
*
* @brief Antenna center variations.
*
* @author Torsten Mayer-Guerr
* @author Sebastian Strasser
* @date 2012-11-22
*
*/
/***********************************************/

#define DOCSTRING_FILEFORMAT_GnssAntennaDefinition

#include "base/import.h"
#include "base/gnssType.h"
#include "inputOutput/fileArchive.h"
#include "inputOutput/logging.h"
#include "files/fileFormatRegister.h"
#include "files/fileGnssAntennaDefinition.h"

GROOPS_REGISTER_FILEFORMAT(GnssAntennaDefinition, "gnssAntennaDefinition")

/***********************************************/

Vector GnssAntennaDefinition::antennaVariations(Angle azimut, Angle elevation, const std::vector<GnssType> &types, NoPatternFoundAction noPatternFoundAction) const
{
  try
  {
    Vector acv(types.size());
    for(UInt idType=0; idType<types.size(); idType++)
      if((types.at(idType) == GnssType::PHASE) || (types.at(idType) == GnssType::RANGE))
      {
        const UInt idPattern = findAntennaPattern(types.at(idType), noPatternFoundAction);
        acv(idType) = (idPattern == NULLINDEX) ? NAN_EXPR : pattern.at(idPattern).antennaVariations(azimut, elevation);
      }

    return acv;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

UInt GnssAntennaDefinition::findAntennaPattern(const GnssType &type, NoPatternFoundAction noPatternFoundAction) const
{
  try
  {
    auto iter = std::find_if(pattern.begin(), pattern.end(), [&](const GnssAntennaPattern &p) {return p.type == type;});
    if(iter == pattern.end()) // no matching pattern found
    {
      if(noPatternFoundAction == IGNORE_OBSERVATION)
        return NULLINDEX;
      if(noPatternFoundAction == USE_NEAREST_FREQUENCY)
      {
        // TODO: there should be the following priorities: phase? nearest phase; code? nearest code or if no code then nearest phase
        iter = std::min_element(pattern.begin(), pattern.end(), [&](const GnssAntennaPattern &p1, const GnssAntennaPattern &p2)
        {
          GnssType t1 = p1.type; if(t1 == GnssType::GLONASS) t1.setFrequencyNumber(0);
          GnssType t2 = p2.type; if(t2 == GnssType::GLONASS) t2.setFrequencyNumber(0);
          return std::fabs(t1.frequency() - type.frequency()) < std::fabs(t2.frequency() - type.frequency());
        });
      }
      if(noPatternFoundAction == THROW_EXCEPTION || iter == pattern.end()) // still no matching pattern found
        throw(Exception("no pattern found for "+type.str()));
    }

    return static_cast<UInt>(std::distance(pattern.begin(), iter));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

UInt GnssAntennaDefinition::find(const std::vector<GnssAntennaDefinitionPtr> &antennaList, const std::string &name, const std::string &serial, const std::string radome)
{
  try
  {
    for(UInt i=antennaList.size(); i-->0;) // backwards search
      if((name == antennaList.at(i)->name) || antennaList.at(i)->name.empty())
        if((serial == antennaList.at(i)->serial) || antennaList.at(i)->serial.empty())
          if((radome == antennaList.at(i)->radome) || antennaList.at(i)->radome.empty())
            return i;

    return NULLINDEX;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/

Double GnssAntennaPattern::antennaVariations(Angle azimut, Angle elevation, Bool applyOffset) const
{
  try
  {
    Double tauA  = std::fmod(Double(azimut)/(2*PI)+1,1) * pattern.rows();
    Double tauZ  = (PI/2-Double(elevation))/Double(dZenit);
    const UInt idxA1 = static_cast<UInt>(std::floor(tauA));
    const UInt idxA2 = (idxA1+1)%pattern.rows();
    const UInt idxZ1 = std::min(static_cast<UInt>(std::floor(tauZ)), pattern.columns()-1);
    const UInt idxZ2 = std::min(idxZ1+1, pattern.columns()-1);
    tauA -= std::floor(tauA);
    tauZ -= std::floor(tauZ);

    // bilinear interpolation
    Double acv = (1-tauA) * (1-tauZ) * pattern(idxA1, idxZ1)
               + (1-tauA) *  (tauZ)  * pattern(idxA1, idxZ2)
               +  (tauA)  * (1-tauZ) * pattern(idxA2, idxZ1)
               +  (tauA)  *  (tauZ)  * pattern(idxA2, idxZ2);

    if(applyOffset)
      acv -= inner(offset, polar(azimut, elevation, 1.));

    return acv;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}


/***********************************************/
/***********************************************/

template<> void save(OutArchive &ar, const GnssAntennaDefinition &x)
{
  try
  {
    ar<<nameValue("name",      x.name);
    ar<<nameValue("serial",    x.serial);
    ar<<nameValue("radome",    x.radome);
    ar<<nameValue("comment",   x.comment);
    ar<<nameValue("pattern",   x.pattern);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

template<> void load(InArchive &ar, GnssAntennaDefinition  &x)
{
  try
  {
    ar>>nameValue("name",      x.name);
    ar>>nameValue("serial",    x.serial);
    ar>>nameValue("radome",    x.radome);
    ar>>nameValue("comment",   x.comment);
    ar>>nameValue("pattern",   x.pattern);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/

template<> void save(OutArchive &ar, const GnssAntennaDefinitionPtr &x)
{
  try
  {
    save(ar, *x.get());
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

template<> void load(InArchive  &ar, GnssAntennaDefinitionPtr &x)
{
  try
  {
    if(!x)
      x = GnssAntennaDefinitionPtr(new GnssAntennaDefinition);
    load(ar, *x.get());
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

template<> void save(OutArchive &ar, const GnssAntennaPattern &x)
{
  try
  {
    ar<<nameValue("type",     x.type);
    ar<<nameValue("offset",   x.offset);
    ar<<nameValue("dZenit",   x.dZenit);
    ar<<nameValue("pattern",  x.pattern);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

template<> void load(InArchive  &ar, GnssAntennaPattern &x)
{
  try
  {
    ar>>nameValue("type",     x.type);
    ar>>nameValue("offset",   x.offset);
    ar>>nameValue("dZenit",   x.dZenit);
    ar>>nameValue("pattern",  x.pattern);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/

void writeFileGnssAntennaDefinition(const FileName &fileName, const GnssAntennaDefinition &x)
{
  try
  {
    OutFileArchive file(fileName, FILE_GNSSANTENNADEFINITION_TYPE);
    file<<nameValue("antennaCount", 1);
    file<<nameValue("antenna", x);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void writeFileGnssAntennaDefinition(const FileName &fileName, const std::vector<GnssAntennaDefinitionPtr> &x)
{
  try
  {
    OutFileArchive file(fileName, FILE_GNSSANTENNADEFINITION_TYPE);
    file<<nameValue("antennaCount", x.size());
    for(UInt i=0; i<x.size(); i++)
      file<<nameValue("antenna", x.at(i));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void readFileGnssAntennaDefinition(const FileName &fileName, GnssAntennaDefinition &x)
{
  try
  {
    InFileArchive file(fileName, FILE_GNSSANTENNADEFINITION_TYPE);
    if(file.version() < 20190304)
      throw(Exception(fileName.str()+": old GnssAntennaDefinition file, definition of reference frames changed"));

    UInt count;
    file>>nameValue("antennaCount", count);
    if(count>1)
      logWarning<<fileName<<" contain more than one antenna, only the first is used"<<Log::endl;
    file>>nameValue("antenna", x);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void readFileGnssAntennaDefinition(const FileName &fileName, std::vector<GnssAntennaDefinitionPtr> &x)
{
  try
  {
    InFileArchive file(fileName, FILE_GNSSANTENNADEFINITION_TYPE);
    if(file.version() < 20190304)
      throw(Exception(fileName.str()+": old GnssAntennaDefinition file, definition of reference frames changed"));

    UInt count;
    file>>nameValue("antennaCount", count);
    x.resize(count);
    for(UInt i=0; i<x.size(); i++)
      file>>nameValue("antenna", x.at(i));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/
