/***********************************************/
/**
* @file fileGnssReceiverDefinition.cpp
*
* @brief GNSS receiver definition.
*
* @author Sebastian Strasser
* @date 2019-08-28
*/
/***********************************************/

#define DOCSTRING_FILEFORMAT_GnssReceiverDefinition

#include "base/import.h"
#include "base/gnssType.h"
#include "inputOutput/fileArchive.h"
#include "inputOutput/logging.h"
#include "files/fileFormatRegister.h"
#include "files/fileGnssReceiverDefinition.h"

GROOPS_REGISTER_FILEFORMAT(GnssReceiverDefinition, "gnssReceiverDefinition")

/***********************************************/

GnssReceiverDefinitionPtr GnssReceiverDefinition::find(const std::vector<GnssReceiverDefinitionPtr> &receivers, const std::string &name, const std::string &serial, const std::string version)
{
  try
  {
    for(UInt i=receivers.size(); i-->0;) // backwards search
      if((name == receivers.at(i)->name) || receivers.at(i)->name.empty())
        if((serial == receivers.at(i)->serial) || receivers.at(i)->serial.empty())
          if((version == receivers.at(i)->version) || receivers.at(i)->version.empty())
            return receivers.at(i);

    return GnssReceiverDefinitionPtr(nullptr);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/

template<> void save(OutArchive &ar, const GnssReceiverDefinition &x)
{
  try
  {
    ar<<nameValue("name",      x.name);
    ar<<nameValue("serial",    x.serial);
    ar<<nameValue("version",   x.version);
    ar<<nameValue("comment",   x.comment);
    ar<<nameValue("types",     x.types);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

template<> void load(InArchive &ar, GnssReceiverDefinition  &x)
{
  try
  {
    ar>>nameValue("name",      x.name);
    ar>>nameValue("serial",    x.serial);
    ar>>nameValue("version",   x.version);
    ar>>nameValue("comment",   x.comment);
    ar>>nameValue("types",     x.types);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/

template<> void save(OutArchive &ar, const GnssReceiverDefinitionPtr &x)
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

template<> void load(InArchive  &ar, GnssReceiverDefinitionPtr &x)
{
  try
  {
    if(!x)
      x = GnssReceiverDefinitionPtr(new GnssReceiverDefinition);
    load(ar, *x.get());
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/

void writeFileGnssReceiverDefinition(const FileName &fileName, const GnssReceiverDefinition &x)
{
  try
  {
    OutFileArchive file(fileName, FILE_GNSSRECEIVERDEFINITION_TYPE);
    file<<nameValue("receiverCount", 1);
    file<<nameValue("receiver", x);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void writeFileGnssReceiverDefinition(const FileName &fileName, const std::vector<GnssReceiverDefinitionPtr> &x)
{
  try
  {
    OutFileArchive file(fileName, FILE_GNSSRECEIVERDEFINITION_TYPE);
    file<<nameValue("receiverCount", x.size());
    for(UInt i=0; i<x.size(); i++)
      file<<nameValue("receiver", x.at(i));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void readFileGnssReceiverDefinition(const FileName &fileName, GnssReceiverDefinition &x)
{
  try
  {
    InFileArchive file(fileName, FILE_GNSSRECEIVERDEFINITION_TYPE);

    UInt count;
    file>>nameValue("receiverCount", count);
    if(count>1)
      logWarning<<fileName<<" contains more than one receiver, only the first is used"<<Log::endl;
    file>>nameValue("receiver", x);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void readFileGnssReceiverDefinition(const FileName &fileName, std::vector<GnssReceiverDefinitionPtr> &x)
{
  try
  {
    InFileArchive file(fileName, FILE_GNSSRECEIVERDEFINITION_TYPE);

    UInt count;
    file>>nameValue("receiverCount", count);
    x.resize(count);
    for(UInt i=0; i<x.size(); i++)
      file>>nameValue("receiver", x.at(i));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/
