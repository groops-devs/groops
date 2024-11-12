/***********************************************/
/**
* @file fileArcList.cpp
*
* @brief Read/write Arc list.
*
* @author Torsten Mayer-Guerr
* @date 2013-02-07
*
*/
/***********************************************/

#define DOCSTRING_FILEFORMAT_ArcList

#include "base/import.h"
#include "inputOutput/fileArchive.h"
#include "files/fileFormatRegister.h"
#include "files/fileArcList.h"

GROOPS_REGISTER_FILEFORMAT(ArcList, "arcList")

/***********************************************/

void writeFileArcList(const FileName &fileName, const std::vector<UInt> &arcsInterval, const std::vector<Time> &timesInterval)
{
  try
  {
    OutFileArchive file(fileName, FILE_ARCLIST_TYPE, FILE_ARCLIST_VERSION);
    file<<nameValue("intervalCount", arcsInterval.size());
    file.comment("time [MJD]               first arc");
    file.comment("==================================");
    for(UInt i=0; i<arcsInterval.size(); i++)
    {
      file<<beginGroup("interval");
      file<<nameValue("time", timesInterval.at(i));
      file<<nameValue("arc",  arcsInterval.at(i));
      file<<endGroup("interval");
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void readFileArcList(const FileName &fileName, std::vector<UInt> &arcsInterval, std::vector<Time> &timesInterval)
{
  try
  {
    InFileArchive file(fileName, FILE_ARCLIST_TYPE, FILE_ARCLIST_VERSION);
    if(file.version() < 20200123)
    {
      file>>nameValue("arcList", arcsInterval);
      file>>nameValue("times",   timesInterval);
      return;
    }
    UInt count;
    file>>nameValue("intervalCount", count);
    arcsInterval.resize(count);
    timesInterval.resize(count);
    for(UInt i=0; i<count; i++)
    {
      file>>beginGroup("interval");
      file>>nameValue("time", timesInterval.at(i));
      file>>nameValue("arc",  arcsInterval.at(i));
      file>>endGroup("interval");
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
