/***********************************************/
/**
* @file fileMeanPolarMotion.cpp
*
* @brief Read/write MeanPolarMotion.
*
* @author Torsten Mayer-Guerr
* @date 2017-12-09
*
*/
/***********************************************/

#define DOCSTRING_FILEFORMAT_MeanPolarMotion

#include "base/import.h"
#include "inputOutput/fileArchive.h"
#include "files/fileFormatRegister.h"
#include "files/fileMeanPolarMotion.h"

GROOPS_REGISTER_FILEFORMAT(MeanPolarMotion, FILE_MEANPOLARMOTION_TYPE)

/***********************************************/

void MeanPolarMotion::compute(const Time &time, Double &xBar, Double &yBar) const
{
  try
  {
    xBar = 0;
    yBar = 0;
    if(timeStart.size() == 0)
      return;

    // find correct interval
    if(time<timeStart.at(0))
      throw(Exception("mean pole model is not valid before "+timeStart.at(0).dateTimeStr()));
    UInt idxInterval = 0;
    while((idxInterval+1 < timeStart.size()) && (time >= timeStart.at(idxInterval+1)))
      idxInterval++;

    // mean pole (IERS2010, eq. (7.25))
    Double tn   = 1;
    const Double t = (time.mjd()-J2000)/365.25; // years after J2000
    for(UInt n=0; n<=degree.at(idxInterval); n++)
    {
      xBar += tn * xp.at(idxInterval)(n);
      yBar += tn * yp.at(idxInterval)(n);
      tn *= t;
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

template<> void save(OutArchive &ar, const MeanPolarMotion &x)
{
  const UInt intervalCount = x.timeStart.size();
  ar<<nameValue("intervalCount", intervalCount);
  if(!intervalCount)
    return;

  //comment
  const UInt maxDegree = *std::max_element(x.degree.begin(), x.degree.end());
  std::string str = "start time [MJD]            degree";
  for(UInt n=0; n<=maxDegree; n++)
    str += "  xp [arcsec/year^"+n%"%i"s+"]      ";
  for(UInt n=0; n<=maxDegree; n++)
    str += "  yp [arcsec/year^"+n%"%i"s+"]      ";
  ar.comment(str);
  ar.comment(std::string(str.size(), '='));

  for(UInt idxInterval=0; idxInterval<intervalCount; idxInterval++)
  {
    ar<<beginGroup("interval");
    ar<<nameValue("timeStart", x.timeStart.at(idxInterval));
    ar<<nameValue("degree",    x.degree.at(idxInterval));
    for(UInt n=0; n<=x.degree.at(idxInterval); n++)
      ar<<nameValue("xp", x.xp.at(idxInterval)(n));
    for(UInt n=0; n<=x.degree.at(idxInterval); n++)
      ar<<nameValue("yp", x.yp.at(idxInterval)(n));
    ar<<endGroup("interval");
  }
}

/***********************************************/

template<> void load(InArchive &ar, MeanPolarMotion &x)
{
  UInt intervalCount;
  ar>>nameValue("intervalCount", intervalCount);
  x.timeStart.resize(intervalCount);
  x.degree.resize(intervalCount);
  x.xp.resize(intervalCount);
  x.yp.resize(intervalCount);
  for(UInt idxInterval=0; idxInterval<intervalCount; idxInterval++)
  {
    ar>>beginGroup("interval");
    ar>>nameValue("timeStart", x.timeStart.at(idxInterval));
    ar>>nameValue("degree",    x.degree.at(idxInterval));
    x.xp.at(idxInterval) = Vector(x.degree.at(idxInterval)+1);
    x.yp.at(idxInterval) = Vector(x.degree.at(idxInterval)+1);
    for(UInt n=0; n<=x.degree.at(idxInterval); n++)
      ar>>nameValue("xp", x.xp.at(idxInterval)(n));
    for(UInt n=0; n<=x.degree.at(idxInterval); n++)
      ar>>nameValue("yp", x.yp.at(idxInterval)(n));
    ar>>endGroup("interval");
  }
}

/***********************************************/

void writeFileMeanPolarMotion(const FileName &fileName, const MeanPolarMotion &x)
{
  try
  {
    OutFileArchive file(fileName, FILE_MEANPOLARMOTION_TYPE);
    file<<nameValue("meanPolarMotion", x);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void readFileMeanPolarMotion(const FileName &fileName, MeanPolarMotion &x)
{
  try
  {
    InFileArchive file(fileName, FILE_MEANPOLARMOTION_TYPE);
    file>>nameValue("meanPolarMotion", x);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
