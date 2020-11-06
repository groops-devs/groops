/***********************************************/
/**
* @file earthRotationFile.cpp
*
* @brief Interpolated values from file.
* @see EarthRotation
*
* @author Torsten Mayer-Guerr
* @date 2017-05-26
*
*/
/***********************************************/

#include "base/import.h"
#include "base/planets.h"
#include "config/config.h"
#include "files/fileMatrix.h"
#include "classes/earthRotation/earthRotation.h"
#include "classes/earthRotation/earthRotationFile.h"

/***********************************************/

EarthRotationFile::EarthRotationFile(Config &config)
{
  try
  {
    FileName fileNameEOP;
    UInt     interpolationDegree;

    readConfig(config, "inputfileEOP",        fileNameEOP,         Config::MUSTSET,  "{groopsDataDir}/earthRotation/timeSeries_EOP_rapid_IAU2000_desai.dat", "");
    readConfig(config, "interpolationDegree", interpolationDegree, Config::DEFAULT,  "3", "for polynomial interpolation");
    if(isCreateSchema(config)) return;

    readFileMatrix(fileNameEOP, EOP);

    times.resize(EOP.rows());
    for(UInt i=0; i<times.size(); i++)
      times.at(i) = mjd2time(EOP(i,0));
    if(!isRegular(times))
      throw(Exception("Time series must be given with constant sampling"));

    // UT1-UTC => UT1-GPS (avoid leap seconds jumps for interpolation)
    for(UInt i=0; i<times.size(); i++)
      EOP(i,4) -= (times.at(i)-timeGPS2UTC(times.at(i))).seconds();

    polynomial.init(interpolationDegree);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void EarthRotationFile::earthOrientationParameter(const Time &timeGPS, Double &xp, Double &yp, Double &sp,
                                                  Double &dUT1, Double &LOD, Double &X, Double &Y, Double &S) const
{
  try
  {
    if((timeGPS<times.at(0)) || (timeGPS>times.back()))
      throw(Exception("No EOPs available: "+timeGPS.dateTimeStr()));

    Matrix eop = polynomial.interpolate({timeGPS}, times, EOP);
    xp   = eop(0,1);
    yp   = eop(0,2);
    sp   = eop(0,3);
    dUT1 = eop(0,4) + (timeGPS-timeGPS2UTC(timeGPS)).seconds();
    LOD  = eop(0,5);
    X    = eop(0,6);
    Y    = eop(0,7);
    S    = eop(0,8);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
