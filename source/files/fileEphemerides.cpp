/***********************************************/
/**
* @file fileEphemerides.cpp
*
* @brief Ephemerides of sun, moon, and planets.
*
* @author Torsten Mayer-Guerr
* @date 2020-08-11
*
*/
/***********************************************/

#define DOCSTRING_FILEFORMAT_Ephemerides

#include "base/import.h"
#include "inputOutput/fileArchive.h"
#include "files/fileFormatRegister.h"
#include "files/fileEphemerides.h"

GROOPS_REGISTER_FILEFORMAT(Ephemerides, FILE_EPHEMERIDES_TYPE)

/***********************************************/

void InFileEphemerides::open(const FileName &name)
{
  try
  {
    close();
    if(name.empty())
      return;

    UInt intervalsCount, bodyCount;

    file.open(name, FILE_EPHEMERIDES_TYPE);
    file>>nameValue("bodyCount",      bodyCount);
    file>>nameValue("intervalsCount", intervalsCount);
    file>>nameValue("earthMoonRatio", earthMoonRatio);

    degree.resize(bodyCount);
    subintervals.resize(bodyCount);
    components.resize(bodyCount);

    for(UInt idBody=0; idBody<bodyCount; idBody++)
    {
      file>>beginGroup("body");
      file>>nameValue("subintervals", subintervals.at(idBody));
      file>>nameValue("components",   components.at(idBody));
      file>>nameValue("degree",       degree.at(idBody));
      file>>endGroup("body");
    }

    times.resize(intervalsCount+1);
    for(UInt i=0; i<times.size(); i++)
      file>>nameValue("time", times.at(i));

    coeff.resize(bodyCount);
    for(UInt idBody=0; idBody<bodyCount; idBody++)
      coeff.at(idBody).resize(subintervals.at(idBody), Matrix(components.at(idBody), degree.at(idBody)+1));

    // If we have a binary file, we can efficiently seek to the needed position in the file.
    if(file.canSeek())
      seekStart = file.position();
    seekSize = 0;
    idInterval = 0;
    readInterval(0);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW_EXTRA("filename=<"+name.str()+">", e)
  }
}

/***********************************************/

void InFileEphemerides::close()
{
  file.close();
}

/***********************************************/

void InFileEphemerides::readInterval(UInt idInterval_)
{
  try
  {
    // must restart?
    if(idInterval > idInterval_)
    {
      if(file.canSeek() && seekSize)
        idInterval = 0;
      else
        open(file.fileName());
    }

    while(idInterval <= idInterval_)
    {
      // Seek to appropriate interval.
      if(file.canSeek() && seekSize)
      {
        file.seek(seekStart + static_cast<std::streamoff>(idInterval_ * seekSize));
        idInterval = idInterval_;
      }

      file>>beginGroup("interval");
      for(UInt idBody=0; idBody<subintervals.size(); idBody++)
      {
        file>>beginGroup("body");
        for(UInt idSub=0; idSub<subintervals.at(idBody); idSub++)
        {
          file>>beginGroup("subinterval");
          for(UInt n=0; n<=degree.at(idBody); n++)
            for(UInt i=0; i<components.at(idBody); i++)
              file>>nameValue("coefficient", coeff.at(idBody).at(idSub)(i, n));
          file>>endGroup("subinterval");
        }
        file>>endGroup("body");
      }
      file>>endGroup("interval");

      if(file.canSeek() && (idInterval == 0))
        seekSize = file.position() - seekStart;
      idInterval++;
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW_EXTRA("filename=<"+file.fileName().str()+">", e)
  }
}

/***********************************************/

Matrix InFileEphemerides::interpolate(const Time &time, Planet planet)
{
  try
  {
    if(file.fileName().empty())
      throw(Exception("no file open"));

    // read new data block?
    if((time < times.at(idInterval-1)) || (time >= times.at(idInterval)))
    {
      if((time < times.front()) || (time >= times.back()))
        throw(Exception(time.dateTimeStr()+" outside interval ["+times.front().dateTimeStr()+", "+times.back().dateTimeStr()+")"));
      readInterval(static_cast<UInt>(std::floor((time-times.front()).mjd()/(times.back()-times.front()).mjd()*(times.size()-1))));
    }

    UInt idBody = 0;
    switch(planet)
    {
      case MERCURY:             idBody =  0; break;
      case VENUS:               idBody =  1; break;
      case EARTH:               idBody =  9; break; // moon relative to earth
      case MARS:                idBody =  3; break;
      case JUPITER:             idBody =  4; break;
      case SATURN:              idBody =  5; break;
      case URANUS:              idBody =  6; break;
      case NEPTUNE:             idBody =  7; break;
      case PLUTO:               idBody =  8; break;
      case MOON:                idBody =  9; break;  // moon relative to earth
      case SUN:                 idBody = 10; break;
      case SOLARBARYCENTER:     return Matrix(3, 2);
      case EARTHMOONBARYCENTER: idBody =  2; break;
      case NUTATION:            idBody = 11; break;
      case LIBRATION:           idBody = 12; break;
      case MOONFROMEARTH:       idBody =  9; break;
    }

    if(idBody >= subintervals.size() || subintervals.at(idBody) == 0)
      throw(Exception(static_cast<UInt>(planet)%"No data for requested planet id %i"s));

    Double dt = (times.at(idInterval)-times.at(idInterval-1)).seconds()/subintervals.at(idBody);
    Double t  = (time-times.at(idInterval-1)).seconds()/dt;
    const UInt idSub = static_cast<UInt>(std::floor(t));
    t = 2.*(t-idSub)-1.; // [-1, 1) in subinterval idSub

    // Chebyshev polynomials and derivative in second column
    Matrix p(degree.at(idBody)+1, 2);
    p(0, 0) = 1.;
    p(1, 0) = t;
    p(1, 1) = 2./dt;
    for(UInt n=2; n<p.rows(); n++)
    {
      p(n, 0) = 2*t*p(n-1, 0) - p(n-2, 0);
      p(n, 1) = 2*t*p(n-1, 1) - p(n-2, 1) + 4./dt*p(n-1, 0);
    }

    Matrix posVel = coeff.at(idBody).at(idSub) * p;
    if(planet == MOON)
    {
      posVel *= earthMoonRatio/(1.+earthMoonRatio);     // relative to earth-moon bary center
      posVel += interpolate(time, EARTHMOONBARYCENTER); // add earth-moon bary center
    }
    if(planet == EARTH)
    {
      posVel *= -1./(1.+earthMoonRatio);                // relative to earth-moon bary center
      posVel += interpolate(time, EARTHMOONBARYCENTER); // add earth-moon bary center
    }
    return posVel;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW_EXTRA("filename=<"+file.fileName().str()+">", e)
  }
}

/***********************************************/

void InFileEphemerides::ephemeris(const Time &time, Planet planet, Planet center, Vector3d &position, Vector3d &velocity)
{
  try
  {
    Matrix posVel;
    if(planet == MOON && center == EARTH)
      posVel = interpolate(time, MOONFROMEARTH);
    else
      posVel = interpolate(time, planet) - interpolate(time, center);
    position = Vector3d(posVel.column(0));
    velocity = Vector3d(posVel.column(1));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW_EXTRA("filename=<"+file.fileName().str()+">", e)
  }
}


/***********************************************/
/***********************************************/

void writeFileEphemerides(const FileName &fileName, Double earthMoonRatio, const std::vector<Time> &times,
                          const std::vector<UInt> &subintervals, const std::vector<UInt> &components, const std::vector<UInt> &degree,
                          const std::vector<std::vector<std::vector<Matrix>>> &coefficients)
{
  try
  {
    OutFileArchive file(fileName, FILE_EPHEMERIDES_TYPE);
    file<<nameValue("bodyCount",      components.size());
    file<<nameValue("intervalsCount", times.size()-1);
    file<<nameValue("earthMoonRatio", earthMoonRatio);
    file.comment("subintervals  components   degree");
    file.comment("=================================");
    for(UInt idBody=0; idBody<components.size(); idBody++)
    {
      file<<beginGroup("body");
      file<<nameValue("subintervals", subintervals.at(idBody));
      file<<nameValue("components",   components.at(idBody));
      file<<nameValue("degree",       degree.at(idBody));
      file<<endGroup("body");
    }
    file.comment("times");
    file.comment("=====");
    for(UInt idEpoch=0; idEpoch<times.size(); idEpoch++)
      file<<nameValue("time", times.at(idEpoch));

    file.comment("Chebyshev polynomial coefficients");
    file.comment("=================================");
    for(UInt idInterval=0; idInterval<coefficients.size(); idInterval++)
    {
      file<<beginGroup("interval");
      for(UInt idBody=0; idBody<subintervals.size(); idBody++)
      {
        file<<beginGroup("body");
        for(UInt idSub=0; idSub<subintervals.at(idBody); idSub++)
        {
          file<<beginGroup("subinterval");
          for(UInt n=0; n<=degree.at(idBody); n++)
            for(UInt i=0; i<components.at(idBody); i++)
              file<<nameValue("coefficient", coefficients.at(idInterval).at(idBody).at(idSub)(i, n));
          file<<endGroup("subinterval");
        }
        file<<endGroup("body");
      }
      file<<endGroup("interval");
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW_EXTRA("filename=<"+fileName.str()+">", e)
  }
}

/***********************************************/
