/***********************************************/
/**
* @file orbitPropagatorFile.h
*
* @brief Read orbit from file.
* @see orbitPropagator
*
* @author Matthias Ellmer
* @date 2017-01-25
*
*/
/***********************************************/

#ifndef __GROOPS_ORBITFILE__
#define __GROOPS_ORBITFILE__

// Latex documentation
#ifdef DOCSTRING_OrbitPropagator
static const char *docstringOrbitPropagatorFile = R"(
\subsection{File}
Reads an orbit from file. If the needed epochs are not given an exception is thrown.
)";
#endif

/***********************************************/

#include "classes/orbitPropagator/orbitPropagator.h"

/***** CLASS ***********************************/

/** @brief Read orbit from file.
* @ingroup orbitPropagatorGroup
* @see orbitPropagator */
class OrbitPropagatorFile : public OrbitPropagator
{
  FileName fileNameOrbit;
  Double   margin;
  Bool     computeForces;

public:
  OrbitPropagatorFile(Config &config);

  OrbitArc integrateArc(OrbitEpoch startEpoch, Time sampling, UInt posCount, ForcesPtr forces, SatelliteModelPtr satellite,
                        EarthRotationPtr earthRotation, EphemeridesPtr ephemerides, Bool timing) const override;
};

/***********************************************/

inline OrbitPropagatorFile::OrbitPropagatorFile(Config &config)
{
  try
  {
    readConfig(config, "inputfileOrbit",  fileNameOrbit, Config::MUSTSET,  "",     "epoch at timeStart is not used");
    readConfig(config, "margin",          margin,        Config::MUSTSET,  "1e-5", "[seconds] to find identical times");
    readConfig(config, "recomputeForces", computeForces, Config::DEFAULT,  "0",    "");
    if(isCreateSchema(config)) return;
  }
  catch (std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline OrbitArc OrbitPropagatorFile::integrateArc(OrbitEpoch startEpoch, Time sampling, UInt posCount, ForcesPtr forces,
                                                  SatelliteModelPtr satellite, EarthRotationPtr earthRotation, EphemeridesPtr ephemerides, Bool timing) const
{
  try
  {
    OrbitArc orbitFile = InstrumentFile::read(fileNameOrbit);

    OrbitArc orbit;
    if(timing) logTimerStart;
    for(UInt k=0; k<posCount; k++)
    {
      if(timing) logTimerLoop(k, posCount);

      const Time time = startEpoch.time + k*sampling;
      OrbitEpoch epoch;
      for(UInt i=0; i<orbitFile.size(); i++)
        if(std::fabs((orbitFile.at(i).time-time).seconds()) < margin)
          epoch = orbitFile.at(i);
      if(std::fabs((epoch.time-time).seconds())>margin)
        throw (Exception("Requested time "+time.dateTimeStr()+" not present in file <"+fileNameOrbit.str()+">"));

      if(computeForces)
        epoch.acceleration = acceleration(epoch, forces, satellite, earthRotation, ephemerides);

      orbit.push_back(epoch);
    }
    if(timing) logTimerLoopEnd(posCount);

    return orbit;
  }
  catch (std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

#endif
