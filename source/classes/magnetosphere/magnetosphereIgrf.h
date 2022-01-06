/***********************************************/
/**
* @file magnetosphereIGRF.h
*
* @brief Magentic field of the Earth.
*
* @author Torsten Mayer-Guerr
* @date 2020-06-08
*
*/
/***********************************************/

#ifndef __GROOPS_MAGNETOSPHEREIGRF__
#define __GROOPS_MAGNETOSPHEREIGRF__

// Latex documentation
#ifdef DOCSTRING_Magnetosphere
static const char *docstringMagnetosphereIgrf = R"(
\subsection{IGRF}
International Geomagnetic Reference Field.
)";
#endif

/***********************************************/

#include "base/polynomial.h"
#include "files/fileInstrument.h"
#include "external/igrf/igrf.h"
#include "classes/magnetosphere/magnetosphere.h"

/***** CLASS ***********************************/

/** @brief Magentic field of the Earth.
* @ingroup magnetosphereGroup
* @see Magnetosphere */
class MagnetosphereIgrf : public Magnetosphere
{
  Polynomial        polynomial;
  Matrix            lonlat;

public:
  MagnetosphereIgrf(Config &config);

  Vector3d geomagneticNorthPole(const Time &time) const override;
  Vector3d magenticFieldVector(const Time &time, const Vector3d &position) const override;
};

/***********************************************/

MagnetosphereIgrf::MagnetosphereIgrf(Config &config)
{
  try
  {
    FileName fileName;
    readConfig(config, "inputfileMagneticNorthPole", fileName, Config::OPTIONAL, "{groopsDataDir}/magnetosphere/magneticNorthPole.txt", "time series of north pole");
    if(isCreateSchema(config)) return;

    if(!fileName.empty())
    {
      MiscValuesArc arc = InstrumentFile::read(fileName);
      lonlat = arc.matrix().column(1, 2);
      polynomial.init(arc.times(), 1);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline Vector3d MagnetosphereIgrf::geomagneticNorthPole(const Time &time) const
{
  try
  {
    if(!lonlat.size())
      throw(Exception("magentic north pole data not provided"));
    const Matrix ll = polynomial.interpolate({time}, lonlat);
    return polar(Angle(ll(0,0)*DEG2RAD), Angle(ll(0,1)*DEG2RAD), 1);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline Vector3d MagnetosphereIgrf::magenticFieldVector(const Time &time, const Vector3d &position) const
{
  try
  {
    Double n,e,u,f;
    igrfSynthesis(0/*main-field*/, time.decimalYear(), 2/*geocentric*/, position.r()/1000, position.theta()*RAD2DEG, position.lambda()*RAD2DEG, n,e,u,f);
    return 1e-9*localNorthEastUp(position).transform(Vector3d(n,e,u));  // nT -> T (Tesla)
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

#endif
