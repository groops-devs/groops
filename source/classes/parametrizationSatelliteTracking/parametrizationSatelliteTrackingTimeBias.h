/***********************************************/
/**
* @file parametrizationSatelliteTrackingTimeBias.h
*
* @brief SST time bias (temporal changing).
*
* @author Torsten Mayer-Guerr
* @date 2015-06-01
*
*/
/***********************************************/

#ifndef __GROOPS_PARAMETRIZATIONSATELLITETRACKINGTIMEBIAS__
#define __GROOPS_PARAMETRIZATIONSATELLITETRACKINGTIMEBIAS__

// Latex documentation
#ifdef DOCSTRING_ParametrizationSatelliteTracking
static const char *docstringParametrizationSatelliteTrackingTimeBias = R"(
\subsection{TimeBias}\label{parametrizationSatelliteTrackingType:timeBias}
Estimate time shift in seconds in SST observations, with defined temporal variation by \configClass{parametrizationTemporal}{parametrizationTemporalType}. The design matrix is computed by taking the derivative of the ranging data w.r.t. time.
)";
#endif

/***********************************************/

#include "base/import.h"
#include "base/polynomial.h"
#include "classes/parametrizationTemporal/parametrizationTemporal.h"
#include "parametrizationSatelliteTracking.h"

/***** CLASS ***********************************/

/** @brief SST time bias.
* @ingroup parametrizationSatelliteTrackingGroup
* @see ParametrizationSatelliteTracking */
class ParametrizationSatelliteTrackingTimeBias : public ParametrizationSatelliteTrackingBase
{
  UInt                       degree;
  ParametrizationTemporalPtr temporal;
  Bool                       perArc;

public:
  ParametrizationSatelliteTrackingTimeBias(Config &config);

  Bool isPerArc() const {return perArc;}
  Bool setInterval(const Time &timeStart, const Time &timeEnd);
  UInt parameterCount() const {return temporal->parameterCount();}
  void parameterName(std::vector<ParameterName> &name) const;
  void compute(UInt sstType, const std::vector<Time> &times, const Vector &sst0,
               const Vector &position1, const Vector &position2, const Vector &velocity1, const Vector &velocity2,
               const std::vector<Rotary3d> &rotSat1, const std::vector<Rotary3d> &rotSat2, MatrixSliceRef A);
};

/***********************************************/

inline ParametrizationSatelliteTrackingTimeBias::ParametrizationSatelliteTrackingTimeBias(Config &config)
{
  try
  {
    readConfig(config, "polynomialDegree",  degree,   Config::MUSTSET,  "10", "polynomial degree");
    readConfig(config, "temporal",          temporal, Config::MUSTSET,  "",  "");
    readConfig(config, "perArc",            perArc,   Config::DEFAULT,  "0", "");
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline Bool ParametrizationSatelliteTrackingTimeBias::setInterval(const Time &timeStart, const Time &timeEnd)
{
  try
  {
    return temporal->setInterval(timeStart, timeEnd, perArc);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline void ParametrizationSatelliteTrackingTimeBias::parameterName(std::vector<ParameterName> &name) const
{
  try
  {
    std::vector<ParameterName> baseName(1);
    baseName.at(0) = ParameterName("satellite1.satellite2", "sstTimeBias");
    temporal->parameterName(baseName, name);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline void ParametrizationSatelliteTrackingTimeBias::compute(UInt /*sstType*/, const std::vector<Time> &times, const Vector &sst0,
                                                              const Vector &/*position1*/, const Vector &/*position2*/, const Vector &/*velocity1*/, const Vector &/*velocity2*/,
                                                              const std::vector<Rotary3d> &/*rotSat1*/, const std::vector<Rotary3d> &/*rotSat2*/, MatrixSliceRef A)
{
  try
  {
    Polynomial polynomial(times, degree);
    const Vector rate = polynomial.derivative(times, sst0); // drho/dt
    for(UInt idEpoch=0; idEpoch<sst0.rows(); idEpoch++)
    {
      std::vector<UInt>   index;
      std::vector<Double> factor;
      temporal->factors(times.at(idEpoch), index, factor);
      for(UInt i=0; i<index.size(); i++)
        A(idEpoch, index.at(i)) = rate(idEpoch)*factor.at(i);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

#endif
