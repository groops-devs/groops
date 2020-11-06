/***********************************************/
/**
* @file parametrizationSatelliteTrackingBias.h
*
* @brief SST bias.
*
* @author Beate Klinger
* @date 2015-05-12
*
*/
/***********************************************/

#ifndef __GROOPS_PARAMETRIZATIONSATELLITETRACKINGBIAS__
#define __GROOPS_PARAMETRIZATIONSATELLITETRACKINGBIAS__

// Latex documentation
#ifdef DOCSTRING_ParametrizationSatelliteTracking
static const char *docstringParametrizationSatelliteTrackingBias = R"(
\subsection{Bias}\label{parametrizationSatelliteTrackingType:bias}
Estimate bias for SST observations. The temporal variation is defined by \configClass{parametrizationTemporal}{parametrizationTemporalType}.
)";
#endif

/***********************************************/

#include "base/import.h"
#include "classes/parametrizationTemporal/parametrizationTemporal.h"
#include "parametrizationSatelliteTracking.h"

/***** CLASS ***********************************/

/** @brief SST bias.
* @ingroup parametrizationSatelliteTrackingGroup
* @see ParametrizationSatelliteTracking */
class ParametrizationSatelliteTrackingBias : public ParametrizationSatelliteTrackingBase
{
  ParametrizationTemporalPtr temporal;
  Bool                       perArc;

public:
  ParametrizationSatelliteTrackingBias(Config &config);

  Bool isPerArc() const {return perArc;}
  void setInterval(const Time &timeStart, const Time &timeEnd);
  UInt parameterCount() const {return temporal->parameterCount();}
  void parameterName(std::vector<ParameterName> &name) const;
  void compute(UInt sstType, const std::vector<Time> &time, const Vector &sst0,
               const Vector &position1, const Vector &position2, const Vector &velocity1, const Vector &velocity2,
               const std::vector<Rotary3d> &rotSat1, const std::vector<Rotary3d> &rotSat2, MatrixSliceRef A);
};

/***********************************************/
/***** Inlines *********************************/
/***********************************************/

inline ParametrizationSatelliteTrackingBias::ParametrizationSatelliteTrackingBias(Config &config)
{
  try
  {
    readConfig(config, "temporal", temporal, Config::MUSTSET,  "",  "");
    readConfig(config, "perArc",   perArc,   Config::DEFAULT,  "0", "");
    if(isCreateSchema(config)) return;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void ParametrizationSatelliteTrackingBias::setInterval(const Time &timeStart, const Time &timeEnd)
{
  try
  {
    temporal->setInterval(timeStart, timeEnd, perArc);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}


/***********************************************/

inline void ParametrizationSatelliteTrackingBias::parameterName(std::vector<ParameterName> &name) const
{
  try
  {
    std::vector<ParameterName> baseName(1);
    baseName.at(0) = ParameterName("satellite1.satellite2", "sstBias");
    temporal->parameterName(baseName, name);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline void ParametrizationSatelliteTrackingBias::compute(UInt /*sstType*/, const std::vector<Time> &time, const Vector &/*sst0*/,
                                                    const Vector &/*position1*/, const Vector &/*position2*/, const Vector &/*velocity1*/, const Vector &/*velocity2*/,
                                                    const std::vector<Rotary3d> &/*rotSat1*/, const std::vector<Rotary3d> &/*rotSat2*/, MatrixSliceRef A)
{
  try
  {
    for(UInt idEpoch=0; idEpoch<time.size(); idEpoch++)
    {
      std::vector<UInt>   index;
      std::vector<Double> factor;
      temporal->factors(time.at(idEpoch), index, factor);
      for(UInt i=0; i<index.size(); i++)
        A(idEpoch, index.at(i)) = factor.at(i);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

#endif
