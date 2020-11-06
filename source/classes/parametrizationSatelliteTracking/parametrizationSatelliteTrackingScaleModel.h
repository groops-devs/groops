/***********************************************/
/**
* @file parametrizationSatelliteTrackingScaleModel.h
*
* @brief Estimate scale factors for model from file.
*
* @author Torsten Mayer-Guerr
* @date 2018-06-29
*
*/
/***********************************************/

#ifndef __GROOPS_ParametrizationSatelliteTrackingScaleModel__
#define __GROOPS_ParametrizationSatelliteTrackingScaleModel__

// Latex documentation
#ifdef DOCSTRING_ParametrizationSatelliteTracking
static const char *docstringParametrizationSatelliteTrackingScaleModel = R"(
\subsection{ScaleModel}\label{parametrizationSatelliteTrackingType:scaleModel}
Estimate scale factors for deterministic signal models from satellite tracking instrument file \configFile{inputfileSatelliteTracking}{instrument}, see \program{EnsembleAveragingScaleModel}.
Amplitude variation of model waveforms is defined by \configClass{parametrizationTemporal}{parametrizationTemporalType}.
)";
#endif

/***********************************************/

#include "base/import.h"
#include "base/polynomial.h"
#include "classes/parametrizationTemporal/parametrizationTemporal.h"
#include "parametrizationSatelliteTracking.h"

/***** CLASS ***********************************/

/** @brief Estimate scale factors for model from file.
* @ingroup parametrizationSatelliteTrackingGroup
* @see ParametrizationSatelliteTracking */
class ParametrizationSatelliteTrackingScaleModel : public ParametrizationSatelliteTrackingBase
{
  SatelliteTrackingArc       sst;
  ParametrizationTemporalPtr temporal;
  Bool                       perArc;
  UInt                       idx;

public:
  ParametrizationSatelliteTrackingScaleModel(Config &config);

  Bool isPerArc() const {return perArc;}
  void setInterval(const Time &timeStart, const Time &timeEnd);
  UInt parameterCount() const {return temporal->parameterCount();}
  void parameterName(std::vector<ParameterName> &name) const;
  void compute(UInt sstType, const std::vector<Time> &time, const Vector &sst0,
               const Vector &position1, const Vector &position2, const Vector &velocity1, const Vector &velocity2,
               const std::vector<Rotary3d> &rotSat1, const std::vector<Rotary3d> &rotSat2, MatrixSliceRef A);
};

/***********************************************/

inline ParametrizationSatelliteTrackingScaleModel::ParametrizationSatelliteTrackingScaleModel(Config &config)
{
  try
  {
    FileName fileName;

    readConfig(config, "inputfileSatelliteTracking", fileName, Config::MUSTSET,  "",  "");
    readConfig(config, "temporal",                   temporal, Config::MUSTSET,  "",  "");
    readConfig(config, "perArc",                     perArc,   Config::DEFAULT,  "0", "");
    if(isCreateSchema(config)) return;

    sst = InstrumentFile::read(fileName);
    idx = 0;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void ParametrizationSatelliteTrackingScaleModel::setInterval(const Time &timeStart, const Time &timeEnd)
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

inline void ParametrizationSatelliteTrackingScaleModel::parameterName(std::vector<ParameterName> &name) const
{
  try
  {
    std::vector<ParameterName> baseName(1);
    baseName.at(0) = ParameterName("satellite1.satellite2", "scaleModel");
    temporal->parameterName(baseName, name);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline void ParametrizationSatelliteTrackingScaleModel::compute(UInt sstType, const std::vector<Time> &time, const Vector &/*sst0*/,
                                                          const Vector &/*position1*/, const Vector &/*position2*/, const Vector &/*velocity1*/, const Vector &/*velocity2*/,
                                                          const std::vector<Rotary3d> &/*rotSat1*/, const std::vector<Rotary3d> &/*rotSat2*/, MatrixSliceRef A)
{
  try
  {
    for(UInt idEpoch=0; idEpoch<time.size(); idEpoch++)
    {
      if((time.at(idEpoch) < sst.at(0).time) || (time.at(idEpoch) > sst.at(sst.size()-1).time))
        throw(Exception("time not given in sst file: "+time.at(idEpoch).dateTimeStr()));

      // find index (interpolations interval)
      if((idx >= sst.size()) || (time.at(idEpoch) < sst.at(idx).time))
        idx = 0;
      while(time.at(idEpoch) > sst.at(idx).time)
        idx++;

      Double value = 0;
      switch(sstType)
      {
        case 0: value = sst.at(idx).range;             break;
        case 1: value = sst.at(idx).rangeRate;         break;
        case 2: value = sst.at(idx).rangeAcceleration; break;
      }

      std::vector<UInt>   index;
      std::vector<Double> factor;
      temporal->factors(time.at(idEpoch), index, factor);
      for(UInt i=0; i<index.size(); i++)
        A(idEpoch, index.at(i)) = value * factor.at(i);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

#endif
