/***********************************************/
/**
* @file parametrizationSatelliteTrackingScale.h
*
* @brief SST scale.
*
* @author Torsten Mayer-Guerr
* @date 2017-08-03
*
*/
/***********************************************/

#ifndef __GROOPS_PARAMETRIZATIONSATELLITETRACKINGSCALE__
#define __GROOPS_PARAMETRIZATIONSATELLITETRACKINGSCALE__

// Latex documentation
#ifdef DOCSTRING_ParametrizationSatelliteTracking
static const char *docstringParametrizationSatelliteTrackingScale = R"(
\subsection{Scale}\label{parametrizationSatelliteTrackingType:scale}
Estimate scale factor for SST observations with respect to reference SST data
\configFile{inputfileSatelliteTracking}{instrument}.
The temporal variation is defined by \configClass{parametrizationTemporal}{parametrizationTemporalType}.

The \file{parameter names}{parameterName} are \verb|satellite1.satellite2:sstScale:<temporal>:<interval>|.
)";
#endif

/***********************************************/

#include "base/import.h"
#include "files/fileInstrument.h"
#include "classes/parametrizationTemporal/parametrizationTemporal.h"
#include "parametrizationSatelliteTracking.h"

/***** CLASS ***********************************/

/** @brief SST scale.
* @ingroup parametrizationSatelliteTrackingGroup
* @see ParametrizationSatelliteTracking */
class ParametrizationSatelliteTrackingScale : public ParametrizationSatelliteTrackingBase
{
  SatelliteTrackingArc       sst;
  ParametrizationTemporalPtr temporal;
  Bool                       perArc;
  UInt                       idx;

public:
  ParametrizationSatelliteTrackingScale(Config &config);

  Bool isPerArc() const {return perArc;}
  Bool setInterval(const Time &timeStart, const Time &timeEnd);
  UInt parameterCount() const {return temporal->parameterCount();}
  void parameterName(std::vector<ParameterName> &name) const;
  void compute(UInt sstType, const std::vector<Time> &time, const Vector &sst0,
               const Vector &position1, const Vector &position2, const Vector &velocity1, const Vector &velocity2,
               const std::vector<Rotary3d> &rotSat1, const std::vector<Rotary3d> &rotSat2, MatrixSliceRef A);
};

/***********************************************/
/***** Inlines *********************************/
/***********************************************/

inline ParametrizationSatelliteTrackingScale::ParametrizationSatelliteTrackingScale(Config &config)
{
  try
  {
    FileName fileNameSst;

    readConfig(config, "inputfileSatelliteTracking", fileNameSst, Config::MUSTSET,  "",  "");
    readConfig(config, "temporal",                   temporal,    Config::MUSTSET,  "",  "");
    readConfig(config, "perArc",                     perArc,      Config::DEFAULT,  "0", "");
    if(isCreateSchema(config)) return;

    sst = InstrumentFile::read(fileNameSst);
    idx = 0;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Bool ParametrizationSatelliteTrackingScale::setInterval(const Time &timeStart, const Time &timeEnd)
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

inline void ParametrizationSatelliteTrackingScale::parameterName(std::vector<ParameterName> &name) const
{
  try
  {
    std::vector<ParameterName> baseName(1);
    baseName.at(0) = ParameterName("satellite1.satellite2", "sstScale");
    temporal->parameterName(baseName, name);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline void ParametrizationSatelliteTrackingScale::compute(UInt sstType, const std::vector<Time> &time, const Vector &/*sst0*/,
                                                     const Vector &/*position1*/, const Vector &/*position2*/, const Vector &/*velocity1*/, const Vector &/*velocity2*/,
                                                     const std::vector<Rotary3d> &/*rotSat1*/, const std::vector<Rotary3d> &/*rotSat2*/, MatrixSliceRef A)
{
  try
  {
    if((time.at(0) < sst.at(0).time) || (time.back() > sst.at(sst.size()-1).time))
      throw(Exception("time not given in satellite tracking file: "+time.at(0).dateTimeStr()));

    for(UInt idEpoch=0; idEpoch<time.size(); idEpoch++)
    {
      // find index
      if((idx >= sst.size()) || (time.at(idEpoch) < sst.at(idx).time))
        idx = 0;
      while(time.at(idEpoch) > sst.at(idx).time)
        idx++;
      if(time.at(idEpoch) != sst.at(idx).time)
        throw(Exception("time not given in sst file: "+time.at(idEpoch).dateTimeStr()));
      Double scale = 0;
           if(sstType == 0) scale = sst.at(idx).range;
      else if(sstType == 1) scale = sst.at(idx).rangeRate;
      else if(sstType == 2) scale = sst.at(idx).rangeAcceleration;

      std::vector<UInt>   index;
      std::vector<Double> factor;
      temporal->factors(time.at(idEpoch), index, factor);
      for(UInt i=0; i<index.size(); i++)
        A(idEpoch, index.at(i)) = scale * factor.at(i);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

#endif
