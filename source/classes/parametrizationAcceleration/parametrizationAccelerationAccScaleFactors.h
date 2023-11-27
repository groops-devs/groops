/***********************************************/
/**
* @file parametrizationAccelerationAccScaleFactors.h
*
* @brief Accelerometer scale factors.
*
* @author Torsten Mayer-Guerr
* @author Beate Klinger
* @date 2015-06-01
*
*/
/***********************************************/

#ifndef __GROOPS_PARAMETRIZATIONACCELERATIONACCSCALEFACTORS__
#define __GROOPS_PARAMETRIZATIONACCELERATIONACCSCALEFACTORS__

// Latex documentation
#ifdef DOCSTRING_ParametrizationAcceleration
static const char *docstringParametrizationAccelerationAccScaleFactors = R"(
\subsection{AccelerometerScaleFactors}\label{parametrizationAccelerationType:accelerometerScaleFactors}
Accelerometer scale factor per axis.

The \file{parameter names}{parameterName} are
\begin{itemize}
\item \verb|*:accScale.x:<temporal>:<interval>|,
\item \verb|*:accScale.y:<temporal>:<interval>|,
\item \verb|*:accScale.z:<temporal>:<interval>|,
\item \verb|*:accScaleCross.xy:<temporal>:<interval>|,
\item \verb|*:accScaleCross.xz:<temporal>:<interval>|,
\item \verb|*:accScaleCross.yz:<temporal>:<interval>|,
\item \verb|*:accScaleRotation.xy:<temporal>:<interval>|,
\item \verb|*:accScaleRotation.xz:<temporal>:<interval>|,
\item \verb|*:accScaleRotation.yz:<temporal>:<interval>|.
\end{itemize}

This parametrization needs the attitude of the satellite.
)";
#endif

/***********************************************/

#include "base/import.h"
#include "files/fileInstrument.h"
#include "classes/parametrizationTemporal/parametrizationTemporal.h"
#include "parametrizationAcceleration.h"

/***** CLASS ***********************************/

/** @brief Accelerometer scale factors.
* @ingroup parametrizationAccelerationGroup
* @see ParametrizationAcceleration */
class ParametrizationAccelerationAccScaleFactors : public ParametrizationAccelerationBase
{
  AccelerometerArc           accelerometer;
  ParametrizationTemporalPtr temporal;
  UInt                       countAxis;
  Bool                       estimateX, estimateY, estimateZ;
  Bool                       estimateCross, estimateRotation;
  Bool                       perArc;
  UInt                       idx;

public:
  ParametrizationAccelerationAccScaleFactors(Config &config);

  Bool isPerArc() const override {return perArc;}
  Bool setInterval(const Time &timeStart, const Time &timeEnd) override {return temporal->setInterval(timeStart, timeEnd, perArc);}
  UInt parameterCount() const override {return countAxis*temporal->parameterCount();}
  void parameterName(std::vector<ParameterName> &name) const override;
  void compute(SatelliteModelPtr satellite, const Time &time, const Vector3d &position, const Vector3d &velocity,
               const Rotary3d &rotSat, const Rotary3d &rotEarth, EphemeridesPtr ephemerides, MatrixSliceRef A) override;
};

/***********************************************/
/***** Inlines *********************************/
/***********************************************/

inline ParametrizationAccelerationAccScaleFactors::ParametrizationAccelerationAccScaleFactors(Config &config)
{
  try
  {
    FileName fileNameAcc;

    readConfig(config, "inputfileAccelerometer", fileNameAcc, Config::MUSTSET,  "",  "");
    readConfig(config, "estimateX",         estimateX,        Config::DEFAULT,  "1", "along");
    readConfig(config, "estimateY",         estimateY,        Config::DEFAULT,  "1", "cross");
    readConfig(config, "estimateZ",         estimateZ,        Config::DEFAULT,  "1", "radial");
    readConfig(config, "estimateCrossTalk", estimateCross,    Config::DEFAULT,  "0", "non-orthognality of axes");
    readConfig(config, "estimateRotation",  estimateRotation, Config::DEFAULT,  "0", "misalignment");
    readConfig(config, "temporal",          temporal,         Config::MUSTSET,  "",  "");
    readConfig(config, "perArc",            perArc,           Config::DEFAULT,  "0", "");
    if(isCreateSchema(config)) return;

    accelerometer = InstrumentFile::read(fileNameAcc);
    countAxis     = estimateX+estimateY+estimateZ;
    if(estimateCross)    countAxis += 3;
    if(estimateRotation) countAxis += 3;
    idx = 0;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline void ParametrizationAccelerationAccScaleFactors::parameterName(std::vector<ParameterName> &name) const
{
  try
  {
    std::vector<ParameterName> baseName;
    if(estimateX) baseName.push_back(ParameterName("satellite", "accScale.x"));
    if(estimateY) baseName.push_back(ParameterName("satellite", "accScale.y"));
    if(estimateZ) baseName.push_back(ParameterName("satellite", "accScale.z"));
    if(estimateCross)
    {
      baseName.push_back(ParameterName("satellite", "accScaleCross.xy"));
      baseName.push_back(ParameterName("satellite", "accScaleCross.xz"));
      baseName.push_back(ParameterName("satellite", "accScaleCross.yz"));
    }
    if(estimateRotation)
    {
      baseName.push_back(ParameterName("satellite", "accScaleRotation.xy"));
      baseName.push_back(ParameterName("satellite", "accScaleRotation.xz"));
      baseName.push_back(ParameterName("satellite", "accScaleRotation.yz"));
    }
    temporal->parameterName(baseName, name);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline void ParametrizationAccelerationAccScaleFactors::compute(SatelliteModelPtr /*satellite*/, const Time &time, const Vector3d &/*position*/, const Vector3d &/*velocity*/,
                                                       const Rotary3d &rotSat, const Rotary3d &rotEarth, EphemeridesPtr /*ephemerides*/, MatrixSliceRef A)
{
  try
  {
    if((time<accelerometer.at(0).time)||(time>accelerometer.at(accelerometer.size()-1).time))
      throw(Exception("time not given in accelerometer file: "+time.dateTimeStr()));

    // find index (interpolations interval)
    if((idx>=accelerometer.size()) || (time<accelerometer.at(idx).time))
      idx = 0;
    while(time>accelerometer.at(idx).time)
      idx++;
    if(time!=accelerometer.at(idx).time)
      throw(Exception("time not given in accelerometer file: "+time.dateTimeStr()));

    const Matrix rotary = (rotEarth*rotSat).matrix();
    Matrix R(3, countAxis);
    UInt idxAxis = 0;
    if(estimateX) axpy(accelerometer.at(idx).acceleration.x(), rotary.column(0), R.column(idxAxis++));
    if(estimateY) axpy(accelerometer.at(idx).acceleration.y(), rotary.column(1), R.column(idxAxis++));
    if(estimateZ) axpy(accelerometer.at(idx).acceleration.z(), rotary.column(2), R.column(idxAxis++));

    if(estimateCross)
    {
      axpy(accelerometer.at(idx).acceleration.y(), rotary.column(0), R.column(idxAxis));
      axpy(accelerometer.at(idx).acceleration.x(), rotary.column(1), R.column(idxAxis++));
      axpy(accelerometer.at(idx).acceleration.z(), rotary.column(0), R.column(idxAxis));
      axpy(accelerometer.at(idx).acceleration.x(), rotary.column(2), R.column(idxAxis++));
      axpy(accelerometer.at(idx).acceleration.z(), rotary.column(1), R.column(idxAxis));
      axpy(accelerometer.at(idx).acceleration.y(), rotary.column(2), R.column(idxAxis++));
    }

    if(estimateRotation)
    {
      axpy(-accelerometer.at(idx).acceleration.y(), rotary.column(0), R.column(idxAxis));
      axpy( accelerometer.at(idx).acceleration.x(), rotary.column(1), R.column(idxAxis++));
      axpy( accelerometer.at(idx).acceleration.z(), rotary.column(0), R.column(idxAxis));
      axpy(-accelerometer.at(idx).acceleration.x(), rotary.column(2), R.column(idxAxis++));
      axpy(-accelerometer.at(idx).acceleration.z(), rotary.column(1), R.column(idxAxis));
      axpy( accelerometer.at(idx).acceleration.y(), rotary.column(2), R.column(idxAxis++));
    }

    temporal->designMatrix(time, R, A);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

#endif
