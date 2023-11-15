/***********************************************/
/**
* @file parametrizationAccelerationAccBias.h
*
* @brief Accelerometer bias (time variable).
*
* @author Torsten Mayer-Guerr
* @date 2015-05-31
*
*/
/***********************************************/

#ifndef __GROOPS_PARAMETRIZATIONACCELERATIONACCBIAS__
#define __GROOPS_PARAMETRIZATIONACCELERATIONACCBIAS__

// Latex documentation
#ifdef DOCSTRING_ParametrizationAcceleration
static const char *docstringParametrizationAccelerationAccBias = R"(
\subsection{AccBias}\label{parametrizationAccelerationType:accBias}
Temporal changing accelerometer bias per axis in $[m/s^2]$ in Satellite Reference Frame (SRF).
If the attitude of the satellite is not provided the Celestial Reference Frame (CRF) is used instead.

The \file{parameter names}{parameterName} are
\begin{itemize}
\item \verb|*:accBias.x:*:*|,
\item \verb|*:accBias.y:*:*|,
\item \verb|*:accBias.z:*:*|.
\end{itemize}
)";
#endif

/***********************************************/

#include "base/import.h"
#include "classes/parametrizationTemporal/parametrizationTemporal.h"
#include "parametrizationAcceleration.h"

/***** CLASS ***********************************/

/** @brief Accelerometer bias (time variable).
* @ingroup parametrizationAccelerationGroup
* @see ParametrizationAcceleration */
class ParametrizationAccelerationAccBias : public ParametrizationAccelerationBase
{
  ParametrizationTemporalPtr temporal;
  UInt countAxis;
  Bool estimateX, estimateY, estimateZ;
  Bool perArc;

public:
  ParametrizationAccelerationAccBias(Config &config);

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

inline ParametrizationAccelerationAccBias::ParametrizationAccelerationAccBias(Config &config)
{
  try
  {
    readConfig(config, "estimateX", estimateX,  Config::DEFAULT,  "1", "along");
    readConfig(config, "estimateY", estimateY,  Config::DEFAULT,  "1", "cross");
    readConfig(config, "estimateZ", estimateZ,  Config::DEFAULT,  "1", "radial");
    readConfig(config, "temporal",  temporal,   Config::MUSTSET,  "",  "");
    readConfig(config, "perArc",    perArc,     Config::DEFAULT,  "0", "");
    if(isCreateSchema(config)) return;

    countAxis = estimateX+estimateY+estimateZ;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline void ParametrizationAccelerationAccBias::parameterName(std::vector<ParameterName> &name) const
{
  try
  {
    std::vector<ParameterName> baseName;
    if(estimateX) baseName.push_back(ParameterName("satellite", "accBias.x"));
    if(estimateY) baseName.push_back(ParameterName("satellite", "accBias.y"));
    if(estimateZ) baseName.push_back(ParameterName("satellite", "accBias.z"));
    temporal->parameterName(baseName, name);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline void ParametrizationAccelerationAccBias::compute(SatelliteModelPtr /*satellite*/, const Time &time, const Vector3d &/*position*/, const Vector3d &/*velocity*/,
                                                        const Rotary3d &rotSat, const Rotary3d &rotEarth, EphemeridesPtr /*ephemerides*/, MatrixSliceRef A)
{
  try
  {
    const Matrix rotary = (rotEarth*rotSat).matrix();
    Matrix R(3, countAxis);
    UInt idxAxis = 0;
    if(estimateX) copy(rotary.column(0), R.column(idxAxis++));
    if(estimateY) copy(rotary.column(1), R.column(idxAxis++));
    if(estimateZ) copy(rotary.column(2), R.column(idxAxis++));

    temporal->designMatrix(time, R, A);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}


/***********************************************/

#endif
