/***********************************************/
/**
* @file parametrizationAccelerationModelScale.h
*
* @brief Model scale factor.
*
* @author Torsten Mayer-Guerr
* @date 2015-06-01
*
*/
/***********************************************/

#ifndef __GROOPS_PARAMETRIZATIONACCELERATIONMODELSCALE__
#define __GROOPS_PARAMETRIZATIONACCELERATIONMODELSCALE__

// Latex documentation
#ifdef DOCSTRING_ParametrizationAcceleration
static const char *docstringParametrizationAccelerationModelScale = R"(
\subsection{ModelScale}\label{parametrizationAccelerationType:modelScale}
Estimate a scale factor for a given model.
)";
#endif

/***********************************************/

#include "base/import.h"
#include "classes/miscAccelerations/miscAccelerations.h"
#include "classes/parametrizationTemporal/parametrizationTemporal.h"
#include "parametrizationAcceleration.h"

/***** CLASS ***********************************/

/** @brief Model scale factor.
* @ingroup parametrizationAccelerationGroup
* @see ParametrizationAcceleration */
class ParametrizationAccelerationModelScale : public ParametrizationAccelerationBase
{
  MiscAccelerationsPtr       model;
  ParametrizationTemporalPtr temporal;
  Bool                       perArc;

public:
  ParametrizationAccelerationModelScale(Config &config);

  Bool isPerArc() const override {return perArc;}
  Bool setInterval(const Time &timeStart, const Time &timeEnd) override {return temporal->setInterval(timeStart, timeEnd, perArc);}
  UInt parameterCount() const override {return temporal->parameterCount();}
  void parameterName(std::vector<ParameterName> &name) const override;
  void compute(SatelliteModelPtr satellite, const Time &time, const Vector3d &position, const Vector3d &velocity,
               const Rotary3d &rotSat, const Rotary3d &rotEarth, EphemeridesPtr ephemerides, MatrixSliceRef A) override;
};

/***********************************************/
/***** Inlines *********************************/
/***********************************************/

inline ParametrizationAccelerationModelScale::ParametrizationAccelerationModelScale(Config &config)
{
  try
  {
    readConfig(config, "miscAccelerations", model,    Config::MUSTSET,  "",  "");
    readConfig(config, "temporal",          temporal, Config::MUSTSET,  "",  "");
    readConfig(config, "perArc",            perArc,   Config::DEFAULT,  "0", "");
    if(isCreateSchema(config)) return;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline void ParametrizationAccelerationModelScale::parameterName(std::vector<ParameterName> &name) const
{
  try
  {
    std::vector<ParameterName> baseName;
    baseName.push_back(ParameterName("satellite", "modelScale"));
    temporal->parameterName(baseName, name);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline void ParametrizationAccelerationModelScale::compute(SatelliteModelPtr satellite, const Time &time, const Vector3d &position, const Vector3d &velocity,
                                                       const Rotary3d &rotSat, const Rotary3d &rotEarth, EphemeridesPtr ephemerides, MatrixSliceRef A)
{
  try
  {
    const Vector a = model->acceleration(satellite, time, position, velocity, rotSat, rotEarth, ephemerides).vector();
    temporal->designMatrix(time, a, A);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

#endif
