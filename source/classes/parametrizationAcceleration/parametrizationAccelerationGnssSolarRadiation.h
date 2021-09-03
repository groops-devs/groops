/***********************************************/
/**
* @file parametrizationAccelerationGnssSolarRadiation.h
*
* @brief GNSS solar radtion pressure model.
*
* @author Sebastian Strasser
* @date 2013-12-18
*
*/
/***********************************************/

#ifndef __GROOPS_PARAMETRIZATIONACCELERATIONGNSSSOLARRADIATION__
#define __GROOPS_PARAMETRIZATIONACCELERATIONGNSSSOLARRADIATION__

// Latex documentation
#ifdef DOCSTRING_ParametrizationAcceleration
static const char *docstringParametrizationAccelerationGnssSolarRadiation = R"(
\subsection{GnssSolarRadiation}\label{parametrizationAccelerationType:gnssSolarRadiation}
GNSS solar radiation pressure model.
)";
#endif


/***********************************************/

#include "base/import.h"
#include "parametrizationAcceleration.h"
#include "inputOutput/logging.h"
#include "classes/eclipse/eclipse.h"

/***** CLASS ***********************************/

/** @brief Oscillation per revoultion.
* @ingroup parametrizationAccelerationGroup
* @see ParametrizationAcceleration */
class ParametrizationAccelerationGnssSolarRadiation : public ParametrizationAccelerationBase
{
  UInt       countParameter;
  Bool       d0, d2, d4;
  Bool       y0;
  Bool       b0, b1, b3;
  Bool       perArc;
  EclipsePtr eclipse;

public:
  ParametrizationAccelerationGnssSolarRadiation(Config &config);

  Bool isPerArc() const override {return perArc;}
  Bool setInterval(const Time &/*timeStart*/, const Time &/*timeEnd*/) override {return FALSE;}
  UInt parameterCount() const override {return countParameter;}
  void parameterName(std::vector<ParameterName> &name) const override;

  void compute(SatelliteModelPtr satellite, const Time &time, const Vector3d &position, const Vector3d &velocity,
               const Rotary3d &rotSat, const Rotary3d &rotEarth, EphemeridesPtr ephemerides, MatrixSliceRef A) override;
};

/***********************************************/

inline ParametrizationAccelerationGnssSolarRadiation::ParametrizationAccelerationGnssSolarRadiation(Config &config)
{
  try
  {
    d0 = d2 = d4 = FALSE;
    y0 = FALSE;
    b0 = b1 = b3 = FALSE;

    readConfig(config, "estimateD0", d0,      Config::DEFAULT,  "1", "constant term along D-axis (sat-sun vector)");
    readConfig(config, "estimateD2", d2,      Config::DEFAULT,  "1", "2-per-rev terms along D-axis");
    readConfig(config, "estimateD4", d4,      Config::DEFAULT,  "0", "4-per-rev terms along D-axis");
    readConfig(config, "estimateY0", y0,      Config::DEFAULT,  "1", "constant term along Y-axis (solar panel axis)");
    readConfig(config, "estimateB0", b0,      Config::DEFAULT,  "1", "constant term along B-axis (cross product D x Y)");
    readConfig(config, "estimateB1", b1,      Config::DEFAULT,  "1", "1-per-rev terms along B-axis");
    readConfig(config, "estimateB3", b3,      Config::DEFAULT,  "0", "3-per-rev terms along B-axis");
    readConfig(config, "perArc",     perArc,  Config::DEFAULT,  "0", "");
    readConfig(config, "eclipse",    eclipse, Config::MUSTSET,  "",  "");
    if(isCreateSchema(config)) return;

    countParameter = 0;
    if(d0) countParameter++;
    if(d2) countParameter += 2;
    if(d4) countParameter += 2;

    if(y0) countParameter++;

    if(b0) countParameter++;
    if(b1) countParameter += 2;
    if(b3) countParameter += 2;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline void ParametrizationAccelerationGnssSolarRadiation::parameterName(std::vector<ParameterName> &name) const
{
  if(d0) name.push_back(ParameterName("satellite", "solarRadiationPressure.ECOM.D0"));
  if(d2) name.push_back(ParameterName("satellite", "solarRadiationPressure.ECOM.DC2"));
  if(d2) name.push_back(ParameterName("satellite", "solarRadiationPressure.ECOM.DS2"));
  if(d4) name.push_back(ParameterName("satellite", "solarRadiationPressure.ECOM.DC4"));
  if(d4) name.push_back(ParameterName("satellite", "solarRadiationPressure.ECOM.DS4"));

  if(y0) name.push_back(ParameterName("satellite", "solarRadiationPressure.ECOM.Y0"));

  if(b0) name.push_back(ParameterName("satellite", "solarRadiationPressure.ECOM.B0"));
  if(b1) name.push_back(ParameterName("satellite", "solarRadiationPressure.ECOM.BC1"));
  if(b1) name.push_back(ParameterName("satellite", "solarRadiationPressure.ECOM.BS1"));
  if(b3) name.push_back(ParameterName("satellite", "solarRadiationPressure.ECOM.BC3"));
  if(b3) name.push_back(ParameterName("satellite", "solarRadiationPressure.ECOM.BS3"));
}

/***********************************************/

inline void ParametrizationAccelerationGnssSolarRadiation::compute(SatelliteModelPtr /*satellite*/, const Time &time, const Vector3d &position, const Vector3d &velocity,
                                                          const Rotary3d &rotSat, const Rotary3d &rotEarth, EphemeridesPtr ephemerides, MatrixSliceRef A)
{
  try
  {
    if(!ephemerides)
      throw(Exception("No ephemerides given"));

    // Positions in TRF
    const Vector3d posSat = rotEarth.rotate(position);
    const Vector3d posSun = rotEarth.rotate(ephemerides->position(time, Ephemerides::SUN));

    // Argument of latitude of satellite relative to sun
    Vector3d z = normalize(crossProduct(posSat, rotEarth.rotate(velocity))); // Orbital plane normal vector
    Vector3d y = normalize(crossProduct(z, posSun));
    Vector3d x = normalize(crossProduct(y, z));                              // Vector sat-sun projected into orbital plane
    const Double du = std::atan2(inner(posSat, y), inner(posSat, x));

    // Satellite DYB coordinate system (in TRF)
    Vector3d dSat = normalize(posSun-posSat);
    Vector3d ySat = rotEarth.rotate(rotSat.rotate(Vector3d(0,1,0)));
    Vector3d bSat = normalize(crossProduct(dSat, ySat));

    // Design matrix elements
    const Double factor = 1e-9 * (eclipse ? eclipse->factor(time, position, ephemerides) : 1);  // Estimate accelerations in [nm/s^2]
    UInt idx = 0;

    if(d0) axpy(factor,                dSat.vector(), A.column(idx++)); // Constant D0  term
    if(d2) axpy(factor*std::cos(2*du), dSat.vector(), A.column(idx++)); // Periodic DC2 term
    if(d2) axpy(factor*std::sin(2*du), dSat.vector(), A.column(idx++)); // Periodic DS2 term
    if(d4) axpy(factor*std::cos(4*du), dSat.vector(), A.column(idx++)); // Periodic DC4 term
    if(d4) axpy(factor*std::sin(4*du), dSat.vector(), A.column(idx++)); // Periodic DS4 term

    if(y0) axpy(factor,                ySat.vector(), A.column(idx++)); // Constant Y0  term

    if(b0) axpy(factor,                bSat.vector(), A.column(idx++)); // Constant B0  term
    if(b1) axpy(factor*std::cos(du),   bSat.vector(), A.column(idx++)); // Periodic BC1 term
    if(b1) axpy(factor*std::sin(du),   bSat.vector(), A.column(idx++)); // Periodic BS1 term
    if(b3) axpy(factor*std::cos(3*du), bSat.vector(), A.column(idx++)); // Periodic BC3 term
    if(b3) axpy(factor*std::sin(3*du), bSat.vector(), A.column(idx++)); // Periodic BS3 term
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

#endif
