/***********************************************/
/**
* @file parametrizationSatelliteTrackingAntennaCenter.h
*
* @brief KBR antenna center parameter.
*
* @author Beate Klinger
* @date 2015-06-01
*
*/
/***********************************************/

#ifndef __GROOPS_PARAMETRIZATIONSATELLITETRACKINGANTENNACENTER__
#define __GROOPS_PARAMETRIZATIONSATELLITETRACKINGANTENNACENTER__

// Latex documentation
#ifdef DOCSTRING_ParametrizationSatelliteTracking
static const char *docstringParametrizationSatelliteTrackingAntennaCenter = R"(
\subsection{AntennaCenter}\label{parametrizationSatelliteTrackingType:antennaCenter}
Estimate the KBR antenna phase centre (APC) coordinates for each spacecraft in satellite reference frame (SRF)
as constant per axis, once per month. The observation equations are computed by taking the derivative
of the antenna offset correction equation w.r.t. the KBR APC coordinates.
)";
#endif

/***********************************************/

#include "base/import.h"
#include "base/polynomial.h"
#include "parametrizationSatelliteTracking.h"

/***** CLASS ***********************************/

/** @brief SST bias.
* @ingroup parametrizationSatelliteTrackingGroup
* @see ParametrizationSatelliteTracking */
class ParametrizationSatelliteTrackingAntennaCenter : public ParametrizationSatelliteTrackingBase
{
  Bool estimate1x, estimate1y, estimate1z;
  Bool estimate2x, estimate2y, estimate2z;
  UInt degree;

public:
  ParametrizationSatelliteTrackingAntennaCenter(Config &config);

  Bool isPerArc() const {return FALSE;}
  void setInterval(const Time &/*timeStart*/, const Time &/*timeEnd*/) {}
  UInt parameterCount() const {return estimate1x+estimate1y+estimate1z+estimate2x+estimate2y+estimate2z;}
  void parameterName(std::vector<ParameterName> &name) const;
  void compute(UInt sstType, const std::vector<Time> &time, const Vector &sst0,
               const Vector &position1, const Vector &position2, const Vector &velocity1, const Vector &velocity2,
               const std::vector<Rotary3d> &rotSat1, const std::vector<Rotary3d> &rotSat2, MatrixSliceRef A);
};

/***********************************************/
/***** Inlines *********************************/
/***********************************************/

inline ParametrizationSatelliteTrackingAntennaCenter::ParametrizationSatelliteTrackingAntennaCenter(Config &config)
{
  try
  {
    readConfig(config, "estimate1X",          estimate1x,  Config::DEFAULT,  "1", "along (satellite 1)");
    readConfig(config, "estimate1Y",          estimate1y,  Config::DEFAULT,  "1", "cross (satellite 1)");
    readConfig(config, "estimate1Z",          estimate1z,  Config::DEFAULT,  "1", "nadir (satellite 1)");
    readConfig(config, "estimate2X",          estimate2x,  Config::DEFAULT,  "1", "along (satellite 2)");
    readConfig(config, "estimate2Y",          estimate2y,  Config::DEFAULT,  "1", "cross (satellite 2)");
    readConfig(config, "estimate2Z",          estimate2z,  Config::DEFAULT,  "1", "nadir (satellite 2)");
    readConfig(config, "interpolationDegree", degree,      Config::DEFAULT,  "2", "differentiation by polynomial approximation of degree n");
    if(isCreateSchema(config)) return;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline void ParametrizationSatelliteTrackingAntennaCenter::parameterName(std::vector<ParameterName> &name) const
{
  try
  {
    if(estimate1x) name.push_back(ParameterName("satellite1", "sstAntennaCenter.x"));
    if(estimate1y) name.push_back(ParameterName("satellite1", "sstAntennaCenter.y"));
    if(estimate1z) name.push_back(ParameterName("satellite1", "sstAntennaCenter.z"));
    if(estimate2x) name.push_back(ParameterName("satellite2", "sstAntennaCenter.x"));
    if(estimate2y) name.push_back(ParameterName("satellite2", "sstAntennaCenter.y"));
    if(estimate2z) name.push_back(ParameterName("satellite2", "sstAntennaCenter.z"));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline void ParametrizationSatelliteTrackingAntennaCenter::compute(UInt sstType, const std::vector<Time> &time, const Vector &/*sst0*/,
                                                             const Vector &position1, const Vector &position2, const Vector &/*velocity1*/, const Vector &/*velocity2*/,
                                                             const std::vector<Rotary3d> &rotSat1, const std::vector<Rotary3d> &rotSat2, MatrixSliceRef A)
{
  try
  {
    for(UInt i=0; i<time.size(); i++)
    {
      Matrix e = (position2.row(3*i,3)-position1.row(3*i,3)).trans();
      e *= 1./norm(e);
      Matrix B(1,6);
      matMult(-1, e, rotSat1.at(i).matrix(), B.column(0, 3));
      matMult( 1, e, rotSat2.at(i).matrix(), B.column(3, 3));
      UInt idx = 0;
      if(estimate1x) A(i, idx++) = B(0, 0);
      if(estimate1y) A(i, idx++) = B(0, 1);
      if(estimate1z) A(i, idx++) = B(0, 2);
      if(estimate2x) A(i, idx++) = B(0, 3);
      if(estimate2y) A(i, idx++) = B(0, 4);
      if(estimate2z) A(i, idx++) = B(0, 5);
    }

    if(sstType == 1) // range rate
    {
      Polynomial polynomial(degree);
      copy(polynomial.derivative((time.at(1)-time.at(0)).seconds(), A), A);
    }
    else if(sstType == 2) // range acceleration
    {
      Polynomial polynomial(degree);
      copy(polynomial.derivative2nd((time.at(1)-time.at(0)).seconds(), A), A);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

#endif
