/***********************************************/
/**
* @file gnssParametrizationLeoDynamicOrbits.h
*
* @brief Orbits by variational equations.
* @see GnssParametrization
*
* @author Torsten Mayer-Guerr
* @date 2021-01-23
*
*/
/***********************************************/

#ifndef __GROOPS_GNSSPARAMETRIZATIONLEODYNAMICORBITS__
#define __GROOPS_GNSSPARAMETRIZATIONLEODYNAMICORBITS__

// Latex documentation
#ifdef DOCSTRING_GnssParametrization
static const char *docstringGnssParametrizationLeoDynamicOrbits = R"(
\subsection{LeoDynamicOrbits}\label{gnssParametrizationType:leoDynamicOrbits}
The estimation of (reduced) dynamic orbits is formulated as variational equations.
It is based on \configFile{inputfileVariational}{variationalEquation} calculated with \program{PreprocessingVariationalEquation}.
Necessary integrations are performed by integrating a moving interpolation polynomial of degree \config{integrationDegree}.
The \configClass{parametrizationAcceleration}{parametrizationAccelerationType} must include at least those
parameters that were estimated in \program{PreprocessingVariationalEquationOrbitFit}.
Additional \configClass{stochasticPulse}{timeSeriesType} parameters can be set up to reduce orbit mismodeling.
If not enough epochs with observations are available (\config{minEstimableEpochsRatio}) the LEO satellite is disabled.
)";
#endif

/***********************************************/

#include "base/import.h"
#include "config/config.h"
#include "gnss/gnss.h"
#include "classes/ephemerides/ephemerides.h"
#include "classes/parametrizationAcceleration/parametrizationAcceleration.h"
#include "classes/platformSelector/platformSelector.h"
#include "gnss/gnssParametrization/gnssParametrization.h"

/***** CLASS ***********************************/

/** @brief Orbits by variational equations.
* @ingroup gnssParametrizationGroup
* @see GnssParametrization */
class GnssParametrizationLeoDynamicOrbits : public GnssParametrizationBase
{
  class Parameter
  {
  public:
    GnssReceiverPtr                recv;
    GnssParameterIndex             index;
    std::vector<ParameterName>     parameterNames;
    Vector                         x;
    // for each arc:
    std::vector<UInt>              startEpoch, endEpoch;
    std::vector<std::vector<Time>> times;
    std::vector<Matrix>            PosDesign, VelDesign;
    std::vector<Vector>            pos, vel;
    std::vector<Polynomial>        polynomial;
  };

  Gnss                          *gnss;
  std::string                    name;
  PlatformSelectorPtr            selectReceivers;
  FileName                       fileNameOrbit, fileNameParameter, fileNameVariational;
  std::vector<Time>              pulses;
  ParametrizationAccelerationPtr parametrizationAcceleration;
  EphemeridesPtr                 ephemerides;
  Double                         minEstimableEpochsRatio;
  UInt                           integrationDegree, interpolationDegree;
  std::vector<Parameter*>        parameters; // for each receiver

public:
  GnssParametrizationLeoDynamicOrbits(Config &config);
 ~GnssParametrizationLeoDynamicOrbits();

  void   init(Gnss *gnss, Parallel::CommunicatorPtr comm) override;
  void   requirements(GnssNormalEquationInfo &normalEquationInfo, std::vector<UInt> &transCount, std::vector<UInt> &transCountEpoch,
                      std::vector<UInt> &recvCount, std::vector<UInt> &recvCountEpoch) override;
  void   initParameter(GnssNormalEquationInfo &normalEquationInfo) override;
  void   aprioriParameter(const GnssNormalEquationInfo &normalEquationInfo, MatrixSliceRef x0) const override;
  void   designMatrix(const GnssNormalEquationInfo &normalEquationInfo, const GnssObservationEquation &eqn, GnssDesignMatrix &A) const override;
  Double updateParameter(const GnssNormalEquationInfo &normalEquationInfo, const_MatrixSliceRef x, const_MatrixSliceRef Wz) override;
  void   writeResults(const GnssNormalEquationInfo &normalEquationInfo, const std::string &suffix) const override;
};

/***********************************************/

#endif
