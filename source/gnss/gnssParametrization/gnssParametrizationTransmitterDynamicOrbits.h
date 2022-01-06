/***********************************************/
/**
* @file gnssParametrizationTransmitterDynamicOrbits.h
*
* @brief Orbits by variational equations.
* @see GnssParametrization
*
* @author Torsten Mayer-Guerr
* @author Sebastian Strasser
* @date 2021-01-23
*
*/
/***********************************************/

#ifndef __GROOPS_GNSSPARAMETRIZATIONTRANSMITTERDYNAMICORBITS__
#define __GROOPS_GNSSPARAMETRIZATIONTRANSMITTERDYNAMICORBITS__

// Latex documentation
#ifdef DOCSTRING_GnssParametrization
static const char *docstringGnssParametrizationTransmitterDynamicOrbits = R"(
\subsection{TransmitterDynamicOrbits}\label{gnssParametrizationType:transmitterDynamicOrbits}
Same as \configClass{leoDynamicOrbits}{gnssParametrizationType:leoDynamicOrbits} but
for transmitting GNSS satellites.
For more details see \reference{orbit integration}{cookbook.gnssNetwork:orbitIntegration}.
)";
#endif

/***********************************************/

#include "base/import.h"
#include "config/config.h"
#include "gnss/gnss.h"
#include "classes/ephemerides/ephemerides.h"
#include "classes/parametrizationAcceleration/parametrizationAcceleration.h"
#include "gnss/gnssTransceiverSelector/gnssTransceiverSelector.h"
#include "gnss/gnssParametrization/gnssParametrization.h"

/***** CLASS ***********************************/

/** @brief Orbits by variational equations.
* @ingroup gnssParametrizationGroup
* @see GnssParametrization */
class GnssParametrizationTransmitterDynamicOrbits : public GnssParametrizationBase
{
  class Parameter
  {
  public:
    GnssTransmitterPtr         trans;
    GnssParameterIndex         index;
    std::vector<ParameterName> parameterNames;
    std::vector<Time>          times;
    Matrix                     PosDesign, VelDesign;
    Vector                     x;
    Polynomial                 polynomial;
  };

  Gnss                          *gnss;
  std::string                    name;
  GnssTransceiverSelectorPtr     selectTransmitters;
  FileName                       fileNameOrbit, fileNameParameter, fileNameVariational;
  std::vector<Time>              pulses;
  ParametrizationAccelerationPtr parametrizationAcceleration;
  EphemeridesPtr                 ephemerides;
  Double                         minEstimableEpochsRatio;
  UInt                           integrationDegree, interpolationDegree;
  std::vector<Parameter*>        parameters;

public:
  GnssParametrizationTransmitterDynamicOrbits(Config &config);
 ~GnssParametrizationTransmitterDynamicOrbits();

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
