/***********************************************/
/**
* @file gnssParametrizationKinematicPositions.h
*
* @brief Position estimation each epoch.
* @see GnssParametrization
*
* @author Torsten Mayer-Guerr
* @date 2021-01-23
*
*/
/***********************************************/

#ifndef __GROOPS_GNSSPARAMETRIZATIONKINEMATICPOSITIONS__
#define __GROOPS_GNSSPARAMETRIZATIONKINEMATICPOSITIONS__

// Latex documentation
#ifdef DOCSTRING_GnssParametrization
static const char *docstringGnssParametrizationKinematicPositions = R"(
\subsection{KinematicPositions}\label{gnssParametrizationType:kinematicPositions}
Estimates the epoch-wise \configFile{outputfilePositions}{instrument}
in an Earth-fixed frame (or in case of LEO satellites in an intertial frame).

The $3\times3$ epoch wise \configFile{outputfileCovarianceEpoch}{instrument}
are computed within
\configClass{GnssProcessing:processingStep:computeCovarianceMatrix}{gnssProcessingStepType:computeCovarianceMatrix}
)";
#endif

/***********************************************/

#include "classes/platformSelector/platformSelector.h"
#include "gnss/gnssParametrization/gnssParametrization.h"

/***** CLASS ***********************************/

/** @brief Position estimation each epoch.
* @ingroup gnssParametrizationGroup
* @see GnssParametrization */
class GnssParametrizationKinematicPositions : public GnssParametrizationBase
{
  Gnss                                        *gnss;
  std::string                                  name;
  PlatformSelectorPtr                          selectReceivers;
  std::vector<Byte>                            selectedReceivers;
  FileName                                     fileNamePositions, fileNameCovariance;
  std::vector<std::vector<GnssParameterIndex>> index; // for each receiver, each epoch
  std::vector<std::vector<Tensor3d>>           cov;

public:
  GnssParametrizationKinematicPositions(Config &config);

  void   init(Gnss *gnss, Parallel::CommunicatorPtr comm) override;
  void   requirements(GnssNormalEquationInfo &normalEquationInfo, std::vector<UInt> &transCount, std::vector<UInt> &transCountEpoch,
                      std::vector<UInt> &recvCount, std::vector<UInt> &recvCountEpoch) override;
  void   initParameter(GnssNormalEquationInfo &normalEquationInfo) override;
  void   aprioriParameter(const GnssNormalEquationInfo &normalEquationInfo, MatrixSliceRef x0) const override;
  void   designMatrix(const GnssNormalEquationInfo &normalEquationInfo, const GnssObservationEquation &eqn, GnssDesignMatrix &A) const override;
  Double updateParameter(const GnssNormalEquationInfo &normalEquationInfo, const_MatrixSliceRef x, const_MatrixSliceRef Wz) override;
  void   updateCovariance(const GnssNormalEquationInfo &normalEquationInfo, const MatrixDistributed &covariance) override;
  void   writeResults(const GnssNormalEquationInfo &normalEquationInfo, const std::string &suffix) const override;
};

/***********************************************/

#endif
