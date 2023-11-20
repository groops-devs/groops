/***********************************************/
/**
* @file gnssParametrizationStaticPositions.h
*
* @brief Position estimation with no-net constraints.
* @see GnssParametrization
*
* @author Torsten Mayer-Guerr
* @author Sebastian Strasser
* @date 2021-01-23
*
*/
/***********************************************/

#ifndef __GROOPS_GNSSPARAMETRIZATIONSTATICPOSITIONS__
#define __GROOPS_GNSSPARAMETRIZATIONSTATICPOSITIONS__

// Latex documentation
#ifdef DOCSTRING_GnssParametrization
static const char *docstringGnssParametrizationStaticPositions = R"(
\subsection{StaticPositions}\label{gnssParametrizationType:staticPositions}
Estimates a static position for all
\configClass{selectReceivers}{platformSelectorType} in the terrestrial frame.

No-net constraints can be applied for a subset of stations,
\configClass{selectNoNetReceivers}{platformSelectorType}, with a
standard deviation of \config{noNetTranslationSigma} and \config{noNetRotationSigma}.
If the template \configFile{inputfileNoNetPositions}{stringList} is provided
the constraints are applied relatively to these positions. Only stations with an existing position file
are considered. Without \configFile{inputfileNoNetPositions}{stringList}
the constraints are applied towards the apriori values from
\configClass{GnssProcessing:receiver}{gnssReceiverGeneratorType}.
As a single corrupted station position can disturb the no-net conditions,
the rotation/translation parameters are estimated in a
\reference{robust least squares adjustment}{fundamentals.robustLeastSquares}
beforehand. The computed weight matrix is used to downweight corrupted stations
in the constraint equations.

In case you want to align to an ITRF/IGS reference frame, precise coordinates can be
generated with \program{Sinex2StationPositions}.

The \file{parameter names}{parameterName} are
\begin{itemize}
\item \verb|<station>:position.x::|,
\item \verb|<station>:position.y::|,
\item \verb|<station>:position.z::|.
\end{itemize}
)";
#endif

/***********************************************/

#include "classes/platformSelector/platformSelector.h"
#include "gnss/gnssParametrization/gnssParametrization.h"

/***** CLASS ***********************************/

/** @brief Position estimation with no-net constraints.
* @ingroup gnssParametrizationGroup
* @see GnssParametrization */
class GnssParametrizationStaticPositions : public GnssParametrizationBase
{
  Gnss                           *gnss;
  std::string                     name, nameConstraint;
  PlatformSelectorPtr             selectReceivers, selectNoNetReceivers;
  std::vector<Byte>               selectedReceivers, selectedNoNetReceivers;
  FileName                        fileNameGrid, fileNamePosition, fileNameNoNetPositions;
  Bool                            applyConstraint;
  Double                          sigmaNoNetTranslation, sigmaNoNetRotation;
  Double                          huber, huberPower;
  std::vector<GnssParameterIndex> index; // for each receiver
  std::vector<Vector3d>           pos, pos0, noNetPos;
  mutable Matrix                  noNetEstimator;

public:
  GnssParametrizationStaticPositions(Config &config);

  void   init(Gnss *gnss, Parallel::CommunicatorPtr comm) override;
  void   initParameter(GnssNormalEquationInfo &normalEquationInfo) override;
  void   aprioriParameter(const GnssNormalEquationInfo &normalEquationInfo, MatrixSliceRef x0) const override;
  void   designMatrix(const GnssNormalEquationInfo &normalEquationInfo, const GnssObservationEquation &eqn, GnssDesignMatrix &A) const override;
  void   constraints(const GnssNormalEquationInfo &normalEquationInfo, MatrixDistributed &normals, std::vector<Matrix> &n, Double &lPl, UInt &obsCount) const override;
  Double updateParameter(const GnssNormalEquationInfo &normalEquationInfo, const_MatrixSliceRef x, const_MatrixSliceRef Wz) override;
  void   writeResults(const GnssNormalEquationInfo &normalEquationInfo, const std::string &suffix) const override;
};

/***********************************************/

#endif
