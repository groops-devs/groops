/***********************************************/
/**
* @file slrParametrizationStaticPositions.h
*
* @brief Position estimation with no-net constraints.
* @see SlrParametrization
*
* @author Torsten Mayer-Guerr
* @date 2022-04-28
*
*/
/***********************************************/

#ifndef __GROOPS_SLRPARAMETRIZATIONSTATICPOSITIONS__
#define __GROOPS_SLRPARAMETRIZATIONSTATICPOSITIONS__

// Latex documentation
#ifdef DOCSTRING_SlrParametrization
static const char *docstringSlrParametrizationStaticPositions = R"(
\subsection{StaticPositions}\label{slrParametrizationType:staticPositions}
Estimates a static position for all
\configClass{selectReceivers}{platformSelectorType} in the terrestrial frame.

No-net constraints can be applied for a subset of stations,
\configClass{selectNoNetReceivers}{platformSelectorType}, with a
standard deviation of \config{noNetTranslationSigma} and \config{noNetRotationSigma} and \config{noNetScaleSigma} and \config{noNetScaleSigma}.
If the template \configFile{inputfileNoNetPositions}{stringList} is provided
the constraints are applied relatively to these positions. Only stations with an existing position file
are considered. Without \configFile{inputfileNoNetPositions}{stringList}
the constraints are applied towards the apriori values from
\configClass{SlrProcessing:station}{slrStationGeneratorType}.
As a single corrupted station position can disturb the no-net conditions,
the rotation/translation parameters are estimated in a
\reference{robust least squares adjustment}{fundamentals.robustLeastSquares}
beforehand. The computed weight matrix is used to downweight corrupted stations
in the constraint equations.

In case you want to align to an ITRF/ILRS reference frame, precise coordinates can be
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
#include "slr/slrParametrization/slrParametrization.h"

/***** CLASS ***********************************/

/** @brief Position estimation with no-net constraints.
* @ingroup slrParametrizationGroup
* @see SlrParametrization */
class SlrParametrizationStaticPositions : public SlrParametrizationBase
{
  Slr                           *slr;
  std::string                    name, nameConstraint;
  PlatformSelectorPtr            selectorStations, selectorNoNetStations;
  std::vector<Byte>              selectedStations, selectedNoNetStations;
  FileName                       fileNameGrid, fileNamePosition, fileNameNoNetPositions;
  Bool                           applyConstraint;
  Double                         sigmaNoNetTranslation, sigmaNoNetRotation, sigmaNoNetScale;
  Double                         huber, huberPower;
  std::vector<SlrParameterIndex> index; // for each station
  std::vector<Vector3d>          pos, pos0, noNetPos;
  mutable Matrix                 noNetEstimator;

public:
  SlrParametrizationStaticPositions(Config &config);

  void   init(Slr *slr, const std::vector<const SlrParametrizationGravityField*> &paramGravityField) override;
  void   initParameter(SlrNormalEquationInfo &normalEquationInfo) override;
  void   aprioriParameter(const SlrNormalEquationInfo &normalEquationInfo, MatrixSliceRef x0) const override;
  void   designMatrix(const SlrNormalEquationInfo &normalEquationInfo, const SlrObservationEquation &eqn, SlrDesignMatrix &A) const override;
  void   constraints(const SlrNormalEquationInfo &normalEquationInfo, MatrixDistributed &normals, std::vector<Matrix> &n, Double &lPl, UInt &obsCount) const override;
  Double updateParameter(const SlrNormalEquationInfo &normalEquationInfo, const_MatrixSliceRef x, const_MatrixSliceRef Wz) override;
  void   writeResults(const SlrNormalEquationInfo &normalEquationInfo, const std::string &suffix) const override;
};

/***********************************************/

#endif
