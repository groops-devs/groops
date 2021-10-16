/***********************************************/
/**
* @file gnssProcessingStepComputeCovarianceMatrix.h
*
* @brief GNSS processing step: ComputeCovarianceMatrix.
*
* @author Torsten Mayer-Guerr
* @date 2021-09-05
*
*/
/***********************************************/

#ifndef __GROOPS_GNSSPROCESSINGSTEPCOMPUTECOVARIANCEMATRIX__
#define __GROOPS_GNSSPROCESSINGSTEPCOMPUTECOVARIANCEMATRIX__

// Latex documentation
#ifdef DOCSTRING_GnssProcessingStep
static const char *docstringGnssProcessingStepComputeCovarianceMatrix = R"(
\subsection{ComputeCovarianceMatrix}\label{gnssProcessingStepType:computeCovarianceMatrix}
Accumulates the normal equations and computes the covariance matrix as inverse of the normal matrix.
It is not the full inverse but only the elements which are set in the normal matrix
(see  \configClass{gnssProcessingStep:selectNormalsBlockStructure}{gnssProcessingStepType:selectNormalsBlockStructure})
are computed. The matrix is passed to the \configClass{parametrizations}{gnssParametrizationType}.
Only used in \configClass{parametrizations:kinematicPositions}{gnssParametrizationType:kinematicPositions}
to get the epoch wise covariance information at the moment.
)";
#endif

/***********************************************/

#include "config/config.h"
#include "gnss/gnssProcessingStep/gnssProcessingStep.h"

/***** CLASS ***********************************/

/** @brief GNSS processing step: ComputeCovarianceMatrix.
* @ingroup gnssProcessingStepGroup
* @see GnssProcessingStep */
class GnssProcessingStepComputeCovarianceMatrix : public GnssProcessingStepBase
{
  Bool   computeResiduals, adjustSigma0, computeWeights;
  Double huber, huberPower;
  Double convergenceThreshold;
  UInt   iterCount;

public:
  GnssProcessingStepComputeCovarianceMatrix(Config &/*config*/) {}
  void process(GnssProcessingStep::State &state) override;
};

/***********************************************/

inline void GnssProcessingStepComputeCovarianceMatrix::process(GnssProcessingStep::State &state)
{
  try
  {
    logStatus<<"=== compute covariance matrix ==============================="<<Log::endl;
    state.buildNormals(FALSE/*constraintsOnly*/, FALSE/*solveEpochParameters*/);
    logStatus<<"cholesky"<<Log::endl;
    state.normals.cholesky(TRUE/*timing*/);
    logStatus<<"sparse inverse"<<Log::endl;
    state.normals.cholesky2SparseInverse();
    state.gnss->updateCovariance(state.normalEquationInfo, state.normals);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

#endif
