/***********************************************/
/**
* @file kalmanSmoother.cpp
*
* @brief Computes time variable gravity fields using Kalman smoother
*
* @author Enrico Kurtenbach
* @author Andreas Kvas
* @date 2008-11-01
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
Apply the Rauch-Tung-Striebel smoother to a gravity field time series computed by \program{KalmanFilter}.
This is the implementation of the approach presented in

Kurtenbach, E., Eicker, A., Mayer-Gürr, T., Holschneider, M., Hayn, M., Fuhrmann, M., and Kusche, J. (2012).
Improved daily GRACE gravity field solutions using a Kalman smoother. Journal of Geodynamics, 59–60, 39–48.
\url{https://doi.org/10.1016/j.jog.2012.02.006}.

The result has zero phase and the squared magnitude response of \configFile{inputfileAutoregressiveModel}{matrix}
(see \reference{autoregressiveModel}{fundamentals.autoregressiveModel} for details).
\configFile{inputfileUpdatedState}{matrix} and \configFile{inputfileUpdatedStateCovariance}{matrix}
are the output of a \program{KalmanFilter} forward sweep.
The matrix files for\configFile{outputfileUpdatedState}{matrix}, \configFile{inputfileUpdatedState}{matrix}
and \configFile{inputfileUpdatedStateCovariance}{matrix} can also be specified using \configClass{loops}{loopType}.

See also \program{KalmanBuildNormals}, \program{KalmanFilter} and \program{KalmanSmootherLeastSquares}.
)";

/***********************************************/

#include "programs/program.h"
#include "parser/dataVariables.h"
#include "files/fileMatrix.h"
#include "classes/timeSeries/timeSeries.h"
#include "misc/kalmanProcessing.h"

/***** CLASS ***********************************/

/** @brief Computes time variable gravity fields using Kalman filter approach
* @ingroup programsGroup */
class KalmanSmoother
{
public:
  void run(Config &config);
};

GROOPS_REGISTER_PROGRAM(KalmanSmoother, SINGLEPROCESS, "Computes time variable gravity fields using Kalman smoother approach", KalmanFilter, NormalEquation)

/***********************************************/

void KalmanSmoother::run(Config &config)
{
  try
  {
    std::vector<FileName> fileNameSmoothedState,  fileNameSmoothedCovariance, fileNameUpdatedState, fileNameUpdatedCovariance;
    FileName fileNameArModel;

    readConfig(config, "outputfileState",                         fileNameSmoothedState,            Config::MUSTSET,   "kalman/smoothedState/x_{loopTime:%D}.txt",                    "estimated parameters (nx1-matrix)");
    readConfig(config, "outputfileStateCovarianceMatrix",         fileNameSmoothedCovariance,       Config::OPTIONAL,  "kalman/smoothedStateCovariance/covariance_{loopTime:%D}.dat", "estimated parameters' covariance matrix");
    readConfig(config, "inputfileUpdatedState",                   fileNameUpdatedState,             Config::MUSTSET,   "kalman/updatedState/x_{loopTime:%D}.txt", "");
    readConfig(config, "inputfileUpdatedStateCovarianceMatrix",   fileNameUpdatedCovariance,        Config::MUSTSET,   "kalman/updatedStateCovariance/covariance_{loopTime:%D}.dat", "");
    readConfig(config, "inputfileAutoregressiveModel",            fileNameArModel,                  Config::MUSTSET,   "", "file name of autoregressive model");
    if(isCreateSchema(config)) return;

    // check input
    // -----------
    const UInt epochCount = fileNameSmoothedState.size();

    if(fileNameUpdatedState.size() != epochCount)
      throw(Exception("Number of update state file names and state vectors does not match (" + fileNameUpdatedState.size()%"%i"s + " vs. " + epochCount%"%i"s + ")."));

    if(fileNameUpdatedCovariance.size() != epochCount)
      throw(Exception("Number of update state covariance file names and state vectors does not match (" + fileNameUpdatedCovariance.size()%"%i"s + " vs. " + epochCount%"%i"s + ")."));

    if( (fileNameSmoothedCovariance.size() > 0) && (fileNameSmoothedCovariance.size() != epochCount))
      throw(Exception("Number of smoothed state covariance file names and state vectors does not match (" + fileNameSmoothedCovariance.size()%"%i"s + " vs. " + epochCount%"%i"s + ")."));

    // compute AR model
    // ----------------
    logStatus<<"read autoregressive model from <"<<fileNameArModel<<">"<<Log::endl;
    Matrix tmp;
    readFileMatrix(fileNameArModel, tmp);
    AutoregressiveModel arModel(tmp);

    Matrix B, Q;
    arModel.orderOneRepresentation(B, Q);

    // Initialize backward smoother:
    Matrix smoothedState, smoothedStateCovariance;
    Matrix updatedState, updatedStateCovariance;
    Matrix predictedState, predictedStateCovariance;

    logStatus <<"initialize state with <"<<fileNameUpdatedState.back()<<"> and <"<<fileNameUpdatedCovariance.back()<<">"<<Log::endl;
    readFileMatrix(fileNameUpdatedState.back(), updatedState);
    readFileMatrix(fileNameUpdatedCovariance.back(), updatedStateCovariance);

    smoothedState = updatedState;
    smoothedStateCovariance = updatedStateCovariance;
    logStatus <<"write smoothed state to <"<<fileNameSmoothedState.back()<<">"<<Log::endl;
    writeFileMatrix(fileNameSmoothedState.back(), smoothedState);
    if(fileNameSmoothedCovariance.size() > 0)
    {
      logStatus <<"write smoothed state to <"<<fileNameSmoothedCovariance.back()<<">"<<Log::endl;
      writeFileMatrix(fileNameSmoothedCovariance.back(), smoothedStateCovariance);
    }

    for(UInt k=epochCount-1; k>0; k--)
    {
      readFileMatrix(fileNameUpdatedState.at(k-1), updatedState);
      readFileMatrix(fileNameUpdatedCovariance.at(k-1), updatedStateCovariance);

      predictedState = B*updatedState;

      // predicted state covariance
      cholesky(updatedStateCovariance);
      Matrix smootherGain = B;
      triangularMult(1.0, updatedStateCovariance, smootherGain.trans());

      predictedStateCovariance = Q;
      rankKUpdate(1.0, smootherGain.trans(), predictedStateCovariance);

      // transposed gain matrix K
      triangularMult(1.0, updatedStateCovariance.trans(), smootherGain.trans());
      solveInPlace(predictedStateCovariance, smootherGain);

      // smoothed state
      smoothedState = updatedState + smootherGain.trans()*(smoothedState-predictedState);

      // smoothed state covariance
      zeroUnusedTriangle(predictedStateCovariance);
      predictedStateCovariance.setType(Matrix::GENERAL);
      rankKUpdate(-1.0, predictedStateCovariance, smoothedStateCovariance);
      fillSymmetric(smoothedStateCovariance);
      smoothedStateCovariance.setType(Matrix::GENERAL);

      smoothedStateCovariance = smootherGain.trans()*smoothedStateCovariance*smootherGain;
      smoothedStateCovariance.setType(Matrix::SYMMETRIC);
      zeroUnusedTriangle(updatedStateCovariance);
      updatedStateCovariance.setType(Matrix::GENERAL);
      rankKUpdate(1.0, updatedStateCovariance, smoothedStateCovariance);

      logStatus <<"write smoothed state to <"<<fileNameSmoothedState.at(k-1)<<">"<<Log::endl;
      writeFileMatrix(fileNameSmoothedState.at(k-1), smoothedState);
      if(!fileNameSmoothedCovariance.empty())
      {
        logStatus <<"write smoothed state covariance to <"<<fileNameSmoothedCovariance.at(k-1)<<">"<<Log::endl;
        writeFileMatrix(fileNameSmoothedCovariance.at(k-1), smoothedStateCovariance);
      }
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
