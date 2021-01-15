/***********************************************/
/**
* @file kalmanFilter.cpp
*
* @brief Computes time variable gravity fields using Kalman filter approach
*
* @author Enrico Kurtenbach
* @author Andreas Kvas
* @date 2008-11-01
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
The program computes time variable gravity fields using the Kalman filter approach of

Kurtenbach, E., Eicker, A., Mayer-Gürr, T., Holschneider, M., Hayn, M., Fuhrmann, M., and Kusche, J. (2012).
Improved daily GRACE gravity field solutions using a Kalman smoother. Journal of Geodynamics, 59–60, 39–48.
\url{https://doi.org/10.1016/j.jog.2012.02.006}.

The updated state $\mathbf{x}_t^+$ is determined by solving the least squares adjustment
\begin{equation}
\mathbf{l}_t = \mathbf{A}_t \mathbf{x}_t + \mathbf{e}_t \hspace{25pt} \mathbf{e}_t \sim \mathcal{N}(0, \mathbf{R}_t)\\
\mathbf{B} \mathbf{x}^+_{t-1} = \mathbf{I} \mathbf{x}_t + \mathbf{v}_t\hspace{25pt} \mathbf{v} \sim \mathcal{N}(0,\mathbf{Q} + \mathbf{B} \mathbf{P}^+_{t-1}\mathbf{B}^T).
\end{equation}
In normal equation form this can be written as
\begin{equation}
\hat{\mathbf{x}}_t = \mathbf{x}^+_t = (\mathbf{N}_t + \mathbf{P}^{-^{-1}}_t)^{-1}(\mathbf{n}_t + \mathbf{P}^{-^{-1}}_t \mathbf{x}^-_t),
\end{equation}
where $\mathbf{x}_t^- = \mathbf{B} \mathbf{x}^+_{t-1}$ and $\mathbf{P}_t^{-} = \mathbf{Q} + \mathbf{B} \mathbf{P}^+_{t-1}\mathbf{B}^T$
are the predicted state and its covariance matrix.

The process dynamic $\mathbf{B}, \mathbf{Q}$ is represented as an \reference{autoregressive model}{fundamentals.autoregressiveModel},
and passed to the program through \configFile{inputfileAutoregressiveModel}{matrix}.
The sequence of normal equations $\mathbf{N}_t, \mathbf{n}_t$ are given as list of \configFile{inputfileNormalEquations}{normalEquation},
which can be generated using \configClass{loops}{loopType}.
In the same way, the \file{matrix files}{matrix} for \config{outputfileUpdatedState} and \config{inputfileUpdatedStateCovariance}
can also be specified using \configClass{loops}{loopType}.

If no \configFile{inputfileInitialState}{matrix} is set, a zero vector with appropriate dimensions is used.
The \configFile{inputfileInitialStateCovarianceMatrix}{matrix} however must be given.

See also \program{KalmanBuildNormals}, \program{KalmanSmoother}.
)";

/***********************************************/

#include "programs/program.h"
#include "parser/dataVariables.h"
#include "files/fileMatrix.h"
#include "files/fileNormalEquation.h"
#include "classes/timeSeries/timeSeries.h"
#include "misc/kalmanProcessing.h"

/***** CLASS ***********************************/

/** @brief Computes time variable gravity fields using Kalman filter approach
* @ingroup programsGroup */
class KalmanFilter
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(KalmanFilter, SINGLEPROCESS, "Computes time variable gravity fields using Kalman filter approach", KalmanFilter, NormalEquation)

/***********************************************/

void KalmanFilter::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    std::vector<FileName> fileNameState, fileNameStateCovarianceMatrix;
    FileName fileNameInitialState, fileNameInitialCovariance;
    FileName fileNameArModel;
    std::vector<FileName> fileNameNormals;

    readConfig(config, "outputfileUpdatedState",                   fileNameState,                          Config::MUSTSET,   "kalman/updatedState/x_{loopTime:%D}.txt",                      "estimated state x+ (nx1-matrix)");
    readConfig(config, "outputfileUpdatedStateCovarianceMatrix",   fileNameStateCovarianceMatrix,          Config::OPTIONAL, "kalman/updatedStateCovariance/covariance_{loopTime:%D}.dat",   "estimated state' s covariance matrix Cov(x+)");
    readConfig(config, "inputfileNormalEquations",                 fileNameNormals,                        Config::MUSTSET,  "", "normal equations input file");
    readConfig(config, "inputfileInitialState",                    fileNameInitialState,                   Config::OPTIONAL,  "", "initial state x0");
    readConfig(config, "inputfileInitialStateCovarianceMatrix",    fileNameInitialCovariance,              Config::MUSTSET,   "", "initial state's covariance matrix Cov(x0)");
    readConfig(config, "inputfileAutoregressiveModel",             fileNameArModel,                        Config::MUSTSET,   "", "file name of autoregressive model");
    if(isCreateSchema(config)) return;

    // check input
    // -----------
    const UInt epochCount = fileNameNormals.size();

    if(fileNameState.size() != epochCount)
      throw(Exception("Number of solution file names and normal equations does not match (" + fileNameState.size()%"%i"s + " vs. " + epochCount%"%i"s + ")."));

    if( (fileNameStateCovarianceMatrix.size() > 0) && (fileNameStateCovarianceMatrix.size() != epochCount))
      throw(Exception("Number of covariance matrix file names and normal equations does not match (" + fileNameStateCovarianceMatrix.size()%"%i"s + " vs. " + epochCount%"%i"s + ")."));

    // set up ar model
    // ---------------
    logStatus<<"read autoregressive model from <"<<fileNameArModel<<">"<<Log::endl;
    Matrix tmp;
    readFileMatrix(fileNameArModel, tmp);
    AutoregressiveModel arModel(tmp);

    Matrix B, Q;
    arModel.orderOneRepresentation(B, Q);

    // load initial state:
    // -------------------
    logStatus<<"initialize state's covariance matrix with <"<<fileNameInitialCovariance<<">"<<Log::endl;
    Matrix updatedStateCovariance;
    readFileMatrix(fileNameInitialCovariance, updatedStateCovariance);

    Matrix updatedState(updatedStateCovariance.rows(),1);
    if(!fileNameInitialState.empty())
    {
      logStatus <<"initialize initial state with <"<<fileNameInitialState<<">"<< Log::endl;
      readFileMatrix(fileNameInitialState, updatedState);
    }

    // Run the filter:
    // ---------------
    for(UInt k = 0; k<fileNameNormals.size(); k++)
    {
      Matrix predictedState = B*updatedState;
      Matrix predictedStateCovariance = Q+B*updatedStateCovariance*B.trans();

      try
      {
        NormalEquationInfo info;
        Matrix N, n;
        readFileNormalEquation(fileNameNormals.at(k), info, N, n);

        inverse(predictedStateCovariance);
        updatedStateCovariance = predictedStateCovariance;
        axpy(1.0, N, updatedStateCovariance.slice(0, 0, N.rows(), n.rows()));
        inverse(updatedStateCovariance);
        fillSymmetric(updatedStateCovariance);

        matMult(-1.0, N, predictedState.row(0, N.rows()), n);
        updatedState = predictedState;
        matMult(1.0, updatedStateCovariance.column(0, n.rows()), n, updatedState);
      }
      catch(std::exception &e)
      {
        logWarning<<e.what()<<Log::endl;
        updatedState = predictedState;
        updatedStateCovariance = predictedStateCovariance;
      }

      logStatus <<"write updated state to <"<<fileNameState.at(k)<<">"<<Log::endl;
      writeFileMatrix(fileNameState.at(k), updatedState);
      if(!fileNameStateCovarianceMatrix.empty())
      {
        logStatus <<"write updated state covariance to <"<<fileNameStateCovarianceMatrix.at(k)<<">"<<Log::endl;
        writeFileMatrix(fileNameStateCovarianceMatrix.at(k), updatedStateCovariance);
      }
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
