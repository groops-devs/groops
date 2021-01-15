/***********************************************/
/**
* @file kalmanSmootherLeastSquares.cpp
*
* @brief Computes smoothed estimates of time variable gravity field using a constraint least squares adjustment.
*
* @author Andreas Kvas
* @date 2015-10-08
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program estimates temporal gravity field variations with a constraint least squares adjustment.
Prior information is introduced by means of a \configClass{autoregressiveModelSequence}{autoregressiveModelSequenceType}
which represent a stationary random process (see the \reference{autoregressive model description}{fundamentals.autoregressiveModel}) for details.

The output files for the estimated gravity field (\configFile{outputfileSolution}{matrix}), the
corresponding standard deviations (\configFile{outputfileSigmax}{matrix}) and the full covariance matrix
(\configFile{outputfileCovariance}{matrix}) can be specified using \configClass{loops}{loopType}.
Similarly, the \configFile{inputfileNormalEquations}{normalEquation}
can also be specified using \configClass{loops}{loopType}.

See also \program{KalmanBuildNormals}, \program{KalmanFilter} and\program{KalmanSmoother}
)";

/***********************************************/

#include "programs/program.h"
#include "parser/dataVariables.h"
#include "parallel/matrixDistributed.h"
#include "files/fileMatrix.h"
#include "files/fileNormalEquation.h"
#include "classes/timeSeries/timeSeries.h"
#include "misc/kalmanProcessing.h"

/***** CLASS ***********************************/

/** @brief Computes smoothed estimates of time variable gravity field using a constraint least squares adjustment.
* RTS smoother implemented as least squares adjustment.
* @ingroup programsGroup */
class KalmanSmootherLeastSquares
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(KalmanSmootherLeastSquares, PARALLEL, "Smoothed time variable gravity field by least squares adjustment", KalmanFilter, NormalEquation)

/***********************************************/

void KalmanSmootherLeastSquares::run(Config &config, Parallel::CommunicatorPtr comm)
{
  try
  {
    std::vector<FileName> fileNameSolution, fileNameSigmax, fileNameCovariance, fileNameNormals;
    AutoregressiveModelSequencePtr arSequence;

    readConfig(config, "outputfileSolution",             fileNameSolution,               Config::MUSTSET,   "kalman/smoothedState/x_{loopTime:%D}.txt", "file name of solution vector (use time tags)");
    readConfig(config, "outputfileSigmax"  ,             fileNameSigmax,                 Config::OPTIONAL,  "kalman/smoothedStateSigma/sigmax_{loopTime:%D}.txt", "file name of sigma vector (use time tags)");
    readConfig(config, "outputfileCovariance"  ,         fileNameCovariance,             Config::OPTIONAL,  "kalman/smoothedStateCovariance/covariance_{loopTime:%D}.dat", "file name of full covariance matrix (use time tags)");
    readConfig(config, "inputfileNormalEquations"  ,     fileNameNormals,                Config::MUSTSET,   "", "input normal equations (loopTime will be expanded)");
    readConfig(config, "autoregressiveModelSequence",    arSequence,                     Config::MUSTSET,   "", "file containing AR model for spatiotemporal constraint");
    if(isCreateSchema(config)) return;

    // check input
    // -----------
    const UInt epochCount = fileNameNormals.size();

    if(fileNameSolution.size() != epochCount)
      throw(Exception("Number of solution file names and normal equations does not match (" + fileNameSolution.size()%"%i"s + " vs. " + epochCount%"%i"s + ")."));

    if( (fileNameSigmax.size() > 0) && (fileNameSigmax.size() != epochCount))
      throw(Exception("Number of standard deviation file names and normal equations does not match (" + fileNameSigmax.size()%"%i"s + " vs. " + epochCount%"%i"s + ")."));

    if( (fileNameCovariance.size() > 0) && (fileNameCovariance.size() != epochCount))
      throw(Exception("Number of covariance matrix file names and normal equations does not match (" + fileNameCovariance.size()%"%i"s + " vs. " + epochCount%"%i"s + ")."));

    // set up normals of process
    // -------------------------
    logInfo<<"process is modelled by an AR("<<arSequence->maximumOrder()<<") model"<<Log::endl;
    const UInt stateCount = arSequence->dimension();
    logInfo<<"dimension of state vector ("<<stateCount<<"x"<<1<<")"<<Log::endl;

    // set up normal equations
    // -----------------------
    logStatus<<"set up normal equations"<<Log::endl;
    std::vector<UInt> blockIndex(1, 0);
    for(UInt k = 0; k<epochCount; k++)
      blockIndex.push_back(blockIndex.back() + stateCount);

    MatrixDistributed normals;
    normals.initEmpty(blockIndex, comm);
    Vector rhs(blockIndex.back());

    // each process reads the normals it requires for filling the main diagonal
    Double lPlSum = 0.0;
    UInt obsCountSum = 0;
    Single::forEach(epochCount, [&](UInt k)
    {
      normals.setBlock(k, k);
      if(normals.isMyRank(k, k))
      {
        NormalEquationInfo info;
        Matrix satelliteNormals, satelliteRightHandSide;
        try
        {
          readFileNormalEquation(fileNameNormals.at(k), info,  satelliteNormals, satelliteRightHandSide);
        }
        catch(std::exception &/*e*/)
        {
          logWarning<<"Unable to read normal equation from <"<<fileNameNormals.at(k)<<">"<<Log::endl;
          return;
        }
        lPlSum += info.lPl(0);
        obsCountSum += info.observationCount;
        copy(satelliteNormals, normals.N(k, k));
        copy(satelliteRightHandSide, rhs.row(k*stateCount, stateCount));
      }
    });

    Parallel::reduceSum(lPlSum, 0, comm);
    Parallel::reduceSum(obsCountSum, 0, comm);
    Parallel::reduceSum(rhs, 0, comm);
    Parallel::broadCast(rhs, 0, comm);

    // add process normals
    // -------------------
    logStatus<<"add normals of pseudo-observations"<<Log::endl;
    auto index = arSequence->distributedNormalsBlockIndex(epochCount);
    Single::forEach(index.size(), [&](UInt k)
    {
      normals.setBlock(index[k].first, index[k].second);
      if(normals.isMyRank(index[k].first, index[k].second))
        arSequence->distributedNormalsBlock(normals.blockCount(), index[k].first, index[k].second, normals.N(index[k].first, index[k].second));
    });

    logStatus<<"solve normal equation system"<<Log::endl;
    Parallel::barrier(comm);
    Matrix x = normals.solve(rhs, TRUE/*timing*/); // normals now holds Cholesky R
    if(Parallel::isMaster(comm))
      logInfo<<"  a posteriori sigma = "<<sqrt((lPlSum-inner(x,rhs))/obsCountSum)<<Log::endl;

    if(!fileNameSigmax.empty() || !fileNameCovariance.empty())
    {
      logStatus<<"compute sparse inverse of normal equations"<<Log::endl;
      normals.cholesky2SparseInverse();
    }
    Parallel::barrier(comm);

    // write results to disk
    // ---------------------
    for(UInt k=0; k<epochCount; k++)
    {
      if(Parallel::isMaster(comm))
      {
        logStatus<<"write solution vector to <"<< fileNameSolution.at(k) <<">."<<Log::endl;
        writeFileMatrix(fileNameSolution.at(k), x.row(k*stateCount, stateCount));
      }

      if((fileNameSigmax.size()>0) && normals.isMyRank(k, k))
      {
        Vector sigmax(normals.N(k,k).rows());
        for(UInt l = 0; l<sigmax.rows(); l++)
          sigmax(l) = std::sqrt(normals.N(k, k)(l, l));
        writeFileMatrix(fileNameSigmax.at(k), sigmax);
      }

      if((fileNameCovariance.size()>0) && normals.isMyRank(k, k))
        writeFileMatrix(fileNameCovariance.at(k), normals.N(k, k));
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
