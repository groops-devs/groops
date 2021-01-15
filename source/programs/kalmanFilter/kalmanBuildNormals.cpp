/***********************************************/
/**
* @file kalmanBuildNormals.cpp
*
* @brief Accumulate normals for sub monthly data sets.
*
* @author Torsten Mayer-Guerr
* @author Andreas Kvas
* @date 2009-10-10
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program sets up normal equations based on \configClass{observation}{observationType}
for short-term gravity field variations.
It computes the normal equations based on the intervals $i \in \{1, ..., N\}$ given in the \configFile{arcList}{arcList}.
It sets up the least squares adjustment
\begin{equation}
    \begin{bmatrix}
    \mathbf{l}_1 \\
    \mathbf{l}_2 \\
    \vdots \\
    \mathbf{l}_N \\
  \end{bmatrix}
  =
  \begin{bmatrix}
    \mathbf{A}_1  &  & & \\
    & \mathbf{A}_2  & &\\
    &  & \ddots & \\
    & & & \mathbf{A}_N \\
  \end{bmatrix}
  \begin{bmatrix}
    \mathbf{x}^{(1)} \\
    \mathbf{x}^{(2)} \\
    \vdots \\
    \mathbf{x}^{(N)} \\
  \end{bmatrix}
  +
  \begin{bmatrix}
    \mathbf{e}_1 \\
    \mathbf{e}_2 \\
    \vdots \\
    \mathbf{e}_N \\
  \end{bmatrix},
\end{equation}
and subsequently computes the normal equations $\mathbf{N}_i, \mathbf{n}_i$ for each interval.
If \config{eliminateNonGravityParameters} is true, all non-gravity parameters are eliminated before the normals
are written to \configFile{outputfileNormalEquation}{normalEquation}.
For each time interval in \config{arcList} a single \file{normal equation file}{normalEquation} is written.

This program computes the input normals for \program{KalmanFilter} and \program{KalmanSmootherLeastSquares}.
)";

/***********************************************/

#include "programs/program.h"
#include "parser/dataVariables.h"
#include "files/fileArcList.h"
#include "files/fileNormalEquation.h"
#include "classes/observation/observation.h"

/***** CLASS ***********************************/

/** @brief Accumulate normals for sub monthly data sets.
* @ingroup programsGroup */
class KalmanBuildNormals
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);

private:
  ObservationPtr    observation;
  std::vector<UInt> arcsInterval;
  std::vector<Time> timesInterval;

  // normal equation system (for each interval)
  std::vector<UInt>    obsCount;
  std::vector<Matrix>  N,  n;
  std::vector<Vector>  lPl;

  void computeArc(UInt arcNo);
};

GROOPS_REGISTER_PROGRAM(KalmanBuildNormals, PARALLEL, "accumulate normals for sub monthly data sets.", KalmanFilter, NormalEquation)

/***********************************************/

void KalmanBuildNormals::run(Config &config, Parallel::CommunicatorPtr comm)
{
  try
  {
    FileName      fileNameNormals;
    FileName      fileNameArcList;
    Bool          eliminatePararmeter;

    renameDeprecatedConfig(config, "inputfileNormalequation", "inputfileNormalEquation", date2time(2020, 6, 3));
    renameDeprecatedConfig(config, "arcList",                 "inputfileArcList",        date2time(2020, 7, 7));

    readConfig(config, "outputfileNormalEquation",      fileNameNormals,     Config::MUSTSET,  "normals/normals_{loopTime:%D}.dat", "outputfile for normal equations");
    readConfig(config, "observation",                   observation,         Config::MUSTSET,  "", "");
    readConfig(config, "inputfileArcList",              fileNameArcList,     Config::MUSTSET,  "",  "list to correspond points of time to arc numbers");
    readConfig(config, "eliminateNonGravityParameters", eliminatePararmeter, Config::DEFAULT,  "1", "eliminate additional parameters from normals, 0: all parameter are saved");
    if(isCreateSchema(config)) return;

    // =======================

    const UInt arcCount = observation->arcCount();

    // read arc list
    // -------------
    logStatus<<"read arc list <"<<fileNameArcList<<">"<<Log::endl;
    readFileArcList(fileNameArcList, arcsInterval, timesInterval);
    if(arcsInterval.back()>arcCount)
      throw(Exception("count of arcs differ in observation and arcList"));

    VariableList fileNameVariableList;
    addTimeVariables(fileNameVariableList);

    // =======================

    // setup observation equations
    // ---------------------------
    logStatus<<"set up observation equations"<<Log::endl;
    // init normals
    N.resize(arcsInterval.size()-1);
    n.resize(arcsInterval.size()-1);
    lPl.resize(arcsInterval.size()-1);
    obsCount.resize(arcsInterval.size()-1, 0);
    Parallel::forEachInterval(arcCount, arcsInterval, [this](UInt arcNo) {return computeArc(arcNo);}, comm);

    // =======================

    // collect system of normal equations
    // ----------------------------------
    if(Parallel::size(comm)>=3)
    {
      logStatus<<"collect system of normal equations"<<Log::endl;
      for(UInt idxInterval=0; idxInterval<timesInterval.size()-1; idxInterval++)
      {
        UInt color = MAX_UINT;
        if(N.at(idxInterval).size())
          color = idxInterval;

        Parallel::CommunicatorPtr commSplit= Parallel::splitCommunicator(color, Parallel::myRank(comm), comm);
        if(commSplit && (Parallel::size(commSplit)>1))
        {
          Parallel::reduceSum(N.at(idxInterval),        0, commSplit);
          Parallel::reduceSum(n.at(idxInterval),        0, commSplit);
          Parallel::reduceSum(lPl.at(idxInterval),      0, commSplit);
          Parallel::reduceSum(obsCount.at(idxInterval), 0, commSplit);
          if(Parallel::myRank(commSplit) != 0)
          {
            N.at(idxInterval) = Matrix();
            n.at(idxInterval) = Matrix();
            lPl.at(idxInterval) = Vector();
            obsCount.at(idxInterval) = 0;
          }
        }
      } // for(idxInterval)
    } // if(Parallel::size()>=3)

    // =======================

    logStatus<<"computing normals"<<Log::endl;
    for(UInt idxInterval=0; idxInterval<timesInterval.size()-1; idxInterval++)
      if(N.at(idxInterval).size())
      {
        NormalEquationInfo normalInfo;
        observation->parameterName(normalInfo.parameterName);

        // eliminate state parameters
        if(eliminatePararmeter && (observation->parameterCount() > observation->gravityParameterCount()))
        {
          // regularize not used parameters
          for(UInt i=observation->gravityParameterCount(); i<N.at(idxInterval).rows(); i++)
            if(N.at(idxInterval)(i,i)==0)
            {
              N.at(idxInterval)(i,i)   += 1.0;
              obsCount.at(idxInterval) += 1;
            }

          Double regul = 1;
          for(;;)
          {
            try
            {
              UInt startElim = observation->gravityParameterCount();
              UInt countElim = observation->parameterCount() - observation->gravityParameterCount();

              // solve for state parameters and update correlation block
              Matrix N22 = N.at(idxInterval).slice(startElim, startElim, countElim, countElim);
              cholesky(N22); // N22 = W^T*W
              triangularSolve(1.0, N22.trans(), N.at(idxInterval).slice(0, startElim, startElim, countElim).trans()); //N12*W^-1

              // update coefficient matrix
              rankKUpdate(-1.0, N.at(idxInterval).slice(0, startElim, startElim, countElim).trans(),
                                N.at(idxInterval).slice(0, 0, startElim, startElim)); // N = N11 - N12 (W^T W)^-1 N12^T

              // update right hand side
              triangularSolve(1.0, N22.trans(), n.at(idxInterval).row(startElim, countElim)); // W^-T * n2
              matMult(-1.0, N.at(idxInterval).slice(0, startElim, startElim, countElim),
                            n.at(idxInterval).row(startElim, countElim),
                            n.at(idxInterval).row(0, startElim)); // n = n1 - N12*W^-1 * W^-T*n2

              // update normals
              obsCount.at(idxInterval) -= countElim;
              for(UInt i=0; i<lPl.at(idxInterval).rows(); i++)
                lPl.at(idxInterval)(i) -= quadsum(n.at(idxInterval).slice(startElim, i, countElim, 1)); // lPl = lPl - n2^T N2^(-1) n2

              N.at(idxInterval) =  N.at(idxInterval).slice(0, 0, startElim, startElim);
              n.at(idxInterval) =  n.at(idxInterval).row(0, startElim);
              break;
            }
            catch(std::exception &e)
            {
              logError<<"error at "<<timesInterval.at(idxInterval).dateStr()<<": "<<e.what()<<" continue..."<<Log::endl;
              for(UInt i=observation->gravityParameterCount(); i<N.at(idxInterval).rows(); i++)
                N.at(idxInterval)(i,i) += regul;
              regul *= 10;
            }
          }

          normalInfo.parameterName.resize(observation->gravityParameterCount());
        }

        // save normals
        // ------------
        normalInfo.lPl = lPl.at(idxInterval);
        normalInfo.observationCount = obsCount.at(idxInterval);
        evaluateTimeVariables(idxInterval, timesInterval.at(idxInterval), timesInterval.at(idxInterval+1), fileNameVariableList);
        logStatus<<"write normal equations to <"<<fileNameNormals(fileNameVariableList)<<">"<<Log::endl;
        writeFileNormalEquation(fileNameNormals(fileNameVariableList), normalInfo, N.at(idxInterval), n.at(idxInterval));
      } // for(idxInterval)
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void KalmanBuildNormals::computeArc(UInt arcNo)
{
  try
  {
    // observation equations
    // ---------------------
    Matrix l,A,B;
    observation->observation(arcNo, l,A,B);
    if(l.rows()==0)
      return;

    // if equations are orthogonaly transformed
    // additional residuals are appended to l
    // ----------------------------------------
    Matrix l2;
    if(l.rows()>A.rows())
    {
      l2 = l.row(A.rows(), l.rows()-A.rows());
      l  = l.row(0, A.rows());
    }

    // eliminate arc related parameters
    // --------------------------------
    if(B.size()!=0)
      eliminationParameter(B,A,l);

    // search time interval
    // --------------------
    UInt idxInterval = 0;
    while(arcsInterval.at(idxInterval+1)<=arcNo)
      idxInterval++;

    // accumulate normal equation system
    // ---------------------------------
    obsCount.at(idxInterval) += l.rows() + l2.rows();
    if(lPl.at(idxInterval).size() == 0)
      lPl.at(idxInterval) = Vector(l.columns());
    for(UInt i=0; i<l.columns(); i++)
      lPl.at(idxInterval)(i) += quadsum(l.column(i)) + quadsum(l2.column(i));
    // right hand side
    if(n.at(idxInterval).size() == 0)
      n.at(idxInterval) = Matrix(A.columns(), l.columns());
    matMult(1., A.trans(), l, n.at(idxInterval));
    // normal matrix
    if(N.at(idxInterval).size() == 0)
      N.at(idxInterval) = Matrix(A.columns(), Matrix::SYMMETRIC);
    rankKUpdate(1., A, N.at(idxInterval));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
