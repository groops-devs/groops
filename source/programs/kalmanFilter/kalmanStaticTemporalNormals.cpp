/***********************************************/
/**
* @file kalmanStaticTemporalNormals.cpp
*
* @brief Combined normal equations with static, temporal and high-frequency gravity field parameters.
*
* @author Torsten Mayer-Guerr
* @date 2012-08-13
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program sets up normal equations based on \configClass{observation}{observationType} for static gravity,
(possibly) multiple long-term \config{temporal} representations and co-estimated high-frequency variations.
It computes the normal equations based on the intervals $i \in \{1, ..., N\}$ given in the \configFile{arcList}{arcList}.
This corresponds to the least-squares adjustment
\begin{equation}
  \label{eq:temporal-least-squares}
  \begin{bmatrix}
    \mathbf{l}_1 \\
    \mathbf{l}_2 \\
    \vdots \\
    \mathbf{l}_N \\
  \end{bmatrix}
  =
  \begin{bmatrix}
    \Phi^{(0)}_1 \mathbf{A}_1  & \Phi^{(1)}_1 \mathbf{A}_1 & \cdots & \Phi^{(m)}_1 \mathbf{A}_1 \\
    \Phi^{(0)}_2 \mathbf{A}_2  & \Phi^{(1)}_2 \mathbf{A}_2 & \cdots & \Phi^{(m)}_2 \mathbf{A}_2 \\
    \vdots  & \vdots &  & \vdots \\
    \Phi^{(0)}_N \mathbf{A}_N  & \Phi^{(1)}_N \mathbf{A}_N & \cdots & \Phi^{(m)}_N \mathbf{A}_N \\
  \end{bmatrix}
  \begin{bmatrix}
    \mathbf{x}^{(0)} \\
    \mathbf{x}^{(1)} \\
    \vdots \\
    \mathbf{x}^{(m)} \\
  \end{bmatrix}
  +
  \begin{bmatrix}
    \mathbf{e}_1 \\
    \mathbf{e}_2 \\
    \vdots \\
    \mathbf{e}_N \\
  \end{bmatrix},
\end{equation}
where $\Phi^{(k)}_i$ ($\Phi^{(0)}_i = 1$) are the temporal basis functions specified by a
\configClass{parametrizationTemporal}{parametrizationTemporalType}, for example, a linear trend $\Phi_i = (t_i - t_0)$.
The normal equations corresponding to \eqref{eq:temporal-least-squares} are
computed by accumulation of each interval which leads to
\begin{equation}
  \label{eq:temporal-normals-full}
  \mathbf{N} = \sum^N_{i=1}
  \begin{bmatrix}
    \Phi^{(0)}_i \mathbf{A}^T_i\\
    \Phi^{(1)}_i \mathbf{A}^T_i \\
    \vdots \\
    \Phi^{(m)}_i \mathbf{A}^T_i \\
  \end{bmatrix}
  \begin{bmatrix}
    \Phi^{(0)}_i \mathbf{A}_i & \Phi^{(1)}_i \mathbf{A}_i & \cdots & \Phi^{(m)}_i \mathbf{A}_i \\
  \end{bmatrix},
  \hspace{15pt}
  \mathbf{n} = \sum^N_{i=1}
  \begin{bmatrix}
    \Phi^{(0)}_i \mathbf{A}^T_i\\
    \Phi^{(1)}_i \mathbf{A}^T_i \\
    \vdots \\
    \Phi^{(m)}_i \mathbf{A}^T_i \\
  \end{bmatrix}
  \mathbf{l}_i.
\end{equation}
For a single  normal equation block $\mathbf{N}^{(k, l)}$ which corresponds to the correlations
between the parameter sets $\mathbf{x}^{(k)}$ and $\mathbf{x}^{(l)}$ follows
\begin{equation}
  \label{eq:temporal-normals-expanded}
  \mathbf{N}^{(k, l)} = \sum_{i=1}^{N} (\Phi^{(k)}_i \mathbf{A}^T_i ) (  \Phi^{(l)}_i  \mathbf{A}_i) =
  \sum_{i=1}^{N} \Phi^{(k)}_i \Phi^{(l)}_i  \mathbf{N}_i,
\end{equation}
where $\mathbf{N}_i = \mathbf{A}_i^T \mathbf{A}_i$.
In analogy to \eqref{eq:temporal-normals-expanded} the right-hand side for each
parameter vector can be computed via
\begin{equation}
  \label{eq:temporal-normals-rhs-expanded}
  \mathbf{n}^{(k)} = \sum_{i=1}^{N} \Phi^{(k)}_i \mathbf{A}^T_i \mathbf{l}_i =  \sum_{i=1}^{N}  \Phi^{(k)}_i \mathbf{n}_i,
\end{equation}
utilizing $\mathbf{n}_i$.
The normal equations are therefore only computed once for each interval and then multiplied with the pairwise product of the temporal factors.
Note that the temporal constituents do not necessarily need to be of the same size as the static part.
The number of parameters of the temporal constituents can be varied with \config{parameterCount}.
This simply slices the first \config{parameterCount} columns from the design matrix $\mathbf{A}_i$
and therefore has to be less or equal than the static parameter count.

If an \configClass{autoregressiveModelSequence}{autoregressiveModelSequenceType} is specified,
high-frequency variations are also set up.
Their purpose is to mitigate temporal aliasing by accounting for short-term gravity field variations.
The parameter count of these high-frequency variations is taken from the dimension of the
\configClass{autoregressiveModelSequence}{autoregressiveModelSequenceType}
and in the same fashion as the temporal variations simply sliced from $\mathbf{A}_i$.
They are modeled as a constant value per interval so the corresponding temporal factor can be expressed as
\begin{equation}
  \label{eq:daily-parameters}
  \Phi_i(t)
  =
  \begin{cases}
    1 &\text{if} \hspace{5pt} t \in [t_i, t_{i+1}) \\
    0 & \text{otherwise}
  \end{cases}.
\end{equation}
A detailed description of the approach is given in:
Kvas, A., Mayer-Gürr, T. GRACE gravity field recovery with background model uncertainties.
J Geod 93, 2543–2552 (2019). \url{https://doi.org/10.1007/s00190-019-01314-1}.

Before writing the normal equations to \configFile{outputfileNormalEquation}{normalEquation}
high-frequency gravity and satellite specific parameters are eliminated.
)";

/***********************************************/

#include "programs/program.h"
#include "parallel/matrixDistributed.h"
#include "files/fileArcList.h"
#include "files/fileNormalEquation.h"
#include "classes/observation/observation.h"
#include "classes/parametrizationTemporal/parametrizationTemporal.h"
#include "misc/kalmanProcessing.h"

/***** CLASS ***********************************/

/** @brief Combined normal equations with static, temporal and high-frequency gravity field parameters.
* @ingroup programsGroup */
class KalmanStaticTemporalNormals
{
public:
  void run(Config &config);

private:
  ObservationPtr    observation;
  std::vector<UInt> arcsInterval;
  std::vector<Time> timesInterval;

  // normal equations
  // ----------------
  MatrixDistributed normals;   // =A'PA, Normal matrix
  Matrix            rhs;       // =A'Pl  right hand side
  Vector            lPl;       // =l'Pl, weighted norm of the observations
  UInt              obsCount;  // number of observations

  // collected normals
  // -----------------
  Bool                computeRightHandSide;
  UInt                turnNo;
  std::vector<UInt>   turnIdx;
  std::vector<UInt>   blockIdxRow, blockIdxCol;
  std::vector<UInt>   idxA;             // start column in design matrix A for each normals block
  std::vector< std::vector<Matrix> > N; // block normals (for each day)
  std::vector<Vector> n;                // =A'Pl, right hand side (for each day)

  void setNormalMatrix(Double factor, Matrix N, UInt rankSource, UInt idxRow, UInt idxCol);
  void reduceIntervalNormals(const std::vector<Double> &factor, const std::vector<Matrix> &N, Bool trans, UInt idxRow, UInt idxCol);
  void computeArc(UInt arcNo);
};

GROOPS_REGISTER_PROGRAM(KalmanStaticTemporalNormals, PARALLEL, "Combined normal equations with static, temporal and high-frequency gravity field parameters.", KalmanFilter, NormalEquation)

/***********************************************/
/***********************************************/

void KalmanStaticTemporalNormals::run(Config &config)
{
  try
  {
    FileName                       fileNameNormals;
    AutoregressiveModelSequencePtr arSequence;
    ParametrizationTemporalPtr     temporal;
    UInt                           countPerTemporalParameter=0;
    FileName                       fileNameArcList;
    UInt                           defaultBlockSize;
    Double                         memSize;

    renameDeprecatedConfig(config, "outputfileNormalequation", "outputfileNormalEquation", date2time(2020, 6, 3));
    renameDeprecatedConfig(config, "arcList",                  "inputfileArcList",         date2time(2020, 7, 7));

    readConfig(config, "outputfileNormalEquation",    fileNameNormals,  Config::MUSTSET,  "", "outputfile for normal equations");
    readConfig(config, "observation",                 observation,      Config::MUSTSET,  "", "");
    readConfig(config, "autoregressiveModelSequence", arSequence,       Config::OPTIONAL, "", "AR model sequence for constraining high frequency variations");
    if(readConfigSequence(config, "temporal", Config::OPTIONAL, "", ""))
    {
      renameDeprecatedConfig(config, "temporalRepresentation", "parametrization", date2time(2020, 6, 3));
      renameDeprecatedConfig(config, "temporalParameterCount", "parameterCount",  date2time(2020, 6, 3));

      readConfig(config, "parametrization", temporal,                  Config::MUSTSET, "", "use state normals to parametrize time variations (trend, annual, ...)");
      readConfig(config, "parameterCount",  countPerTemporalParameter, Config::MUSTSET, "", "use state normals to parametrize time variations (trend, annual, ...)");
      endSequence(config);
    }
    readConfig(config, "inputfileArcList",            fileNameArcList,  Config::OPTIONAL, "",     "list to correspond points of time to arc numbers");
    readConfig(config, "normalsBlockSize",            defaultBlockSize, Config::DEFAULT,  "2048", "block size for distributing the normal equations, 0: one block");
    readConfig(config, "memorySizePerNodeInGByte",    memSize,          Config::DEFAULT,  "0",    "avaiable memory for the normal equations (without observation equations)");
    if(isCreateSchema(config)) return;

    // =======================

    UInt arcCount = observation->arcCount();
    memSize *= 1024*1024*1024; // GByte -> Byte
    if(memSize <= 0) // enough memory avaiable?
      memSize = 1e99;

    // read arc list
    // -------------
    if(!fileNameArcList.empty())
    {
      logStatus<<"read arc list <"<<fileNameArcList<<">"<<Log::endl;
      readFileArcList(fileNameArcList, arcsInterval, timesInterval);
      if(arcsInterval.back()!=arcCount)
        throw(Exception("Different arc count in observation and arcList ("+arcCount%"%i"s+" vs. "+arcsInterval.back()%"%i"s+")."));
    }
    else if(arSequence != nullptr)
      throw(Exception("arcList must be given, if autoregressiveModelSequence is set"));

    // normal equations of high-frequency variations
    // ---------------------------------------------
    UInt countHighFrequencyParameters = 0;
    std::vector< std::vector< std::vector<Matrix> > > normalsHighFrequency;
    if(arSequence != nullptr)
    {
      logStatus<<"initialize normal equations for high-frequency gravity field parameters"<<Log::endl;
      countHighFrequencyParameters = arSequence->dimension();
      normalsHighFrequency = arSequence->normalEquationSequence();
    }

    // temporal variations (trend, annual, ...)
    // ----------------------------------------
    std::vector< std::vector<Double> > factorTemporal;
    UInt temporalParameterCount = 0;
    if(temporal)
    {
      temporalParameterCount = temporal->parameterCount();
      factorTemporal.resize(temporalParameterCount);
      for(UInt k=0; k<factorTemporal.size(); k++)
        factorTemporal.at(k).resize(timesInterval.size()-1);
      for(UInt idxInterval=0; idxInterval<timesInterval.size()-1; idxInterval++)
      {
        const Vector f = temporal->factors(0.5*(timesInterval.at(idxInterval)+timesInterval.at(idxInterval+1)));
        for(UInt k=0; k<factorTemporal.size(); k++)
          factorTemporal.at(k).at(idxInterval) = f(k);
      }
    }

    // test parameters
    // ---------------
    if(countHighFrequencyParameters > observation->gravityParameterCount())
      throw(Exception("More high-frequency gravity field parameters than static gravity field parameters"));
    if(temporalParameterCount)
    {
      if(countPerTemporalParameter==0)
        throw(Exception("No parameters are set for temporal representation"));
      if(countPerTemporalParameter > observation->gravityParameterCount())
        throw(Exception("More temporal parameters than static gravity field parameters"));
      if(countPerTemporalParameter < countHighFrequencyParameters)
        throw(Exception("More process dynamic parameters than temporal gravity field parameters"));
    }
    else
      countPerTemporalParameter = 0;

    // =======================

    // init distributed normal equation system
    // ---------------------------------------
    idxA.clear(); // start column in design matrix A for each normals block
    UInt countParameter = 0;
    std::vector<UInt> blockIndex;
    blockIndex.push_back(0);

    // daily normals
    std::vector<UInt> idxBlockHighFrequency;
    if(countHighFrequencyParameters)
    {
      countParameter += (timesInterval.size()-1) * countHighFrequencyParameters;
      for(UInt i=1; i<timesInterval.size(); i++)
      {
        idxBlockHighFrequency.push_back(blockIndex.size()-1);
        idxA.push_back(0);
        blockIndex.push_back( blockIndex.back() + countHighFrequencyParameters);
      }
    }

    // state parameter (satellite state, calibration parameter, ...)
    const UInt idxBlockState = blockIndex.size()-1;
    countParameter += observation->parameterCount() - observation->gravityParameterCount();
    while(blockIndex.back() < countParameter)
    {
      idxA.push_back( observation->gravityParameterCount() + blockIndex.back() - blockIndex.at(idxBlockState) );
      blockIndex.push_back( (defaultBlockSize==0) ? countParameter : std::min(blockIndex.back()+defaultBlockSize, countParameter));
    }

    // static gravity field
    const UInt idxBlockStatic = blockIndex.size()-1;
    if(countHighFrequencyParameters)
    {
      idxA.push_back( 0 );
      blockIndex.push_back( blockIndex.back()+countHighFrequencyParameters );
    }
    countParameter += countPerTemporalParameter;
    while(blockIndex.back() < countParameter)
    {
      idxA.push_back( blockIndex.back() - blockIndex.at(idxBlockStatic) );
      blockIndex.push_back( (defaultBlockSize==0) ? countParameter : std::min(blockIndex.back()+defaultBlockSize, countParameter));
    }
    countParameter += observation->gravityParameterCount() - countPerTemporalParameter;
    while(blockIndex.back() < countParameter)
    {
      idxA.push_back( blockIndex.back() - blockIndex.at(idxBlockStatic) );
      blockIndex.push_back( (defaultBlockSize==0) ? countParameter : std::min(blockIndex.back()+defaultBlockSize, countParameter));
    }

    // temporal gravity field (annual, trend, ...)
    std::vector<UInt> idxBlockTemporal;
    for(UInt i=0; i<temporalParameterCount; i++)
    {
      idxBlockTemporal.push_back(blockIndex.size()-1);
      if(countHighFrequencyParameters)
      {
        idxA.push_back( 0 );
        blockIndex.push_back( blockIndex.back()+countHighFrequencyParameters );
      }
      countParameter += countPerTemporalParameter;
      while(blockIndex.back() < countParameter)
      {
        idxA.push_back( blockIndex.back() - blockIndex.at(idxBlockTemporal.at(i)) );
        blockIndex.push_back( (defaultBlockSize==0) ? countParameter : std::min(blockIndex.back()+defaultBlockSize, countParameter));
      }
    }
    idxBlockTemporal.push_back(blockIndex.size()-1);

    UInt countBlockTemporal = 0;
    if(temporalParameterCount)
      countBlockTemporal = idxBlockTemporal.at(1)-idxBlockTemporal.at(0);

    // =======================

    // reserve memory
    normals.initEmpty(blockIndex);

    // =======================

    // parameter names
    std::vector<ParameterName> paraNameAll;
    observation->parameterName(paraNameAll);

    std::vector<ParameterName> paraName;
    for(UInt i=idxBlockStatic; i<idxBlockTemporal.at(0); i++)
      for(UInt k=0; k<normals.blockSize(i); k++)
        paraName.push_back( paraNameAll.at(idxA.at(i)+k) );
    if(temporal)
    {
      std::vector<ParameterName> baseName;
      for(UInt k=0; k<countPerTemporalParameter; k++)
        baseName.push_back( paraNameAll.at(idxA.at(idxBlockStatic)+k) );
      temporal->parameterName(baseName, paraName);
    }

    // =======================

    // compute block distribution
    blockIdxRow.clear();
    blockIdxCol.clear();
    for(UInt i=idxBlockState; i<idxBlockTemporal.at(0); i++)
      for(UInt k=i; k<idxBlockTemporal.at(0); k++)
      {
        blockIdxRow.push_back( i );
        blockIdxCol.push_back( k );
      }

    turnIdx.clear();
    turnIdx.push_back(0);
    turnIdx.push_back(0);
    Double memLeft = memSize;
    for(UInt i=0; i<blockIdxRow.size(); i++)
    {
      Double size = sizeof(Double)*normals.blockSize(blockIdxRow.at(i))*normals.blockSize(blockIdxCol.at(i));
      if(size>memLeft)
      {
        memLeft = memSize;
        turnIdx.push_back( turnIdx.back());
      }
      memLeft -= size;
      turnIdx.back()++;
    }

    //=============================================================
    //=============================================================

    if(turnIdx.size()-1 < 2)
      logStatus<<"accumulate normals from observation equations"<<Log::endl;
    else
      logStatus<<"accumulate normals from observation equations in "<<turnIdx.size()-1<<" turns"<<Log::endl;
    // normal equation matrix
    N.resize(blockIdxRow.size());
    for(UInt i=0; i<N.size(); i++)
      N.at(i).resize(timesInterval.size()-1);
    // right hand side
    n.resize(timesInterval.size()-1, Matrix(observation->parameterCount(), observation->rightSideCount()));
    computeRightHandSide = TRUE;
    lPl                  = Vector(observation->rightSideCount());
    obsCount             = 0;

    for(turnNo=0; turnNo<turnIdx.size()-1; turnNo++)
    {
      // compute observation equations
      // -----------------------------
      Parallel::forEachInterval(arcCount, arcsInterval, [this](UInt arcNo) {computeArc(arcNo);});

      logStatus<<"compute right hand sides"<<Log::endl;
      if(computeRightHandSide)
      {
        Parallel::reduceSum(lPl);
        Parallel::reduceSum(obsCount);
        for(UInt idxInterval=0; idxInterval<timesInterval.size()-1; idxInterval++)
          Parallel::reduceSum(n.at(idxInterval));

        if(Parallel::isMaster())
        {
          rhs = Matrix(normals.parameterCount(), lPl.rows());

          // interval
          if(countHighFrequencyParameters)
            for(UInt idxInterval=0; idxInterval<timesInterval.size()-1; idxInterval++)
            {
              axpy(1., n.at(idxInterval).row(idxA.at(idxInterval), normals.blockSize(idxInterval)),
                   rhs.row(normals.blockIndex(idxInterval), normals.blockSize(idxInterval)));
            }

          // state, static
          for(UInt idxBlock=idxBlockState; idxBlock<idxBlockTemporal.at(0); idxBlock++)
            for(UInt idxInterval=0; idxInterval<timesInterval.size()-1; idxInterval++)
            {
              axpy(1., n.at(idxInterval).row(idxA.at(idxBlock), normals.blockSize(idxBlock)),
                   rhs.row(normals.blockIndex(idxBlock), normals.blockSize(idxBlock)));
            }

          // temporal
          for(UInt k=0; k<temporalParameterCount; k++)
            for(UInt idxBlock=idxBlockTemporal.at(k); idxBlock<idxBlockTemporal.at(k+1); idxBlock++)
              for(UInt idxInterval=0; idxInterval<timesInterval.size()-1; idxInterval++)
              {
                axpy(factorTemporal.at(k).at(idxInterval), n.at(idxInterval).row(idxA.at(idxBlock), normals.blockSize(idxBlock)),
                    rhs.row(normals.blockIndex(idxBlock), normals.blockSize(idxBlock)));
              }
        } // if(isMaster())

        n.clear();
        computeRightHandSide = FALSE;
      }

      // ==================================

      // collect system of normal equations
      // ----------------------------------
      if(Parallel::size()>=3)
      {
        logStatus<<"collect system of normal equations"<<Log::endl;
        // collect normals separated in each interval
        for(UInt idxInterval=0; idxInterval<timesInterval.size()-1; idxInterval++)
        {
          UInt color = NULLINDEX;
          if(N.at(turnIdx.at(turnNo)).at(idxInterval).size())
            color = 1;
          Parallel::CommunicatorPtr comm = Parallel::splitCommunicator(color, Parallel::myRank());
          if((comm!=nullptr) && (Parallel::size(comm)>1))
          {
            for(UInt i=turnIdx.at(turnNo); i<turnIdx.at(turnNo+1); i++)
            {
              Parallel::reduceSum(N.at(i).at(idxInterval), 0, comm);
              if(Parallel::myRank(comm) != 0)
                N.at(i).at(idxInterval) = Matrix();
            } // for(i=turnIdx)
          } // if(comm!=nullptr)
        }
      }

      // rank of the reduced normals
      // ---------------------------
      logStatus<<"rank of the reduced normals"<<Log::endl;
      std::vector<UInt> rankInterval(timesInterval.size()-1, 0);
      for(UInt idxInterval=0; idxInterval<timesInterval.size()-1; idxInterval++)
      {
        if(N.at(turnIdx.at(turnNo)).at(idxInterval).size())
          rankInterval.at(idxInterval) = Parallel::myRank();
        Parallel::reduceSum(rankInterval.at(idxInterval));
        Parallel::broadCast(rankInterval.at(idxInterval));
      }

      // ==================================

      // setup normal equations
      // ----------------------
      logStatus<<"setup normal equations"<<Log::endl;
      logTimerStart;
      for(UInt i=turnIdx.at(turnNo); i<turnIdx.at(turnNo+1); i++)
      {
        logTimerLoop(i-turnIdx.at(turnNo), turnIdx.at(turnNo+1)-turnIdx.at(turnNo));

        for(UInt idxInterval=0; idxInterval<timesInterval.size()-1; idxInterval++)
          if(Parallel::myRank() == rankInterval.at(idxInterval))
          {
            if(N.at(i).at(idxInterval).size() == 0)
              N.at(i).at(idxInterval) = Matrix(normals.blockSize(blockIdxRow.at(i)), normals.blockSize(blockIdxCol.at(i)));
            if(N.at(i).at(idxInterval).getType() == Matrix::SYMMETRIC)
              fillSymmetric(N.at(i).at(idxInterval));
            N.at(i).at(idxInterval).setType(Matrix::GENERAL);
          }

        if(countHighFrequencyParameters)
        {
          // interval
          // --------
          if((blockIdxRow.at(i) == idxBlockStatic) && (blockIdxCol.at(i) == idxBlockStatic))
            for(UInt idxInterval=0; idxInterval<timesInterval.size()-1; idxInterval++)
              setNormalMatrix(1., N.at(i).at(idxInterval), rankInterval.at(idxInterval), idxInterval, idxInterval);

          // interval <-> state
          // ------------------
          if((blockIdxRow.at(i) >= idxBlockState) && (blockIdxRow.at(i) < idxBlockStatic) && (blockIdxCol.at(i) == idxBlockStatic))
            for(UInt idxInterval=0; idxInterval<timesInterval.size()-1; idxInterval++)
              setNormalMatrix(1., N.at(i).at(idxInterval).trans(), rankInterval.at(idxInterval), idxInterval, blockIdxRow.at(i));

          // interval <-> static
          // -------------------
          if(blockIdxRow.at(i) == idxBlockStatic)
            for(UInt idxInterval=0; idxInterval<timesInterval.size()-1; idxInterval++)
              setNormalMatrix(1., N.at(i).at(idxInterval), rankInterval.at(idxInterval), idxInterval, blockIdxCol.at(i));

          // interval <-> temporal
          // ---------------------
          if((blockIdxRow.at(i) == idxBlockStatic) &&
             (blockIdxCol.at(i) >= idxBlockStatic) && (blockIdxCol.at(i) < idxBlockStatic + countBlockTemporal))
          {
            for(UInt idxInterval=0; idxInterval<timesInterval.size()-1; idxInterval++)
              for(UInt k=0; k<temporalParameterCount; k++)
                setNormalMatrix(factorTemporal.at(k).at(idxInterval), N.at(i).at(idxInterval), rankInterval.at(idxInterval), idxInterval, blockIdxCol.at(i) - idxBlockStatic + idxBlockTemporal.at(k));
          }
        } // if(countHighFrequencyParameter)

        // state, static
        // -------------
        std::vector<Double> factor(timesInterval.size()-1, 1.);
        reduceIntervalNormals(factor, N.at(i), FALSE, blockIdxRow.at(i), blockIdxCol.at(i));

        // state,static <-> temporal
        // -------------------------
        if((blockIdxRow.at(i) >= idxBlockState)  && (blockIdxRow.at(i) < idxBlockStatic + countBlockTemporal) &&
           (blockIdxCol.at(i) >= idxBlockStatic) && (blockIdxCol.at(i) < idxBlockStatic + countBlockTemporal))
        {
          for(UInt k=0; k<temporalParameterCount; k++)
            reduceIntervalNormals(factorTemporal.at(k), N.at(i), FALSE, blockIdxRow.at(i), blockIdxCol.at(i) - idxBlockStatic + idxBlockTemporal.at(k));
        }

        // static <-> temporal
        // -------------------
        if((blockIdxRow.at(i) >= idxBlockStatic) && (blockIdxRow.at(i) < idxBlockStatic + countBlockTemporal) &&
           (blockIdxCol.at(i) >= idxBlockStatic + countBlockTemporal))
        {
          for(UInt k=0; k<temporalParameterCount; k++)
            reduceIntervalNormals(factorTemporal.at(k), N.at(i), TRUE, blockIdxCol.at(i), blockIdxRow.at(i)-idxBlockStatic+idxBlockTemporal.at(k));
        }

        // temporal <-> temporal
        // ---------------------
        if((blockIdxRow.at(i) >= idxBlockStatic) && (blockIdxRow.at(i) < idxBlockStatic + countBlockTemporal) &&
           (blockIdxCol.at(i) >= idxBlockStatic) && (blockIdxCol.at(i) < idxBlockStatic + countBlockTemporal))
        {
          for(UInt k=0; k<temporalParameterCount; k++)
            for(UInt l=k; l<temporalParameterCount; l++)
            {
              std::vector<Double> factor(timesInterval.size()-1);
              std::transform(factorTemporal.at(k).begin(), factorTemporal.at(k).end(), factorTemporal.at(l).begin(), factor.begin(), std::multiplies<Double>());
              reduceIntervalNormals(factor, N.at(i), FALSE, idxBlockTemporal.at(k) + blockIdxRow.at(i)-idxBlockStatic, idxBlockTemporal.at(l) + blockIdxCol.at(i)-idxBlockStatic);
            }
        } // if(block is temporal)

        // free memory
        N.at(i).clear();
      } // for(i=turnIdx.at(turnNo))
      logTimerLoopEnd(turnIdx.at(turnNo+1)-turnIdx.at(turnNo));
    } // for(turnNo)

    for(UInt i=0; i<normals.blockCount(); i++)
      if(normals.isMyRank(i,i))
        normals.N(i,i).setType(Matrix::SYMMETRIC);

    observation = ObservationPtr(nullptr);

    // ========================================================

    // fill symmetric
    // --------------
    logStatus<<"fill symmetric"<<Log::endl;
    for(UInt idxRow=0; idxRow<countBlockTemporal; idxRow++)
      for(UInt idxCol=idxRow+1; idxCol<countBlockTemporal; idxCol++)
      {
        for(UInt k=0; k<temporalParameterCount; k++)
        {
          // static <-> temporal
          setNormalMatrix(1., normals.N(idxBlockStatic+idxRow, idxBlockTemporal.at(k)+idxCol).trans(),
                          normals.rank(idxBlockStatic+idxRow, idxBlockTemporal.at(k)+idxCol),
                          idxBlockStatic+idxCol, idxBlockTemporal.at(k)+idxRow);

          // temporal <-> temporal
          for(UInt l=k+1; l<temporalParameterCount; l++)
            setNormalMatrix(1., normals.N(idxBlockTemporal.at(k)+idxRow, idxBlockTemporal.at(l)+idxCol).trans(),
                            normals.rank(idxBlockTemporal.at(k)+idxRow, idxBlockTemporal.at(l)+idxCol),
                            idxBlockTemporal.at(k)+idxCol, idxBlockTemporal.at(l)+idxRow);
        }
      }

    // ========================================================

    if(countHighFrequencyParameters)
    {
      logStatus<<"normals of process dynamic"<<Log::endl;
      if(Parallel::isMaster())
        obsCount += idxBlockHighFrequency.size() * countHighFrequencyParameters;

      for(UInt id=0; id<idxBlockHighFrequency.size(); id++)
      {
        const UInt idx = std::min(id, normalsHighFrequency.size()-1);
        for(UInt i=0; i<normalsHighFrequency.at(idx).size(); i++)
          for(UInt k=i; k<normalsHighFrequency.at(idx).at(i).size(); k++)
          {
            normals.setBlock(idxBlockHighFrequency.at(id+i-idx), idxBlockHighFrequency.at(id+k-idx));
            if(normals.isMyRank(idxBlockHighFrequency.at(id+i-idx), idxBlockHighFrequency.at(id+k-idx)))
              axpy(1.0, normalsHighFrequency.at(idx).at(i).at(k), normals.N(idxBlockHighFrequency.at(id+i-idx), idxBlockHighFrequency.at(id+k-idx)));
          }
      }
    }

    // Regularize not used parameters
    // ------------------------------
    UInt countRegul = 0;
    for(UInt i=0; i<idxBlockStatic; i++)
    {
      normals.setBlock(i,i);
      if(normals.isMyRank(i,i))
      {
        Matrix &N = normals.N(i,i);
        for(UInt k=0; k<normals.blockSize(i); k++)
          if(N(k,k) == 0)
          {
            N(k,k) = 1.;
            countRegul++;
          }
      }
    }
    Parallel::reduceSum(countRegul);
    if(Parallel::isMaster() && countRegul)
      logWarning<<countRegul<<" parameters are not used"<<Log::endl;

    // eliminate interval & state parameters
    // -------------------------------------
    if(idxBlockStatic>0)
    {
      logStatus<<"eliminate interval parameters from normal equations"<<Log::endl;
      normals.cholesky(TRUE/*timing*/, 0, idxBlockStatic, TRUE/*collect*/);
      normals.triangularTransSolve(rhs, 0, idxBlockStatic);
      if(Parallel::isMaster())
      {
        obsCount -= normals.blockIndex(idxBlockStatic);
        // lPl = lPl - n2^T N2^(-1) n2
        for(UInt i=0; i<lPl.rows(); i++)
          lPl(i) -= quadsum(rhs.slice(0,i,normals.blockIndex(idxBlockStatic),1));
        // remove additional parameters
        rhs = rhs.row(normals.blockIndex(idxBlockStatic), normals.parameterCount() - normals.blockIndex(idxBlockStatic));
      }
      normals.eraseBlocks(0, idxBlockStatic);
    }

    // Write normal equations
    // ----------------------
    logStatus<<"write normal equations to <"<<fileNameNormals<<">"<<Log::endl;
    writeFileNormalEquation(fileNameNormals, NormalEquationInfo(paraName, lPl, obsCount), normals, rhs);
    Parallel::barrier();
    logStatus<<"finish"<<Log::endl;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void KalmanStaticTemporalNormals::setNormalMatrix(Double factor, Matrix N, UInt rankSource, UInt idxRow, UInt idxCol)
{
  try
  {
    normals.setBlock(idxRow, idxCol);

    if(rankSource != normals.rank(idxRow, idxCol))
    {
      if(Parallel::myRank() == rankSource)
        Parallel::send(N, normals.rank(idxRow, idxCol));
      if(Parallel::myRank() == normals.rank(idxRow, idxCol))
        Parallel::receive(N, rankSource);
    }

    if(normals.isMyRank(idxRow, idxCol))
    {
      normals.N(idxRow, idxCol) = N;
      if(factor!=1)
        normals.N(idxRow, idxCol) *= factor;
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void KalmanStaticTemporalNormals::reduceIntervalNormals(const std::vector<Double> &factor, const std::vector<Matrix> &N, Bool trans, UInt idxRow, UInt idxCol)
{
  try
  {
    normals.setBlock(idxRow, idxCol);
    for(UInt idxInterval=0; idxInterval<N.size(); idxInterval++)
      if(N.at(idxInterval).size())
      {
        if(normals.N(idxRow, idxCol).size()==0)
          normals.N(idxRow, idxCol) = ((idxRow==idxCol) ? Matrix(normals.blockSize(idxRow), Matrix::SYMMETRIC, Matrix::UPPER)
                                                        : Matrix(normals.blockSize(idxRow), normals.blockSize(idxCol)));
        axpy(factor.at(idxInterval), ((trans) ? N.at(idxInterval).trans() : N.at(idxInterval)), normals.N(idxRow, idxCol));
      }

    // reduce normals
    UInt color = ((normals.N(idxRow, idxCol).size()) ? 1 : NULLINDEX);
    UInt key   = Parallel::myRank()+1;
    if(normals.isMyRank(idxRow, idxCol))
    {
      color = 1;
      key   = 0;
    }
    Parallel::CommunicatorPtr comm = Parallel::splitCommunicator(color, key);
    if((comm != nullptr) && (Parallel::size(comm)>1))
      Parallel::reduceSum(normals.N(idxRow, idxCol), 0, comm);
    if(!normals.isMyRank(idxRow, idxCol))
      normals.N(idxRow, idxCol) = Matrix();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void KalmanStaticTemporalNormals::computeArc(UInt arcNo)
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
    Matrix l2;
    if(l.rows()>A.rows())
    {
      l2 = l.row(A.rows(), l.rows()-A.rows());
      l  = l.row(0, A.rows());
    }

    // eliminate arc dependent parameters
    // ----------------------------------
    Vector tau;
    if(B.size())
    {
      tau = QR_decomposition(B);
      QTransMult(B, tau, l); // transform observations: l:= Q'l
      QTransMult(B, tau, A); // transform design matrix A:=Q'A
    }
    // use only nullspace of design matrix B
    MatrixSlice A_bar( A.row(B.columns(), A.rows()-B.columns()) );
    MatrixSlice l_bar( l.row(B.columns(), l.rows()-B.columns()) );

    // search time interval
    // --------------------
    UInt idxInterval = 0;
    while(arcsInterval.at(idxInterval+1)<=arcNo)
      idxInterval++;

    // right hand side
    // ---------------
    if(computeRightHandSide)
    {
      obsCount += l_bar.rows() + l2.rows();
      for(UInt i=0; i<l_bar.columns(); i++)
        lPl(i) += quadsum(l_bar.column(i)) + quadsum(l2.column(i));
      if(n.at(idxInterval).size() == 0)
        n.at(idxInterval) = Matrix(A_bar.columns(), l_bar.columns());
      matMult(1., A_bar.trans(), l_bar, n.at(idxInterval));
    }

    // accumulate normals
    // ------------------
    for(UInt i=turnIdx.at(turnNo); i<turnIdx.at(turnNo+1); i++)
    {
      const UInt rows = normals.blockSize(blockIdxRow.at(i));
      const UInt cols = normals.blockSize(blockIdxCol.at(i));

      // N = A'A
      if(blockIdxRow.at(i) == blockIdxCol.at(i))
      {
        if(N.at(i).at(idxInterval).size() == 0)
          N.at(i).at(idxInterval) = Matrix(rows, Matrix::SYMMETRIC, Matrix::UPPER);
        rankKUpdate(1., A_bar.column(idxA.at(blockIdxRow.at(i)), rows), N.at(i).at(idxInterval));
      }
      else
      {
        if(N.at(i).at(idxInterval).size() == 0)
          N.at(i).at(idxInterval) = Matrix(rows, cols);
        matMult(1., A_bar.column(idxA.at(blockIdxRow.at(i)), rows).trans(), A_bar.column(idxA.at(blockIdxCol.at(i)), cols), N.at(i).at(idxInterval));
      }
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
