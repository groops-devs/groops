/***********************************************/
/**
* @file gnssProcessing.cpp
*
* @brief GNSS/LEO satellite orbit determination, station network analysis, PPP.
*
* @author Torsten Mayer-Guerr
* @author Sebastian Strasser
* @date 2010-08-04
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program processes GNSS observations. The primary use cases are:
\begin{itemize}
  \item \reference{GNSS satellite orbit determination and station network analysis}{cookbook.gnssNetwork}
  \item \reference{Kinematic orbit determination of LEO satellites}{cookbook.kinematicOrbit}
  \item \reference{GNSS precise point positioning (PPP)}{cookbook.gnssPpp}
\end{itemize}

One or more GNSS constellations must be defined via \configClass{transmitter}{gnssParametrizationTransmitterType}.
Receivers such as ground station networks or Low Earth Orbit (LEO) satellites can be defined via \configClass{receiver}{gnssParametrizationReceiverType}.

Earth rotation parameters (ERPs) can be estimated via \configClass{earthRotation}{gnssParametrizationEarthRotationType}.

Constraints can be added to the normal equation system via \configClass{constraints}{gnssParametrizationConstraintsType} in case there are no direct constraint options
in the respective config element (e.g. \config{noNetRotationSigma} in \configClass{receiver:stationNetwork}{gnssParametrizationReceiverType:stationNetwork}).
For details on possible additional constraints see \configClass{receiver:stationNetwork}{gnssParametrizationReceiverType:stationNetwork} (tropospheric parameters),
\configClass{transmitter}{gnssParametrizationTransmitterType} (antenna center offsets), or \configClass{earthRotation}{gnssParametrizationEarthRotationType} (dUT1).

The processing flow is controlled via a list of \configClass{processingSteps}{gnssProcessingStepType}.
Each step is processed consecutively. Some steps allow the selection of parameters, epochs, or the normal equation structure,
which affects all subsequent steps.

If the program is run on multiple processes in \reference{parallel}{general.parallelization}, \config{parallelIntervals}=\verb|yes| distributes the
defined \configClass{intervals}{timeSeriesType} over all processes. Otherwise the intervals are processed consecutively while the
\configClass{receiver}{gnssParametrizationReceiverType}s (stations or LEO satellites) are distributed over the processes.

See also \program{GnssSimulateReceiver}.
)";

/***********************************************/

#include "programs/program.h"
#include "config/configRegister.h"
#include "base/planets.h"
#include "parallel/matrixDistributed.h"
#include "parser/dataVariables.h"
#include "files/fileStringTable.h"
#include "files/fileParameterName.h"
#include "files/fileNormalEquation.h"
#include "files/fileMatrix.h"
#include "misc/varianceComponentEstimation.h"
#include "gnss/gnss.h"
#include "gnss/gnssObservation.h"
#include "gnss/gnssDesignMatrix.h"
#include "gnss/gnssParametrizationTransmitter.h"
#include "gnss/gnssParametrizationReceiver.h"
#include "gnss/gnssParametrizationIonosphere.h"
#include "gnss/gnssParametrizationEarthRotation.h"
#include "gnss/gnssParametrizationGravityField.h"
#include "gnss/gnssParametrizationAmbiguities.h"
#include "gnss/gnssParametrizationConstraints.h"
#include "classes/timeSeries/timeSeries.h"
#include "classes/parameterSelector/parameterSelector.h"

/***** TYPES ***********************************/

class GnssProcessingStep;
typedef std::shared_ptr<GnssProcessingStep> GnssProcessingStepPtr;

/***** CLASS ***********************************/

/** @brief GNSS/LEO satellite orbit determination, station network analysis, PPP.
* @ingroup programsGroup */
class GnssProcessing
{
public:
  FileName           fileNameNormalsInfo;
  std::vector<Time>  timeSeries, timesInterval;
  VariableList       fileNameVariableList;
  Double             marginSeconds;
  Bool               parallelIntervals;

  // Gnss
  // ----
  Gnss gnss;
  GnssParametrizationIonospherePtr gnssIonosphere;

  // System of normal equations
  // --------------------------
  MatrixDistributed   normals;
  std::vector<Matrix> n;        // at master (after solve)
  Vector              lPl;      // at master (after solve)
  UInt                obsCount; // at master (after solve)

  // processing steps
  // ----------------
  std::vector<GnssProcessingStepPtr> processingSteps;

  void   computeInterval(UInt idInterval, Parallel::CommunicatorPtr comm);
  void   initParameter(Gnss::NormalEquationInfo &normalEquationInfo);
  void   regularizeNotUsedParameters(const Gnss::NormalEquationInfo &normalEquationInfo, UInt blockStart, UInt blockCount);
  void   collectNormalsBlocks(UInt blockStart, UInt blockCount);
  void   buildNormals(Gnss::NormalEquationInfo &normalEquationInfo, Bool constraintsOnly, Bool solveEpochParameters, Bool removeEpochParameters);
  Double estimateSolution(Gnss::NormalEquationInfo &normalEquationInfo,
                          Bool resolveAmbiguities, Bool dryRun, Bool computeResiduals,
                          Gnss::Receiver::WeightingType computeWeights, Gnss::Receiver::WeightingType adjustSigma0);

  void run(Config &config);
};

GROOPS_REGISTER_PROGRAM(GnssProcessing, PARALLEL, "GNSS/LEO satellite orbit determination, station network analysis, PPP", Gnss)

/***********************************************/

template<> Bool readConfig(Config &config, const std::string &name, Gnss::Receiver::WeightingType &type, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation)
{
  try
  {
    std::string choice;
    if(!readConfigChoice(config, name, choice, mustSet, defaultValue, annotation))
      return FALSE;
    if(readConfigChoiceElement(config, "individual",    choice, "individual observations (e.g. L1P, L2P, C1P, C2P)"))
      type = Gnss::Receiver::WeightingType::INDIVIDUAL;
    if(readConfigChoiceElement(config, "grouped",       choice, "grouped by phase and code observation (e.g. L**, C**)"))
      type = Gnss::Receiver::WeightingType::GROUPED;
    if(readConfigChoiceElement(config, "groupedPhases", choice, "grouped phases and individual code (e.g. L**, C1P, C2P)"))
      type = Gnss::Receiver::WeightingType::GROUPEDPHASES;
    endChoice(config);
    return TRUE;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***** CLASS ***********************************/
/***********************************************/

/**
* @defgroup gnssProcessingStepGroup GnssProcessingStep
* @brief Processing steps in @ref GnssProcessing.
* @ingroup gnssGroup
* The interface is given by @ref GnssProcessingStep.
* An Instance can be created by @ref readConfig. */
/// @{

static const char *docstringGnssProcessingStep = R"(
\section{GnssProcessingStep}\label{gnssProcessingStepType}
Processing step in \program{GnssProcessing}.

Processing steps enable a dynamic definition of the consecutive steps performed during any kind of GNSS processing.
The most common steps are \configClass{estimate}{gnssProcessingStepType:estimate}, which performs an iterative least
squares adjustment, and \configClass{writeResults}{gnssProcessingStepType:writeResults}, which writes all output files
defined in \program{GnssProcessing} and is usually the last step.
Some steps such as \configClass{selectParameters}{gnssProcessingStepType:selectParameters},
\configClass{selectEpochs}{gnssProcessingStepType:selectEpochs},
\configClass{selectNormalsBlockStructure}{gnssProcessingStepType:selectNormalsBlockStructure}, and
\configClass{selectReceivers}{gnssProcessingStepType:selectReceivers} affect all subsequent steps.
In case these steps are used within a \configClass{group}{gnssProcessingStepType:group} or
\configClass{forEachReceiverSeparately}{gnssProcessingStepType:forEachReceiverSeparately} step,
they only affect the steps within this level.

For usage examples see cookbooks on \reference{GNSS satellite orbit determination and network analysis}{cookbook.gnssNetwork:processing}
or \reference{Kinematic orbit determination of LEO satellites}{cookbook.kinematicOrbit}.
)";

/** @brief Processing steps in GnssProcessing.
* An Instance of this class can be created by @ref readConfig. */
class GnssProcessingStep
{
protected:
  std::vector<Byte> selectReceivers(const Gnss &gnss, const Vector &estimableReceiver,
                                    const std::vector<std::vector<std::string>> &selectNames,
                                    const std::vector<std::vector<std::string>> &excludeNames) const;
public:
  virtual ~GnssProcessingStep() {}

  /** @brief Execute the processing step. */
  virtual void process(GnssProcessing &gnssProcessing, const Vector &estimableReceiver,
                       Gnss::NormalEquationInfo &normalEquationInfo, Bool &changedNormalEquationInfo) const = 0;

  /** @brief creates an derived instance of this class. */
  static GnssProcessingStepPtr create(Config &config, const std::string &name);
};

/***********************************************/

static const char *docstringGnssProcessingStepEstimate = R"(
\subsection{Estimate}\label{gnssProcessingStepType:estimate}
Iterative least squares adjustment.

Iterates until either \config{convergenceThreshold} or \config{maxIterationCount} are reached.
If \config{computeResiduals}=\verb|yes| and \config{computeWeights} or \config{adjustSigma0} are set,
observation weights and/or $\sigma_0$ are adjusted after each iteration.
)";

/** @brief Iterative least squares adjustment.
* @see GnssProcessingStep */
class GnssProcessingStepEstimate : public GnssProcessingStep
{
  Bool                          computeResiduals;
  Gnss::Receiver::WeightingType adjustSigma0;
  Gnss::Receiver::WeightingType computeWeights;
  Double                        convergenceThreshold;
  UInt                          iterCount;

public:
  GnssProcessingStepEstimate(Config &config)
  {
    computeWeights = Gnss::Receiver::WeightingType::NONE;
    adjustSigma0   = Gnss::Receiver::WeightingType::NONE;

    readConfig(config, "computeResiduals",     computeResiduals,     Config::DEFAULT,  "1",    "");
    readConfig(config, "adjustSigma0",         adjustSigma0,         Config::OPTIONAL, "",     "adjust sigma0 by scale factor (per receiver)");
    readConfig(config, "computeWeights",       computeWeights,       Config::OPTIONAL, "",     "downweight outliers");
    readConfig(config, "convergenceThreshold", convergenceThreshold, Config::DEFAULT,  "0.01", "[m] stop iteration once full convergence is reached");
    readConfig(config, "maxIterationCount",    iterCount,            Config::DEFAULT,  "3",    "maximum number of iterations");
  }

  void process(GnssProcessing &gnssProcessing, const Vector &/*estimableReceiver*/,
               Gnss::NormalEquationInfo &normalEquationInfo, Bool &changedNormalEquationInfo) const override
  {
    try
    {
      logStatus<<"=== estimate ================================================"<<Log::endl;
      if(changedNormalEquationInfo)
        gnssProcessing.initParameter(normalEquationInfo);
      changedNormalEquationInfo = FALSE;
      for(UInt iter=0; iter<iterCount; iter++)
      {
        logStatus<<iter+1<<". iteration  --------------------------"<<Log::endl;
        if(convergenceThreshold > gnssProcessing.estimateSolution(normalEquationInfo, FALSE/*resolveAmbiguities*/, FALSE/*dryRun*/, computeResiduals, computeWeights, adjustSigma0))
          break;
      }
    }
    catch(std::exception &e)
    {
      GROOPS_RETHROW(e)
    }
  }
};

/***********************************************/

static const char *docstringGnssProcessingStepResolveAmbiguities = R"(
\subsection{ResolveAmbiguities}\label{gnssProcessingStepType:resolveAmbiguities}
Least squares adjustment with ambiguity resolution.

After this step all resolved ambiguities are removed from the normal equation system, except if \config{dryRun}=\verb|yes|.
If \config{computeResiduals}=\verb|yes| and \config{computeWeights} or \config{adjustSigma0} are set,
observation weights and/or $\sigma_0$ are adjusted after resolving the ambiguities.

For a description of the ambiguity resolution algorithm see \configClass{GnssParametrizationAmbiguities}{gnssParametrizationAmbiguitiesType}.
)";

/** @brief Least squares adjustment with ambiguity resolution.
* @see GnssProcessingStep */
class GnssProcessingStepResolveAmbiguities : public GnssProcessingStep
{
  Bool                          computeResiduals;
  Gnss::Receiver::WeightingType adjustSigma0;
  Gnss::Receiver::WeightingType computeWeights;
  Bool                          dryRun;

public:
  GnssProcessingStepResolveAmbiguities(Config &config)
  {
    computeWeights = Gnss::Receiver::WeightingType::NONE;
    adjustSigma0   = Gnss::Receiver::WeightingType::NONE;

    readConfig(config, "computeResiduals", computeResiduals, Config::DEFAULT,  "1",    "");
    readConfig(config, "adjustSigma0",     adjustSigma0,     Config::OPTIONAL, "",     "adjust sigma0 by scale factor (per receiver)");
    readConfig(config, "computeWeights",   computeWeights,   Config::OPTIONAL, "",     "downweight outliers");
    readConfig(config, "dryRun",           dryRun,           Config::DEFAULT,  "0",    "compute integer ambiguities but don't mark as solved");
  }

  void process(GnssProcessing &gnssProcessing, const Vector &/*estimableReceiver*/,
               Gnss::NormalEquationInfo &normalEquationInfo, Bool &changedNormalEquationInfo) const override
  {
    try
    {
      logStatus<<"=== resolve ambiguities  ===================================="<<Log::endl;
      if(changedNormalEquationInfo)
        gnssProcessing.initParameter(normalEquationInfo);
      changedNormalEquationInfo = FALSE;
      gnssProcessing.estimateSolution(normalEquationInfo, TRUE/*resolveAmbiguities*/, dryRun, computeResiduals, computeWeights, adjustSigma0);
      gnssProcessing.initParameter(normalEquationInfo);
    }
    catch(std::exception &e)
    {
      GROOPS_RETHROW(e)
    }
  }
};

/***********************************************/

static const char *docstringGnssProcessingStepComputeCovarianceMatrix = R"(
\subsection{ComputeCovarianceMatrix}\label{gnssProcessingStepType:computeCovarianceMatrix}
Accumulate normal equations and compute the sparse inverse.
)";

/** @brief Accumulate normal equations and compute the sparse inverse.
* @see GnssProcessingStep */
class GnssProcessingStepComputeCovarianceMatrix : public GnssProcessingStep
{
public:
  GnssProcessingStepComputeCovarianceMatrix(Config &/*config*/) {}

  void process(GnssProcessing &gnssProcessing, const Vector &/*estimableReceiver*/,
               Gnss::NormalEquationInfo &normalEquationInfo, Bool &changedNormalEquationInfo) const override
  {
    try
    {
      logStatus<<"=== compute covariance matrix ==============================="<<Log::endl;
      if(changedNormalEquationInfo)
        gnssProcessing.initParameter(normalEquationInfo);
      changedNormalEquationInfo = FALSE;
      gnssProcessing.buildNormals(normalEquationInfo, FALSE/*constraintsOnly*/, FALSE/*solveEpochParameters*/, FALSE/*removeEpochParameters*/);
      logStatus<<"cholesky"<<Log::endl;
      gnssProcessing.normals.cholesky(TRUE/*timing*/);
      logStatus<<"sparse inverse"<<Log::endl;
      gnssProcessing.normals.cholesky2SparseInverse();
      gnssProcessing.gnss.updateCovariance(normalEquationInfo, gnssProcessing.normals);
    }
    catch(std::exception &e)
    {
      GROOPS_RETHROW(e)
    }
  }
};

/***********************************************/

static const char *docstringGnssProcessingStepWriteNormalEquations = R"(
\subsection{WriteNormalEquations}\label{gnssProcessingStepType:writeNormalEquations}
Write unconstrained normal equations and constraint equations.

Epoch parameters (e.g. clock errors) are eliminated before writing the normal equations
if \config{eliminateEpochParameters}=\verb|yes|. If \configClass{remainingParameters}{parameterSelectorType}
is set only the selected parameters are written to the normal equations and all other parameters are eliminated beforehand.
)";

/** @brief Write unconstrained normal equations and constraint equations.
* @see GnssProcessingStep */
class GnssProcessingStepWriteNormalEquations : public GnssProcessingStep
{
  FileName             fileNameNormals, fileNameConstraints, fileNameApriori, fileNameParameterNames;
  std::string          variableReceiver;
  Bool                 eliminateEpochParameters;
  ParameterSelectorPtr parameterSelector;
  UInt                 defaultBlockSize;

public:
  GnssProcessingStepWriteNormalEquations(Config &config)
  {
    defaultBlockSize = NULLINDEX;
    readConfig(config, "outputfileNormalEquations",     fileNameNormals,          Config::OPTIONAL, "output/normals_{loopTime:%D}.dat",        "unconstrained normals");
    readConfig(config, "outputfileConstraintEquations", fileNameConstraints,      Config::OPTIONAL, "output/constraints_{loopTime:%D}.dat",    "applied constraints");
    readConfig(config, "outputfileAprioriSolution",     fileNameApriori,          Config::OPTIONAL, "output/x0_{loopTime:%D}.txt",             "a priori parameters");
    readConfig(config, "outputfileParameterName",       fileNameParameterNames,   Config::OPTIONAL, "output/parameterNames_{loopTime:%D}.txt", "parameter names");
    readConfig(config, "variableReceiver",              variableReceiver,         Config::OPTIONAL, "station", "variable in single receiver loop");
    readConfig(config, "eliminateEpochParameters",      eliminateEpochParameters, Config::DEFAULT,  "1",    "");
    readConfig(config, "remainingParameters",           parameterSelector,        Config::OPTIONAL, "",     "parameter order/selection of output normal equations");
    readConfig(config, "defaultNormalsBlockSize",       defaultBlockSize,         Config::OPTIONAL, "2048", "block size for distributing the normal equations, 0: one block, empty: original block size");
  }

  void process(GnssProcessing &gnssProcessing, const Vector &/*estimableReceiver*/,
               Gnss::NormalEquationInfo &normalEquationInfo, Bool &changedNormalEquationInfo) const override
  {
    try
    {
      logStatus<<"=== write normal equations =================================="<<Log::endl;
      if(changedNormalEquationInfo)
        gnssProcessing.initParameter(normalEquationInfo);
      changedNormalEquationInfo = FALSE;

      // set variable receiver
      auto fileNameVariableList2 = gnssProcessing.fileNameVariableList;
      if(!variableReceiver.empty())
      {
        const UInt idRecv = std::distance(normalEquationInfo.estimateReceiver.begin(), std::find(normalEquationInfo.estimateReceiver.begin(), normalEquationInfo.estimateReceiver.end(), TRUE));
        addVariable(variableReceiver, gnssProcessing.gnss.receiver.at(idRecv)->name(), fileNameVariableList2);
      }

      // parameterNames
      std::vector<ParameterName> parameterNamesOld = normalEquationInfo.parameterNames();
      if(eliminateEpochParameters)
        parameterNamesOld.erase(parameterNamesOld.begin(), parameterNamesOld.begin()+normalEquationInfo.blockIndex(normalEquationInfo.blockInterval()));

      // compute index vectors and block structure for remaining parameters
      Bool mustReorder = FALSE;
      std::vector<UInt> indexVector(parameterNamesOld.size());
      std::iota(indexVector.begin(), indexVector.end(), 0);
      if(parameterSelector)
      {
        indexVector = parameterSelector->indexVector(parameterNamesOld);
        mustReorder = TRUE;
      }
      const UInt parameterCountOld  = parameterNamesOld.size();
      const UInt parameterCountNew  = indexVector.size();

      UInt defaultBlockSize = this->defaultBlockSize;
      if(defaultBlockSize != NULLINDEX)
        mustReorder = TRUE;
      if(mustReorder && (defaultBlockSize == NULLINDEX))
        defaultBlockSize = 2048;

      // create list of remaining parameter names
      std::vector<ParameterName> parameterNames(indexVector.size());
      for(UInt i=0; i<indexVector.size(); i++)
        parameterNames.at(i) = (indexVector.at(i) != NULLINDEX) ? parameterNamesOld.at(indexVector.at(i)) : ParameterName();

      if(!fileNameParameterNames.empty() && Parallel::isMaster(normalEquationInfo.comm))
      {
        logStatus<<"write parameter name file <"<<fileNameParameterNames(fileNameVariableList2)<<">"<<Log::endl;
        writeFileParameterName(fileNameParameterNames(fileNameVariableList2), parameterNames);
      }

      if(!fileNameApriori.empty())
      {
        logStatus<<"write apriori solution to <"<<fileNameApriori(fileNameVariableList2)<<">"<<Log::endl;
        Vector x0 = gnssProcessing.gnss.aprioriParameter(normalEquationInfo);
        if(eliminateEpochParameters)
          x0 = x0.row(normalEquationInfo.blockIndex(normalEquationInfo.blockInterval()), x0.rows()-normalEquationInfo.blockIndex(normalEquationInfo.blockInterval()));
        if(Parallel::isMaster(normalEquationInfo.comm))
          writeFileMatrix(fileNameApriori(fileNameVariableList2), reorder(x0, indexVector));
      }

      // compute index vectors and block structure for to-be-eliminated parameters
      std::vector<UInt> eliminationIndexVector = ParameterSelector::indexVectorComplement(indexVector, parameterCountOld);
      const UInt eliminationCount = eliminationIndexVector.size();
      std::vector<UInt> blockIndex            = MatrixDistributed::computeBlockIndex(parameterCountNew, defaultBlockSize);
      std::vector<UInt> eliminationBlockIndex = MatrixDistributed::computeBlockIndex(eliminationCount,  defaultBlockSize);
      const UInt eliminationBlocks = eliminationBlockIndex.size()-1;

      if(eliminationCount > 0) // prepend to-be-eliminated parameters to (remaining) index vector and block structure
      {
        for(auto &index : blockIndex)
          index += eliminationCount;
        eliminationBlockIndex.pop_back();
        blockIndex.insert(blockIndex.begin(),   eliminationBlockIndex.begin(),  eliminationBlockIndex.end());
        indexVector.insert(indexVector.begin(), eliminationIndexVector.begin(), eliminationIndexVector.end());
      }

      // setup constraint normals
      MatrixDistributed  normalsConstraints;
      Matrix             rhsConstraints;
      NormalEquationInfo infoConstraints;
      if(!fileNameConstraints.empty() || eliminationCount)
      {
        auto normalEquationInfo2 = normalEquationInfo;
        if(eliminateEpochParameters)
          normalEquationInfo2.estimationType &= ~(Gnss::NormalEquationInfo::MASK_CONSTRAINT & Gnss::NormalEquationInfo::MASK_EPOCH);
        gnssProcessing.buildNormals(normalEquationInfo2, TRUE/*constraintsOnly*/, eliminateEpochParameters/*solveEpochParameters*/, eliminateEpochParameters/*removeEpochParameters*/);
        infoConstraints = NormalEquationInfo(parameterNames, gnssProcessing.lPl, gnssProcessing.obsCount);
        // merge right hand side to one matrix
        rhsConstraints = Matrix(gnssProcessing.normals.parameterCount(), gnssProcessing.lPl.rows());
        if(Parallel::isMaster(normalEquationInfo2.comm))
          for(UInt i=0; i<gnssProcessing.normals.blockCount(); i++)
            copy(gnssProcessing.n.at(i), rhsConstraints.row(gnssProcessing.normals.blockIndex(i), gnssProcessing.normals.blockSize(i)));
        if(mustReorder)
        {
          gnssProcessing.normals.reorder(indexVector, blockIndex);
          rhsConstraints = reorder(rhsConstraints, indexVector);
        }
        normalsConstraints = gnssProcessing.normals;
      }

      // setup unconstrained normals
      if(!fileNameNormals.empty())
      {
        auto normalEquationInfo2 = normalEquationInfo;
        normalEquationInfo2.estimationType &= ~(Gnss::NormalEquationInfo::MASK_CONSTRAINT & ((eliminateEpochParameters) ? ~Gnss::NormalEquationInfo::MASK_EPOCH : Gnss::NormalEquationInfo::MASK_ALL));
        gnssProcessing.buildNormals(normalEquationInfo2, FALSE/*constraintsOnly*/, eliminateEpochParameters/*solveEpochParameters*/, eliminateEpochParameters/*removeEpochParameters*/);
        // merge right hand side to one matrix
        Matrix rhs(gnssProcessing.normals.parameterCount(), gnssProcessing.lPl.rows());
        if(Parallel::isMaster(normalEquationInfo2.comm))
          for(UInt i=0; i<gnssProcessing.normals.blockCount(); i++)
            copy(gnssProcessing.n.at(i), rhs.row(gnssProcessing.normals.blockIndex(i), gnssProcessing.normals.blockSize(i)));

        if(mustReorder)
        {
          gnssProcessing.normals.reorder(indexVector, blockIndex);
          rhs = reorder(rhs, indexVector);
        }

        if(eliminationCount > 0)
        {
          // add constraints
          Bool added = FALSE;
          for(UInt i=0; i<eliminationBlocks; i++)
            for(UInt k=i; k<eliminationBlocks; k++)
              if(normalsConstraints.isBlockUsed(i,k))
              {
                gnssProcessing.normals.setBlock(i, k);
                gnssProcessing.normals.N(i,k) += normalsConstraints.N(i,k);
                added = TRUE;
              }
          if(added)
          {
            rhs.row(0, eliminationCount) += rhsConstraints.row(0, eliminationCount);
            gnssProcessing.lPl      += infoConstraints.lPl;
            gnssProcessing.obsCount += infoConstraints.observationCount;

            // some checks
            for(UInt i=0; i<eliminationBlocks; i++)
              for(UInt k=eliminationBlocks; k<normalsConstraints.blockCount(); k++)
                if(normalsConstraints.isMyRank(i,k) && !isStrictlyZero(normalsConstraints.N(i,k)))
                  throw(Exception("parameter elimination not possible: constraints between eliminated and remaining parameters"));
            for(UInt i=eliminationBlocks; i<normalsConstraints.blockCount(); i++)
              if(normalsConstraints.isMyRank(i,i) && !isStrictlyZero(normalsConstraints.N(i,i)))
              {
                logWarning<<"constraints on remaining parameters: lPl and obsCount of unconstrained normals not fully consistent"<<Log::endl;
                break;
              }
          }

          gnssProcessing.regularizeNotUsedParameters(normalEquationInfo, 0, eliminationBlocks);
          gnssProcessing.normals.cholesky(TRUE, 0, eliminationBlocks, TRUE);
          gnssProcessing.normals.triangularTransSolve(rhs, 0, eliminationBlocks);
          gnssProcessing.normals.eraseBlocks(0, eliminationBlocks);
          gnssProcessing.obsCount -= eliminationCount;
          for(UInt i=0; i<gnssProcessing.lPl.rows(); i++)
            gnssProcessing.lPl(i) -= quadsum(rhs.slice(0, i, eliminationCount, 1)); // lPl = lPl - n1' N1^(-1) n1
          rhs = rhs.row(eliminationCount, rhs.rows()-eliminationCount);
        }

        logStatus<<"write unconstrained normal equations to <"<<fileNameNormals(fileNameVariableList2)<<">"<<Log::endl;
        writeFileNormalEquation(fileNameNormals(fileNameVariableList2), NormalEquationInfo(parameterNames, gnssProcessing.lPl, gnssProcessing.obsCount), gnssProcessing.normals, rhs);
      }

      if(!fileNameConstraints.empty())
      {
        logStatus<<"write constraint normal equations to <"<<fileNameConstraints(fileNameVariableList2)<<">"<<Log::endl;
        if(eliminationCount > 0)
        {
          normalsConstraints.eraseBlocks(0, eliminationBlocks);
          rhsConstraints = rhsConstraints.row(eliminationCount, rhsConstraints.rows()-eliminationCount);
        }
        writeFileNormalEquation(fileNameConstraints(fileNameVariableList2), infoConstraints, normalsConstraints, rhsConstraints);
      }
    }
    catch(std::exception &e)
    {
      GROOPS_RETHROW(e)
    }
  }
};

/***********************************************/

static const char *docstringGnssProcessingStepWriteResults = R"(
\subsection{WriteResults}\label{gnssProcessingStepType:writeResults}
Write estimated solution to output files.

All output files defined in \program{GnssProcessing} are written in this step.
It is usually the last processing step, but can also be used at other points in the
processing in combination with \config{suffix} to write intermediate results, for example
before \configClass{resolveAmbiguities}{gnssProcessingStepType:resolveAmbiguities} to
output the float solution.
)";

/** @brief Write estimated solution to output files.
* @see GnssProcessingStep */
class GnssProcessingStepWriteResults : public GnssProcessingStep
{
  std::string suffix;

public:
  GnssProcessingStepWriteResults(Config &config)
  {
    readConfig(config, "suffix", suffix, Config::OPTIONAL, "", "appended to every output file name (e.g. orbit.G01.<suffix>.dat)");
  }

  void process(GnssProcessing &gnssProcessing, const Vector &/*estimableReceiver*/,
               Gnss::NormalEquationInfo &normalEquationInfo, Bool &changedNormalEquationInfo) const override
  {
    try
    {
      logStatus<<"=== write results  =========================================="<<Log::endl;
      if(changedNormalEquationInfo)
        gnssProcessing.initParameter(normalEquationInfo);
      changedNormalEquationInfo = FALSE;
      gnssProcessing.gnss.writeResults(normalEquationInfo, !suffix.empty() ? "."+suffix : "");
    }
    catch(std::exception &e)
    {
      GROOPS_RETHROW(e)
    }
  }
};

/***********************************************/

static const char *docstringGnssProcessingStepPrintResidualStatistics = R"(
\subsection{PrintResidualStatistics}\label{gnssProcessingStepType:printResidualStatistics}
Print residual statistics.
)";

/** @brief Print residual statistics.
* @see GnssProcessingStep */
class GnssProcessingStepPrintResidualStatistics : public GnssProcessingStep
{
public:
  GnssProcessingStepPrintResidualStatistics(Config &/*config*/) {}

  void process(GnssProcessing &gnssProcessing, const Vector &/*estimableReceiver*/,
               Gnss::NormalEquationInfo &normalEquationInfo, Bool &changedNormalEquationInfo) const override
  {
    try
    {
      logStatus<<"=== print residual statistics  =============================="<<Log::endl;
      if(changedNormalEquationInfo)
        gnssProcessing.initParameter(normalEquationInfo);
      changedNormalEquationInfo = FALSE;
      logInfo<<"################ receiver #####################"<<Log::endl;
      for(UInt idRecv=0; idRecv<gnssProcessing.gnss.receiver.size(); idRecv++)
        if(normalEquationInfo.estimateReceiver.at(idRecv))
        {
          std::vector<GnssType> type = gnssProcessing.gnss.types(~GnssType::FREQ_NO);
          std::vector<Double>   ePe(type.size(), 0), redundancy(type.size(), 0);
          std::vector<UInt>     obsCount(type.size(), 0), outlierCount(type.size(), 0);
          gnssProcessing.gnss.receiver.at(idRecv)->residualsStatistics(NULLINDEX, normalEquationInfo.idEpochs, type, ePe, redundancy, obsCount, outlierCount);
          for(UInt i=0; i<type.size(); i++)
          {
            Parallel::reduceSum(obsCount.at(i), 0, normalEquationInfo.comm);
            Parallel::broadCast(obsCount.at(i), 0, normalEquationInfo.comm);
            if(obsCount.at(i))
            {
              Double factor = 0;
              if(gnssProcessing.gnss.receiver.at(idRecv)->useable())
                factor = gnssProcessing.gnss.receiver.at(idRecv)->sigmaFactor(type.at(i));
              Parallel::reduceSum(factor, 0, normalEquationInfo.comm);
              Parallel::reduceSum(ePe.at(i),          0, normalEquationInfo.comm);
              Parallel::reduceSum(redundancy.at(i),   0, normalEquationInfo.comm);
              Parallel::reduceSum(outlierCount.at(i), 0, normalEquationInfo.comm);
              if(Parallel::isMaster(normalEquationInfo.comm))
                logInfo<<"  "<<gnssProcessing.gnss.receiver.at(idRecv)->name()<<": "<<type.at(i).str()
                      <<": factor = "    <<factor%"%5.2f"s
                      <<", sigma0 = "    <<Vce::standardDeviation(ePe.at(i), redundancy.at(i), 2.5/*huber*/, 1.5/*huberPower*/)%"%4.2f"s
                      <<", count = "     <<obsCount.at(i)%"%5i"s
                      <<", outliers = "  <<outlierCount.at(i)%"%5i"s<<" ("<<(100.*outlierCount.at(i)/obsCount.at(i))%"%4.2f"s<<" %)"<<Log::endl;
            }
          }
        } // for(idRecv)

      logInfo<<"################ transmitter #####################"<<Log::endl;
      for(UInt idTrans=0; idTrans<gnssProcessing.gnss.transmitter.size(); idTrans++)
      {
        Matrix A;
        std::vector<GnssType> typesStation = gnssProcessing.gnss.types(GnssType::ALL);
        std::vector<GnssType> types;
        gnssProcessing.gnss.defaultSignalComposition(typesStation, types, A);
        std::vector<Double>   ePe(types.size(), 0), redundancy(types.size(), 0);
        std::vector<UInt>     obsCount(types.size(), 0), outlierCount(types.size(), 0), stationCount(types.size(), 0);
        for(UInt idRecv=0; idRecv<gnssProcessing.gnss.receiver.size(); idRecv++)
          if(normalEquationInfo.estimateReceiver.at(idRecv))
          {
            std::vector<Double>   ePe2(typesStation.size(), 0), redundancy2(typesStation.size(), 0);
            std::vector<UInt>     obsCount2(typesStation.size(), 0), outlierCount2(typesStation.size(), 0);
            gnssProcessing.gnss.receiver.at(idRecv)->residualsStatistics(idTrans, normalEquationInfo.idEpochs, typesStation, ePe2, redundancy2, obsCount2, outlierCount2);
            // transform to transmitter types
            for(UInt i=0; i<types.size(); i++)
            {
              Bool found = FALSE;
              for(UInt k=0; k<typesStation.size(); k++)
                if(obsCount2.at(k) && A(k,i))
                {
                  found = TRUE;
                  obsCount.at(i)     += obsCount2.at(k);
                  ePe.at(i)          += ePe2.at(k);
                  redundancy.at(i)   += redundancy2.at(k);
                  outlierCount.at(i) += outlierCount2.at(k);
                }
              if(found)
                stationCount.at(i)++;
            }
          }
        for(UInt i=0; i<types.size(); i++)
        {
          Parallel::reduceSum(obsCount.at(i), 0, normalEquationInfo.comm);
          Parallel::broadCast(obsCount.at(i), 0, normalEquationInfo.comm);
          if(obsCount.at(i)!=0)
          {
            Parallel::reduceSum(ePe.at(i),          0, normalEquationInfo.comm);
            Parallel::reduceSum(redundancy.at(i),   0, normalEquationInfo.comm);
            Parallel::reduceSum(outlierCount.at(i), 0, normalEquationInfo.comm);
            Parallel::reduceSum(stationCount.at(i), 0, normalEquationInfo.comm);
            if(Parallel::isMaster(normalEquationInfo.comm))
              logInfo<<"  "<<gnssProcessing.gnss.transmitter.at(idTrans)->PRN().str()<<": "<<types.at(i).str()
                    <<": sigma0 = "           <<Vce::standardDeviation(ePe.at(i), redundancy.at(i), 2.5/*huber*/, 1.5/*huberPower*/)%"%4.2f"s
                    <<", observing stations = "<<stationCount.at(i)%"%4i"s
                    <<", count = "             <<obsCount.at(i)%"%6i"s
                    <<", outliers = "          <<outlierCount.at(i)%"%4i"s<<" ("<<(100.*outlierCount.at(i)/obsCount.at(i))%"%4.2f"s<<" %)"<<Log::endl;
          }
        }
      } // for(idTrans)
    }
    catch(std::exception &e)
    {
      GROOPS_RETHROW(e)
    }
  }
};

/***********************************************/

static const char *docstringGnssProcessingStepSelectParameters = R"(
\subsection{SelectParameters}\label{gnssProcessingStepType:selectParameters}
Select parameters for subsequent steps.

This step can be used to control which parameters are estimated and which constraints are applied in subsequent
\configClass{estimate}{gnssProcessingStepType:estimate} steps. An example would be to estimate only the TEC bias
in an intermediate step as seen in \reference{GNSS satellite orbit determination and network analysis}{cookbook.gnssNetwork:processing}.
Another example is to process at a 5-minute sampling using \configClass{selectEpochs}{gnssProcessingStepType:selectEpochs} and then at
the end to densify the clock parameters to the full 30-second observation sampling while keeping all other parameters fixed.
)";

/** @brief Select parameters for subsequent steps.
* @see GnssProcessingStep */
class GnssProcessingStepSelectParameters : public GnssProcessingStep
{
  Gnss::NormalEquationInfo::EstimationType estimationType;

public:
  GnssProcessingStepSelectParameters(Config &config)
  {
    estimationType = 0;
    for(const auto &estimationType : Gnss::NormalEquationInfo::estimationTypes)
    {
      Bool flag;
      readConfig(config, std::get<1>(estimationType), flag, Config::DEFAULT, std::get<2>(estimationType), std::get<3>(estimationType));
      if(flag) this->estimationType |=  std::get<0>(estimationType);
    }
  }

  void process(GnssProcessing &/*gnssProcessing*/, const Vector &estimableReceiver,
               Gnss::NormalEquationInfo &normalEquationInfo, Bool &changedNormalEquationInfo) const override
  {
    try
    {
      logStatus<<"=== select parameters ======================================="<<Log::endl;
      if(normalEquationInfo.estimationType != estimationType)
        changedNormalEquationInfo = TRUE;
      normalEquationInfo.estimationType = estimationType;
      // for single receiver: allow receiver parameters only
      if(sum(estimableReceiver) == 1)
        normalEquationInfo.estimationType &= Gnss::NormalEquationInfo::MASK_SINGLERECEIVER;
      for(const auto &estimationType : Gnss::NormalEquationInfo::estimationTypes)
        if(normalEquationInfo.estimationType & std::get<0>(estimationType))
          logInfo<<"  "<<std::get<1>(estimationType)<<Log::endl;
    }
    catch(std::exception &e)
    {
      GROOPS_RETHROW(e)
    }
  }
};

/***********************************************/

static const char *docstringGnssProcessingStepSelectEpochs = R"(
\subsection{SelectEpochs}\label{gnssProcessingStepType:selectEpochs}
Select epochs for subsequent steps.

This step can be used to reduce the processing sampling while keeping the original observation sampling
for all preprocessing steps (e.g. outlier and cycle slip detection).
Another example is to process at a 5-minute sampling by setting \config{nthEpoch}=\verb|10| and then
at the end to densify only the clock parameters (see \configClass{selectParameters}{gnssProcessingStepType:selectParameters})
to the full 30-second observation sampling by setting \config{nthEpoch}=\verb|1| while keeping all other parameters fixed.
)";

/** @brief Select epochs for subsequent steps.
* @see GnssProcessingStep */
class GnssProcessingStepSelectEpochs : public GnssProcessingStep
{
  UInt nthEpoch;

public:
  GnssProcessingStepSelectEpochs(Config &config)
  {
    readConfig(config, "nthEpoch", nthEpoch, Config::MUSTSET, "1", "use only every nth epoch in all subsequent processing steps");
  }

  void process(GnssProcessing &gnssProcessing, const Vector &/*estimableReceiver*/,
               Gnss::NormalEquationInfo &normalEquationInfo, Bool &changedNormalEquationInfo) const override
  {
    try
    {
      logStatus<<"=== select epochs ==========================================="<<Log::endl;
      normalEquationInfo.idEpochs.clear();
      for(UInt idEpoch=0; idEpoch<gnssProcessing.gnss.times.size(); idEpoch+=nthEpoch)
        normalEquationInfo.idEpochs.push_back(idEpoch);
      logInfo<<"  "<<normalEquationInfo.idEpochs.size()<<" of "<<gnssProcessing.gnss.times.size()<<" epochs selected"<<Log::endl;
      changedNormalEquationInfo = TRUE;
    }
    catch(std::exception &e)
    {
      GROOPS_RETHROW(e)
    }
  }
};

/***********************************************/

static const char *docstringGnssProcessingStepSelectNormalsBlockStructure = R"(
\subsection{SelectNormalsBlockStructure}\label{gnssProcessingStepType:selectNormalsBlockStructure}
Select block structure of sparse normal equations for subsequent steps.

This step can be used to define the structure of the different parts of the normal equation system,
which can have a major impact on computing performance and memory consumption depending on the processing setup.

\fig{!hb}{0.4}{gnss_normals_structure}{fig:gnss_normals_structure}{Structure of normal equations in GNSS processing}

The normal equation system is divided into three parts for epoch, interval, and ambiguity parameters.
The epoch part is subdivided further into one subpart per epoch. Each part is divided into blocks and only non-zero
blocks are stored in memory to reduce memory consumption and to prevent unnecessary matrix computations.
\config{defaultBlockSizeEpoch}, \config{defaultBlockSizeInterval}, and \config{defaultBlockSizeAmbiguity} control
the size of the blocks within each part of the normal equations. \config{defaultBlockReceiverCount} can be set to group
a number of receivers into one block within the epoch and interval parts.

If \config{keepEpochNormalsInMemory}=\verb|no| epoch blocks are eliminated after they are set up to reduce the number
of parameters in the normal equation system. \config{defaultBlockCountReduction} controls after how many epoch blocks
an elimination step is performed. For larger processing setups or high sampling rates epoch block elimination is recommended
as the large number of clock parameters require a lot of memory.
)";

/** @brief Select block structure of sparse normal equations for subsequent steps.
* @see GnssProcessingStep */
class GnssProcessingStepSelectNormalsBlockStructure : public GnssProcessingStep
{
  UInt defaultBlockSizeEpoch, defaultBlockSizeInterval, defaultBlockSizeAmbiguity;
  UInt defaultBlockReceiverCount;
  UInt defaultBlockCountReduction;
  Bool keepEpochNormalsInMemory;
  Bool accumulateEpochObservations;

public:
  GnssProcessingStepSelectNormalsBlockStructure(Config &config)
  {
    readConfig(config, "defaultBlockSizeEpoch",       defaultBlockSizeEpoch,       Config::DEFAULT,  "0",  "block size of epoch parameters, 0: one block");
    readConfig(config, "defaultBlockSizeInterval",    defaultBlockSizeInterval,    Config::DEFAULT,  "64", "block size of interval parameters, 0: one block");
    readConfig(config, "defaultBlockSizeAmbiguity",   defaultBlockSizeAmbiguity,   Config::DEFAULT,  "64", "block size of ambiguity parameters, 0: one block");
    readConfig(config, "defaultBlockReceiverCount",   defaultBlockReceiverCount,   Config::DEFAULT,  "0",  "number of receivers to group into one block for epoch and interval");
    readConfig(config, "defaultBlockCountReduction",  defaultBlockCountReduction,  Config::DEFAULT,  "32", "minimum number of blocks for epoch reduction");
    readConfig(config, "keepEpochNormalsInMemory",    keepEpochNormalsInMemory,    Config::DEFAULT,  "1",  "speeds up processing but uses much more memory");
    readConfig(config, "accumulateEpochObservations", accumulateEpochObservations, Config::DEFAULT,  "0",  "set up all observations per epoch and receiver at once");
  }

  void process(GnssProcessing &/*gnssProcessing*/, const Vector &/*estimableReceiver*/,
               Gnss::NormalEquationInfo &normalEquationInfo, Bool &changedNormalEquationInfo) const override
  {
    try
    {
      logStatus<<"=== select block structure =================================="<<Log::endl;
      normalEquationInfo.defaultBlockSizeEpoch       = defaultBlockSizeEpoch;
      normalEquationInfo.defaultBlockSizeInterval    = defaultBlockSizeInterval;
      normalEquationInfo.defaultBlockSizeAmbiguity   = defaultBlockSizeAmbiguity;
      normalEquationInfo.defaultBlockReceiverCount   = defaultBlockReceiverCount;
      normalEquationInfo.defaultBlockCountReduction  = defaultBlockCountReduction;
      normalEquationInfo.keepEpochNormalsInMemory    = keepEpochNormalsInMemory;
      normalEquationInfo.accumulateEpochObservations = accumulateEpochObservations;
      logInfo<<"  blockSizeEpoch              = "<<normalEquationInfo.defaultBlockSizeEpoch    <<Log::endl;
      logInfo<<"  blockSizeInterval           = "<<normalEquationInfo.defaultBlockSizeInterval <<Log::endl;
      logInfo<<"  blockSizeAmbiguity          = "<<normalEquationInfo.defaultBlockSizeAmbiguity <<Log::endl;
      logInfo<<"  blockReceiverCount          = "<<normalEquationInfo.defaultBlockReceiverCount<<Log::endl;
      logInfo<<"  blockCountReduction         = "<<normalEquationInfo.defaultBlockCountReduction<<Log::endl;
      logInfo<<"  keepEpochNormalsInMemory    = "<<(normalEquationInfo.keepEpochNormalsInMemory ? "yes" : "no")<<Log::endl;
      logInfo<<"  accumulateEpochObservations = "<<(normalEquationInfo.accumulateEpochObservations ? "yes" : "no")<<Log::endl;
      changedNormalEquationInfo = TRUE;
    }
    catch(std::exception &e)
    {
      GROOPS_RETHROW(e)
    }
  }
};

/***********************************************/

static const char *docstringGnssProcessingStepSelectReceivers = R"(
\subsection{SelectReceivers}\label{gnssProcessingStepType:selectReceivers}
Select receivers of a \configClass{stationNetwork}{gnssParametrizationReceiverType:stationNetwork} to process in subsequent steps.

This step can be used to process only a subset of stations in subsequent processing steps.
The most common use is to start the processing with a well-distributed network of core stations as seen in
\reference{GNSS satellite orbit determination and network analysis}{cookbook.gnssNetwork:processing}.
In this case \configFile{inputfileStationList}{stringList} should contain one line per station with alternative
stations on the same line separated by whitespace. To later process all other stations individually, use
the processing step \configClass{forEachReceiverSeparately}{gnssProcessingStepType:forEachReceiverSeparately}
and provide the same station list as \configFile{inputfileExcludeStationList}{stringList} in that step.
)";

/** @brief Select receivers to process in subsequent steps.
* @see GnssProcessingStep */
class GnssProcessingStepSelectReceivers : public GnssProcessingStep
{
  std::vector<std::vector<std::string>> receiverName, receiverExcludeName;

public:
  GnssProcessingStepSelectReceivers(Config &config)
  {
    FileName fileNameStationList, fileNameExcludeList;
    readConfig(config, "inputfileStationList",        fileNameStationList, Config::OPTIONAL, "", "empty = all stations");
    readConfig(config, "inputfileExcludeStationList", fileNameExcludeList, Config::OPTIONAL, "", "exclude these stations");

    if(!fileNameStationList.empty()) readFileStringTable(fileNameStationList, receiverName);
    if(!fileNameExcludeList.empty()) readFileStringTable(fileNameExcludeList, receiverExcludeName);
  }

  void process(GnssProcessing &gnssProcessing, const Vector &estimableReceiver,
               Gnss::NormalEquationInfo &normalEquationInfo, Bool &changedNormalEquationInfo) const override
  {
    try
    {
      logStatus<<"=== select receivers ========================================"<<Log::endl;
      normalEquationInfo.estimateReceiver = selectReceivers(gnssProcessing.gnss, estimableReceiver, receiverName, receiverExcludeName);
      logInfo<<"  "<<std::count(normalEquationInfo.estimateReceiver.begin(), normalEquationInfo.estimateReceiver.end(), TRUE)<<" receivers selected"<<Log::endl;
      changedNormalEquationInfo = TRUE;
    }
    catch(std::exception &e)
    {
      GROOPS_RETHROW(e)
    }
  }
};

/***********************************************/

static const char *docstringGnssProcessingStepForEachReceiverSeparately = R"(
\subsection{ForEachReceiverSeparately}\label{gnssProcessingStepType:forEachReceiverSeparately}
Perform these processing steps for each receiver of a \configClass{stationNetwork}{gnssParametrizationReceiverType:stationNetwork} separately.
All \configClass{transmitter}{gnssParametrizationTransmitterType} parameters are disabled in these processing steps.

This step can be used for individual precise point positioning (PPP) of all stations.
During \reference{GNSS satellite orbit determination and network analysis}{cookbook.gnssNetwork:processing} this step is used after the
initial processing of the core network to process all other stations individually. In that case provide the same station list as
\configFile{inputfileExcludeStationList}{stringList} in this step that was used as \configFile{inputfileStationList}{stringList} in the
\configClass{selectReceivers}{gnssProcessingStepType:selectReceivers} step where the core network was selected.
)";

/** @brief Perform these processing steps for each receiver separately.
* @see GnssProcessingStep */
class GnssProcessingStepForEachReceiverSeparately : public GnssProcessingStep
{
  std::vector<std::vector<std::string>> receiverName, receiverExcludeName;
  std::vector<GnssProcessingStepPtr>    processingSteps;

public:
  GnssProcessingStepForEachReceiverSeparately(Config &config)
  {
    FileName fileNameStationList, fileNameExcludeList;
    readConfig(config, "inputfileStationList",        fileNameStationList, Config::OPTIONAL, "", "empty = all stations");
    readConfig(config, "inputfileExcludeStationList", fileNameExcludeList, Config::OPTIONAL, "", "exclude these stations");
    readConfig(config, "processingStep",              processingSteps,     Config::MUSTSET,  "", "steps are processed consecutively");

    if(!fileNameStationList.empty()) readFileStringTable(fileNameStationList, receiverName);
    if(!fileNameExcludeList.empty()) readFileStringTable(fileNameExcludeList, receiverExcludeName);
  }

  void process(GnssProcessing &gnssProcessing, const Vector &estimableReceiver,
               Gnss::NormalEquationInfo &normalEquationInfo, Bool &changedNormalEquationInfo) const override
  {
    try
    {
      logStatus<<"=== for each receiver separatley ============================"<<Log::endl;
      auto estimateSingleReceiver = selectReceivers(gnssProcessing.gnss, estimableReceiver, receiverName, receiverExcludeName);
      logInfo<<"  "<<std::count(estimateSingleReceiver.begin(), estimateSingleReceiver.end(), TRUE)<<" receivers selected"<<Log::endl;
      if(Parallel::size(normalEquationInfo.comm) > 1)
        logInfo<<"  Only results of a subset of stations are displayed in the following"<<Log::endl;

      for(UInt idRecv=0; idRecv<gnssProcessing.gnss.receiver.size(); idRecv++)
        if(estimateSingleReceiver.at(idRecv) && gnssProcessing.gnss.receiver.at(idRecv)->useable())
        {
          logStatus<<"=== select single receiver ("<<gnssProcessing.gnss.receiver.at(idRecv)->name()<<") ==========================="<<Log::endl;
          Vector estimableReceiverSingle(estimableReceiver.size());
          estimableReceiverSingle(idRecv) = TRUE;

          Gnss::NormalEquationInfo normalEquationInfo2 = normalEquationInfo;
          std::fill(normalEquationInfo2.estimateReceiver.begin(), normalEquationInfo2.estimateReceiver.end(), FALSE);
          normalEquationInfo2.estimationType &= Gnss::NormalEquationInfo::MASK_SINGLERECEIVER;
          normalEquationInfo2.estimateReceiver.at(idRecv) = TRUE;
          normalEquationInfo2.comm = Parallel::selfCommunicator();
          changedNormalEquationInfo = TRUE;

          try
          {
            for(const auto &processingStep : processingSteps)
              processingStep->process(gnssProcessing, estimableReceiverSingle, normalEquationInfo2, changedNormalEquationInfo);
          }
          catch(std::exception &e)
          {
            logWarning<<gnssProcessing.gnss.receiver.at(idRecv)->name()<<": disabled due to exception in single receiver loop:\n"<<e.what()<<Log::endl;
            gnssProcessing.gnss.receiver.at(idRecv)->disable();
          }
        } // for(idRecv)

      changedNormalEquationInfo = TRUE;
    }
    catch(std::exception &e)
    {
      GROOPS_RETHROW(e)
    }
  }
};

/***********************************************/

static const char *docstringGnssProcessingStepGroup = R"(
\subsection{Group}\label{gnssProcessingStepType:group}
Perform these processing steps. This step can be used to structure complex processing flows and has no further effect on the processing itself.
)";

/** @brief Perform these processing steps. This step can be used to structure complex processing flows.
* @see GnssProcessingStep */
class GnssProcessingStepGroup : public GnssProcessingStep
{
  std::vector<GnssProcessingStepPtr> processingSteps;

public:
  GnssProcessingStepGroup(Config &config)
  {
    readConfig(config, "processingStep", processingSteps, Config::MUSTSET,   "",     "steps are processed consecutively");
  }

  void process(GnssProcessing &gnssProcessing, const Vector &estimableReceiver,
               Gnss::NormalEquationInfo &normalEquationInfo, Bool &changedNormalEquationInfo) const override
  {
    try
    {
      logStatus<<"=== process group ==========================================="<<Log::endl;
      Gnss::NormalEquationInfo normalEquationInfo2 = normalEquationInfo;
      for(auto &processingStep : processingSteps)
        processingStep->process(gnssProcessing, estimableReceiver, normalEquationInfo2, changedNormalEquationInfo);
      changedNormalEquationInfo =  (normalEquationInfo2.estimationType   != normalEquationInfo.estimationType)
                                || (normalEquationInfo2.estimateReceiver != normalEquationInfo.estimateReceiver)
                                || (normalEquationInfo2.idEpochs         != normalEquationInfo.idEpochs);
    }
    catch(std::exception &e)
    {
      GROOPS_RETHROW(e)
    }
  }
};

/// @}

/***********************************************/

GROOPS_REGISTER_CLASS(GnssProcessingStep, "gnssProcessingStepType",
                      GnssProcessingStepEstimate,
                      GnssProcessingStepResolveAmbiguities,
                      GnssProcessingStepComputeCovarianceMatrix,
                      GnssProcessingStepWriteNormalEquations,
                      GnssProcessingStepWriteResults,
                      GnssProcessingStepPrintResidualStatistics,
                      GnssProcessingStepSelectParameters,
                      GnssProcessingStepSelectEpochs,
                      GnssProcessingStepSelectNormalsBlockStructure,
                      GnssProcessingStepSelectReceivers,
                      GnssProcessingStepForEachReceiverSeparately,
                      GnssProcessingStepGroup)

GROOPS_READCONFIG_CLASS(GnssProcessingStep, "gnssProcessingStepType")

/***********************************************/

GnssProcessingStepPtr GnssProcessingStep::create(Config &config, const std::string &name)
{
  try
  {
    GnssProcessingStepPtr ptr;
    std::string           choice;

    readConfigChoice(config, name, choice, Config::MUSTSET, "", "");
    if(readConfigChoiceElement(config, "estimate",                    choice, "least squares adjustment"))
      ptr = std::make_shared<GnssProcessingStepEstimate>(config);
    if(readConfigChoiceElement(config, "resolveAmbiguities",          choice, "resolve integer ambiguities"))
      ptr = std::make_shared<GnssProcessingStepResolveAmbiguities>(config);
    if(readConfigChoiceElement(config, "computeCovarianceMatrix",     choice, "compute covariance matrix"))
      ptr = std::make_shared<GnssProcessingStepComputeCovarianceMatrix>(config);
    if(readConfigChoiceElement(config, "writeNormalEquations",        choice, "write unconstrained and constraint normal equations"))
      ptr = std::make_shared<GnssProcessingStepWriteNormalEquations>(config);
    if(readConfigChoiceElement(config, "writeResults",                choice, "write all estimated parameters"))
      ptr = std::make_shared<GnssProcessingStepWriteResults>(config);
    if(readConfigChoiceElement(config, "printResidualStatistics",     choice, "print residual statistics"))
      ptr = std::make_shared<GnssProcessingStepPrintResidualStatistics>(config);
    if(readConfigChoiceElement(config, "selectParameters",            choice, "select parameters for all subsequent processing steps"))
      ptr = std::make_shared<GnssProcessingStepSelectParameters>(config);
    if(readConfigChoiceElement(config, "selectEpochs",                choice, "select epochs to be used in all subsequent processing steps"))
      ptr = std::make_shared<GnssProcessingStepSelectEpochs>(config);
    if(readConfigChoiceElement(config, "selectNormalsBlockStructure", choice, "select block structure of the distributed normal equation matrix for all subsequent processing steps"))
      ptr = std::make_shared<GnssProcessingStepSelectNormalsBlockStructure>(config);
    if(readConfigChoiceElement(config, "selectReceivers",             choice, "use this subset of stations in all subsequent processing steps"))
      ptr = std::make_shared<GnssProcessingStepSelectReceivers>(config);
    if(readConfigChoiceElement(config, "forEachReceiverSeparately",   choice, "process receiver by receiver, transmitter parameters disabled"))
      ptr = std::make_shared<GnssProcessingStepForEachReceiverSeparately>(config);
    if(readConfigChoiceElement(config, "group",                       choice, "group processing steps"))
      ptr = std::make_shared<GnssProcessingStepGroup>(config);
    endChoice(config);

    return ptr;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

std::vector<Byte> GnssProcessingStep::selectReceivers(const Gnss &gnss, const Vector &estimableReceiver,
                                                      const std::vector<std::vector<std::string>> &selectNames,
                                                      const std::vector<std::vector<std::string>> &excludeNames) const
{
  try
  {
    // get all receiver names
    std::vector<std::string> receiverName(gnss.receiver.size());
    for(UInt idRecv=0; idRecv<gnss.receiver.size(); idRecv++)
      receiverName.at(idRecv) = gnss.receiver.at(idRecv)->name();

    std::vector<Byte> estimateReceivers(estimableReceiver.rows());
    for(UInt idRecv=0; idRecv<estimateReceivers.size(); idRecv++)
      estimateReceivers.at(idRecv) = (static_cast<Byte>(estimableReceiver(idRecv)) != 0);

    // station list
    if(selectNames.size())
    {
      std::fill(estimateReceivers.begin(), estimateReceivers.end(), FALSE);
      for(UInt k=0; k<selectNames.size(); k++)
        for(UInt l=0; l<selectNames.at(k).size(); l++) // alternatives
        {
          const UInt idRecv = std::distance(receiverName.begin(), std::find(receiverName.begin(), receiverName.end(), selectNames.at(k).at(l)));
          if((idRecv >= receiverName.size()) || !estimableReceiver(idRecv))
            continue;
          estimateReceivers.at(idRecv) = TRUE;
          break; // skip alternative stations
        }
    }

    // exclude station list
    for(UInt k=0; k<excludeNames.size(); k++)
      for(UInt l=0; l<excludeNames.at(k).size(); l++) // alternatives
      {
        const UInt idRecv = std::distance(receiverName.begin(), std::find(receiverName.begin(), receiverName.end(), excludeNames.at(k).at(l)));
        if((idRecv >= receiverName.size()) || !estimableReceiver(idRecv))
          continue;
        estimateReceivers.at(idRecv) = FALSE;
        break; // skip alternative stations
      }

    return estimateReceivers;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/

void GnssProcessing::run(Config &config)
{
  try
  {
    TimeSeriesPtr                                  timeSeriesPtr, timesIntervalPtr;
    std::vector<GnssParametrizationTransmitterPtr> gnssTransmitters;
    std::vector<GnssParametrizationReceiverPtr>    gnssReceivers;
    GnssParametrizationEarthRotationPtr            gnssEarthRotation;
    GnssParametrizationGravityFieldPtr             gnssGravityField;
    GnssParametrizationAmbiguitiesPtr              gnssAmbiguities;
    std::vector<GnssParametrizationConstraintsPtr> gnssConstraints;

    readConfig(config, "outputfileParameterName",  fileNameNormalsInfo, Config::OPTIONAL, "output/parameterNames_{loopTime:%D}.txt", "parameter names");
    readConfig(config, "intervals",                timesIntervalPtr,    Config::OPTIONAL, "",     "each interval is processed independently");
    readConfig(config, "parallelIntervals",        parallelIntervals,   Config::DEFAULT,  "1",    "parallelize intervals instead of receivers per interval");
    readConfig(config, "timeSeries",               timeSeriesPtr,       Config::MUSTSET,  "",     "defines epochs within intervals");
    readConfig(config, "timeMargin",               marginSeconds,       Config::DEFAULT,  "0.1",  "[seconds] margin to consider two times identical");
    readConfig(config, "transmitter",              gnssTransmitters,    Config::MUSTSET,  "",     "constellation of GNSS satellites");
    readConfig(config, "receiver",                 gnssReceivers,       Config::MUSTSET,  "",     "ground station network or LEO satellite");
    readConfig(config, "ionosphere",               gnssIonosphere,      Config::MUSTSET,  "",     "ionosphere settings");
    readConfig(config, "earthRotation",            gnssEarthRotation,   Config::MUSTSET,  "",     "Earth rotation and EOP estimation");
    readConfig(config, "gravityField",             gnssGravityField,    Config::OPTIONAL, "",     "gravity field estimation");
    readConfig(config, "ambiguities",              gnssAmbiguities,     Config::MUSTSET,  "",     "ambiguity resolution settings");
    readConfig(config, "constraints",              gnssConstraints,     Config::OPTIONAL, "",     "additional constraints");
    readConfig(config, "processingStep",           processingSteps,     Config::MUSTSET,  "",     "steps are processed consecutively");
    if(isCreateSchema(config)) return;

    // ============================

    // init time intervals
    // -------------------
    timeSeries = timeSeriesPtr->times();
    if(timesIntervalPtr)
      timesInterval = timesIntervalPtr->times();
    if(timesInterval.size()==0)
    {
      timesInterval.push_back(timeSeries.at(0));
      timesInterval.push_back(timeSeries.back()+seconds2time(1));
    }

    addTimeVariables(fileNameVariableList);

    // ============================

    // init
    // ----
    std::vector<Gnss::TransmitterPtr>     transmitters;
    std::vector<Gnss::ReceiverPtr>        receivers;
    std::vector<Gnss::ParametrizationPtr> parametrizations;

    for(auto &gnssTransmitter : gnssTransmitters)
    {
      auto t = gnssTransmitter->transmitters();
      auto p = gnssTransmitter->parametrizations();
      transmitters.insert(transmitters.end(), t.begin(), t.end());
      parametrizations.insert(parametrizations.end(), p.begin(), p.end());
    }

    for(auto &gnssReceiver : gnssReceivers)
    {
      auto r = gnssReceiver->receivers();
      auto p = gnssReceiver->parametrizations();
      receivers.insert(receivers.end(), r.begin(), r.end());
      parametrizations.insert(parametrizations.end(), p.begin(), p.end());
    }

    // constraints
    parametrizations.insert(parametrizations.end(), gnssConstraints.begin(), gnssConstraints.end());

    logInfo<<"Init GNSS"<<Log::endl;
    gnss.init(receivers, transmitters, parametrizations, gnssIonosphere, gnssEarthRotation, gnssGravityField, gnssAmbiguities);
    logInfo<<"  transmitter: "<<transmitters.size()<<Log::endl;
    logInfo<<"  receiver:    "<<receivers.size()<<Log::endl;

    // ============================

    logStatus<<"compute intervals"<<Log::endl;
    if(parallelIntervals)
    {
      Parallel::forEach(timesInterval.size()-1, [this](UInt idInterval) {computeInterval(idInterval, Parallel::selfCommunicator());});
    }
    else
    {
      for(UInt idInterval=0; idInterval<timesInterval.size()-1; idInterval++)
        computeInterval(idInterval, nullptr); // use all processes for each interval
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssProcessing::computeInterval(UInt idInterval, Parallel::CommunicatorPtr comm)
{
  try
  {
    logStatus<<"================================================================"<<Log::endl;
    logStatus<<idInterval+1<<". "<<timesInterval.at(idInterval).dateTimeStr()<<" - "<<timesInterval.at(idInterval+1).dateTimeStr()<<Log::endl;
    evaluateTimeVariables(idInterval, timesInterval.at(idInterval), timesInterval.at(idInterval+1), fileNameVariableList);

    // init time series
    // ----------------
    std::vector<Time> times;
    for(UInt i=0; i<timeSeries.size(); i++)
      if(timeSeries.at(i).isInInterval(timesInterval.at(idInterval), timesInterval.at(idInterval+1)))
        times.push_back(timeSeries.at(i));

    if(!times.size() && Parallel::isMaster(comm))
    {
      logWarning<<idInterval+1<<". "<<timesInterval.at(idInterval).dateTimeStr()<<" - "<<timesInterval.at(idInterval+1).dateTimeStr()<<": empty interval"<<Log::endl;
      return;
    }

    // init the GNSS system
    // --------------------
    Gnss::AnalysisType analysisType = Gnss::ANALYSIS_RAW;
    gnss.initInterval(analysisType, times, seconds2time(marginSeconds), comm);
    if(!std::any_of(gnss.transmitter.begin(), gnss.transmitter.end(), [](auto trans){return trans->useable();}) && Parallel::isMaster(comm))
    {
      logWarning<<idInterval+1<<". "<<timesInterval.at(idInterval).dateTimeStr()<<" - "<<timesInterval.at(idInterval+1).dateTimeStr()<<": no useable transmitter"<<Log::endl;
      return;
    }
    UInt useable = std::any_of(gnss.receiver.begin(), gnss.receiver.end(), [](auto recv){return recv->useable();});
    Parallel::reduceSum(useable, 0, comm);
    Parallel::broadCast(useable, 0, comm);
    if(!useable && Parallel::isMaster(comm))
    {
      logWarning<<idInterval+1<<". "<<timesInterval.at(idInterval).dateTimeStr()<<" - "<<timesInterval.at(idInterval+1).dateTimeStr()<<": no useable receiver"<<Log::endl;
      return;
    }

    // count observation types
    // -----------------------
    logInfo<<"types and number of observations:"<<Log::endl;
    std::vector<GnssType> types = gnss.types(~GnssType::FREQ_NO);
    Vector countTypes(types.size());
    UInt   countTracks = 0;
    for(UInt idRecv=0; idRecv<gnss.receiver.size(); idRecv++)
      if(gnss.receiver.at(idRecv)->useable())
      {
        std::vector<UInt> count(types.size());
        gnss.receiver.at(idRecv)->countObservationsPerType(NULLINDEX, NULLINDEX, NULLINDEX, types, count);
        for(UInt i=0; i<count.size(); i++)
          countTypes(i) += count.at(i);
        countTracks += gnss.receiver.at(idRecv)->track.size();
      }
    Parallel::reduceSum(countTypes, 0, comm);
    Parallel::reduceSum(countTracks, 0, comm);

    for(UInt idType=0; idType<types.size(); idType++)
      logInfo<<"  "<<types.at(idType).str()<<":"<<countTypes(idType)%"%10i"s<<Log::endl;
    logInfo<<"          ========="<<Log::endl;
    logInfo<<"  total:"<<sum(countTypes)%"%11i"s<<Log::endl;
    logInfo<<"number of tracks: "<<countTracks<<Log::endl;

    Gnss::NormalEquationInfo normalEquationInfo(times.size(), gnss.receiver.size(), gnss.transmitter.size(), comm);
    normalEquationInfo.analysisType = analysisType;
    initParameter(normalEquationInfo);
    if(!fileNameNormalsInfo.empty() && Parallel::isMaster(comm))
    {
      logStatus<<"write parameter name file <"<<fileNameNormalsInfo(fileNameVariableList)<<">"<<Log::endl;
      writeFileParameterName(fileNameNormalsInfo(fileNameVariableList), normalEquationInfo.parameterNames());
    }

    // =======================================================

    // Processing steps
    // ----------------
    Vector estimableReceiver(gnss.receiver.size());
    for(UInt idRecv=0; idRecv<gnss.receiver.size(); idRecv++)
      if(gnss.receiver.at(idRecv)->useable())
        estimableReceiver(idRecv) = TRUE;
    Parallel::reduceSum(estimableReceiver, 0, comm);
    Parallel::broadCast(estimableReceiver, 0, comm);
    for(UInt idRecv=0; idRecv<gnss.receiver.size(); idRecv++)
      if(estimableReceiver(idRecv))
        estimableReceiver(idRecv) = 1.;

    Bool changedNormalEquationInfo = FALSE;
    for(const auto &processingStep : processingSteps)
      processingStep->process(*this, estimableReceiver, normalEquationInfo, changedNormalEquationInfo);
  }
  catch(std::exception &e)
  {
    if(parallelIntervals)
    {
      logError<<timesInterval.at(idInterval).dateTimeStr()<<": Cannot solve (continue): "<<e.what()<<Log::endl;
      return;
    }
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssProcessing::initParameter(Gnss::NormalEquationInfo &normalEquationInfo)
{
  try
  {
    gnss.initParameter(normalEquationInfo);
    logInfo<<"  parameter count:        "<<normalEquationInfo.parameterCount()<<Log::endl;
    logInfo<<"  - epoch     parameters: "<<normalEquationInfo.blockIndex(normalEquationInfo.blockInterval())<<Log::endl;
    logInfo<<"  - interval  parameters: "<<normalEquationInfo.blockIndex(normalEquationInfo.blockAmbiguity())-normalEquationInfo.blockIndex(normalEquationInfo.blockInterval())<<Log::endl;
    logInfo<<"  - ambiguity parameters: "<<normalEquationInfo.parameterCount()-normalEquationInfo.blockIndex(normalEquationInfo.blockAmbiguity())<<Log::endl;
    logInfo<<"  block count:            "<<normalEquationInfo.blockCount()<<Log::endl;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssProcessing::regularizeNotUsedParameters(const Gnss::NormalEquationInfo &normalEquationInfo, UInt blockStart, UInt blockCount)
{
  try
  {
    UInt countZeros = 0;
    for(UInt i=blockStart; i<blockStart+blockCount; i++)
      if(normals.isMyRank(i,i))
      {
        Matrix &N = normals.N(i,i);
        for(UInt k=0; k<N.rows(); k++)
          if(N(k,k)==0.)
          {
            N(k,k) += 1.0;
            countZeros++;
            logWarning<<"    "<<normalEquationInfo.parameterNames().at(normals.blockIndex(i)+k).str()<<" has zero diagonal element"<<Log::endl;
          }
      }
    Parallel::reduceSum(countZeros, 0, normals.communicator());
    if(countZeros && Parallel::isMaster(normals.communicator()))
    {
      logWarning<<"  "<<countZeros<<" parameters have zero diagonal elements -> set to one"<<Log::endl;
      obsCount -= countZeros;
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssProcessing::collectNormalsBlocks(UInt blockStart, UInt blockCount)
{
  try
  {
    if(Parallel::size(normals.communicator()) <= 1)
      return;

    // synchronize used blocks
    for(UInt idBlock=blockStart; idBlock<blockStart+blockCount; idBlock++)
      for(UInt idProcess=0; idProcess<Parallel::size(normals.communicator()); idProcess++)
      {
        std::vector<UInt> columns; columns.reserve(normals._column.at(idBlock).size());
        std::vector<UInt> ranks;   ranks.reserve(normals._column.at(idBlock).size());
        normals.loopBlockRow(idBlock, {idBlock, normals.blockCount()}, [&](UInt k, UInt ik)
        {
          columns.push_back(k);
          ranks.push_back(normals._rank[ik]);
        });
        Parallel::broadCast(columns, idProcess, normals.communicator());
        Parallel::broadCast(ranks,   idProcess, normals.communicator());
        for(UInt ik=0; ik<columns.size(); ik++)
          normals.setBlock(idBlock, columns[ik], ranks[ik]);
      }

    // reduce sum
    for(UInt idBlock=blockStart; idBlock<blockStart+blockCount; idBlock++)
    {
      normals.loopBlockRow(idBlock, {idBlock, normals.blockCount()}, [&](UInt k, UInt /*ik*/) {normals.reduceSum(idBlock, k);});
      Parallel::reduceSum(n.at(idBlock), 0, normals.comm);
      if(!Parallel::isMaster(normals.comm))
        n.at(idBlock).setNull();
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssProcessing::buildNormals(Gnss::NormalEquationInfo &normalEquationInfo, Bool constraintsOnly, Bool solveEpochParameters, Bool removeEpochParameters)
{
  try
  {
    normals.initEmpty(normalEquationInfo.blockIndices(), normalEquationInfo.comm);
    n.resize(normals.blockCount());
    for(UInt i=0; i<normals.blockCount(); i++)
      n.at(i) = Vector(normals.blockSize(i));
    lPl = Vector(1);
    obsCount = 0;

    // Loop over all epochs
    // --------------------
    Parallel::barrier(normalEquationInfo.comm);
    logStatus<<"accumulate normals"<<Log::endl;
    if(normalEquationInfo.estimationType & Gnss::NormalEquationInfo::MASK_CONSTRAINT & Gnss::NormalEquationInfo::MASK_EPOCH)
    {
      logStatus<<"- additional epoch observations and constraints"<<Log::endl;
      UInt idLoop = 0;
      logTimerStart;
      for(UInt idEpoch : normalEquationInfo.idEpochs)
      {
        logTimerLoop(idLoop++, normalEquationInfo.idEpochs.size());
        gnss.observationEquationEpoch(normalEquationInfo, idEpoch, normals, n, lPl(0), obsCount);
      } // for(idEpoch)
      Parallel::barrier(normalEquationInfo.comm);
      logTimerLoopEnd(normalEquationInfo.idEpochs.size());
    }

    UInt blockStart = 0; // first block, which is not regularized and reduced
    if(!constraintsOnly && (normalEquationInfo.estimationType & ~Gnss::NormalEquationInfo::MASK_CONSTRAINT))
    {
      logStatus<<"- observation equations"<<Log::endl;
      Gnss::DesignMatrix A(normalEquationInfo);
      UInt idLoop     = 0;
      UInt blockCount = 0;
      logTimerStart;
      std::vector<Gnss::ObservationEquation> eqns(gnss.transmitter.size());
      for(UInt idEpoch : normalEquationInfo.idEpochs)
      {
        logTimerLoop(idLoop++, normalEquationInfo.idEpochs.size());

        // loop over all receivers
        for(UInt idRecv=0; idRecv<gnss.receiver.size(); idRecv++)
          if(normalEquationInfo.estimateReceiver.at(idRecv) && gnss.receiver.at(idRecv)->isEpochEstimable(normalEquationInfo.analysisType, idEpoch))
          {
            // all observation equations for this epoch
            UInt countEqn = 0;
            for(UInt idTrans=0; idTrans<gnss.receiver.at(idRecv)->idTransmitterSize(idEpoch); idTrans++)
              if(gnss.basicObservationEquations(normalEquationInfo, idRecv, idTrans, idEpoch, eqns.at(countEqn)))
                gnssIonosphere->eliminateTecParameter(normalEquationInfo, eqns.at(countEqn++));

            if(!normalEquationInfo.accumulateEpochObservations)
            {
              for(UInt i=0; i<countEqn; i++)
              {
                A.init(eqns.at(i).l);
                gnss.designMatrix(normalEquationInfo, eqns.at(i), A);
                A.accumulateNormals(normals, n, lPl(0), obsCount);
              }
            }
            else
            {
              // copy all observations to a single vector
              Vector l(std::accumulate(eqns.begin(), eqns.end(), UInt(0), [](UInt count, auto &e) {return count+e.l.rows();}));
              UInt idx=0;
              for(UInt i=0; i<countEqn; i++)
              {
                copy(eqns.at(i).l, l.row(idx, eqns.at(i).l.rows()));
                idx += eqns.at(i).l.rows();
              }
              A.init(l);

              idx=0;
              for(UInt i=0; i<countEqn; i++)
              {
                gnss.designMatrix(normalEquationInfo, eqns.at(i), A.selectRows(idx, eqns.at(i).l.rows()));
                idx += eqns.at(i).l.rows();
              }

              A.selectRows(0, 0); // select all
              A.accumulateNormals(normals, n, lPl(0), obsCount);
            }
          }

        // perform following steps not every epoch
        blockCount += normalEquationInfo.blockCountEpoch(idEpoch);
        if((blockCount < normalEquationInfo.defaultBlockCountReduction) && (idEpoch != normalEquationInfo.idEpochs.back()))
          continue;

        collectNormalsBlocks(blockStart, blockCount);

        if(solveEpochParameters)
        {
          regularizeNotUsedParameters(normalEquationInfo, blockStart, blockCount);
          normals.cholesky(FALSE, blockStart, blockCount, FALSE);
          normals.triangularTransSolve(n, blockStart, blockCount, FALSE);
          if(Parallel::isMaster(normalEquationInfo.comm))
            for(UInt idBlock=blockStart; idBlock<blockStart+blockCount; idBlock++)
            {
              lPl(0)   -= quadsum(n.at(idBlock)); // lPl = lPl - n1' N1^(-1) n1
              obsCount -= normals.blockSize(idBlock);
            }

          // remove N12 (epoch <-> other parameters)
          for(UInt idBlock=blockStart; idBlock<blockStart+blockCount; idBlock++)
            normals.loopBlockRow(idBlock, {removeEpochParameters ? idBlock : normalEquationInfo.blockInterval(), normals.blockCount()}, [&](UInt /*k*/, UInt ik)
            {
              if(normals._N[ik].size())
                normals._N[ik] = Matrix();
            });
        }

        blockStart += blockCount;
        blockCount  = 0;
      } // for(idEpoch)
      Parallel::barrier(normalEquationInfo.comm);
      logTimerLoopEnd(normalEquationInfo.idEpochs.size());
    } // if(!constraintsOnly)

    // other observations and constraints
    // ----------------------------------
    if(normalEquationInfo.estimationType & Gnss::NormalEquationInfo::MASK_CONSTRAINT & ~Gnss::NormalEquationInfo::MASK_EPOCH)
      gnss.observationEquation(normalEquationInfo, normals, n, lPl(0), obsCount);

    if(removeEpochParameters)
    {
      const UInt countBlocksEpochs = normalEquationInfo.blockInterval();
      normals.eraseBlocks(0, countBlocksEpochs);
      n.erase(n.begin(), n.begin()+countBlocksEpochs);
      normalEquationInfo.removeEpochParameters();
      blockStart = 0;
    }

    // collect normal equations
    // ------------------------
    if(Parallel::size(normals.comm)>1)
    {
      logStatus<<"collect normal equations"<<Log::endl;
      collectNormalsBlocks(blockStart, normals.blockCount()-blockStart);
      Parallel::reduceSum(lPl,      0, normalEquationInfo.comm);
      Parallel::reduceSum(obsCount, 0, normalEquationInfo.comm);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Double GnssProcessing::estimateSolution(Gnss::NormalEquationInfo &normalEquationInfo,
                                     Bool resolveAmbiguities, Bool dryRun, Bool computeResiduals,
                                     Gnss::Receiver::WeightingType computeWeights, Gnss::Receiver::WeightingType adjustSigma0)
{
  try
  {
    // setup normal equations
    const Bool solveEpochParameters = !normalEquationInfo.keepEpochNormalsInMemory && (normalEquationInfo.blockInterval() < normalEquationInfo.blockCount());
    buildNormals(normalEquationInfo, FALSE/*constraintsOnly*/, solveEpochParameters, FALSE/*removeEpochParameters*/);

    // eliminate all other parameters from the ambiguity normals
    // ---------------------------------------------------------
    Parallel::barrier(normalEquationInfo.comm);
    logStatus<<"solve"<<Log::endl;
    const UInt blockStart = (solveEpochParameters) ? normalEquationInfo.blockInterval() : 0; // epoch parameters already eliminated?
    const UInt blockCount = ((resolveAmbiguities) ? normalEquationInfo.blockAmbiguity() : normals.blockCount()) - blockStart;
    regularizeNotUsedParameters(normalEquationInfo, blockStart, normals.blockCount()-blockStart);
    normals.cholesky(TRUE/*timing*/, blockStart, blockCount, TRUE/*collect*/);
    normals.triangularTransSolve(n, blockStart, blockCount, TRUE/*collect*/); // forward step
    if(Parallel::isMaster(normalEquationInfo.comm))
      for(UInt z=blockStart; z<blockStart+blockCount; z++)
      {
        lPl(0)   -= quadsum(n.at(z)); // lPl = lPl - n1' N1^(-1) n1
        obsCount -= normals.blockSize(z);
      }

    // resolve integer ambiguities
    // ---------------------------
    if(resolveAmbiguities)
    {
      logStatus<<"Resolve integer ambiguities (may take a while)"<<Log::endl;
      Double sigmaFloat = gnss.ambiguityResolve(normalEquationInfo, normals, n, lPl(0), obsCount, dryRun);
      logInfo<<"  sigma(float) = "<<sigmaFloat%"%.2f"s<<Log::endl;
      logInfo<<"  sigma(fixed) = "<<sqrt(lPl(0)/obsCount)%"%.2f"s<<Log::endl;
      Parallel::barrier(normalEquationInfo.comm);
    } // if(parameterCountAmbiguities)

    // Compute sigma
    // -------------
    if(Parallel::isMaster(normalEquationInfo.comm))
    {
      const Double sigma = sqrt(lPl(0)/obsCount);
      logInfo<<"  sigma = "<<sigma%"%.2f"s<<Log::endl;
      if((sigma!=sigma) || (sigma<=0))
        logWarning<<"  Cannot compute sigma = sqrt("<<lPl(0)<<"/"<<obsCount<<")"<<Log::endl;
    }

    std::vector<Matrix> monteCarlo; // Monte Carlo Vector for redundancy computation

    // ========================================================================================

    // solve (Backward step)
    // ---------------------
    if(!solveEpochParameters)
    {
      normals.triangularSolve(n);
      if(computeResiduals)
      {
        monteCarlo.resize(normals.blockCount());
        for(UInt i=0; i<normals.blockCount(); i++)
          monteCarlo.at(i) = Matrix(normals.blockSize(i), 100);
        if(Parallel::isMaster(normalEquationInfo.comm))
          for(UInt i=0; i<blockCount; i++) // possible without resolved ambiguity parameters
            monteCarlo.at(i) = Vce::monteCarlo(monteCarlo.at(i).rows(), monteCarlo.at(i).columns());
        normals.triangularSolve(monteCarlo, 0, blockCount);
      }
    }

    // ========================================================================================

    // Reconstruct epoch parameters
    // ----------------------------
    if(solveEpochParameters && (normalEquationInfo.estimationType & ~Gnss::NormalEquationInfo::MASK_CONSTRAINT))
    {
      logStatus<<"Reconstruct epoch parameters"<<Log::endl;

      // solve other parameters (Backward step)
      // --------------------------------------
      normals.triangularSolve(n, blockStart, normals.blockCount()-blockStart);

      if(computeResiduals)
      {
        monteCarlo.resize(normals.blockCount());
        for(UInt i=blockStart; i<normals.blockCount(); i++)
          monteCarlo.at(i) = Matrix(normals.blockSize(i), 100);
        if(Parallel::isMaster(normalEquationInfo.comm))
          for(UInt i=blockStart; i<blockStart+blockCount; i++) // possible without resolved ambiguity parameters
            monteCarlo.at(i) = Vce::monteCarlo(monteCarlo.at(i).rows(), monteCarlo.at(i).columns());
        normals.triangularSolve(monteCarlo, blockStart, blockCount);
      }

      // free N22
      // --------
      for(UInt idBlock=blockStart; idBlock<normals.blockCount(); idBlock++)
        normals.loopBlockRow(idBlock, {idBlock, normals.blockCount()}, [&](UInt /*k*/, UInt ik)
        {
          if(normals._N[ik].size())
            normals._N[ik] = Matrix();
        });

      // broadcast n, Wz
      // ---------------
      for(UInt i=0; i<blockStart; i++)
        n.at(i).setNull();
      for(UInt i=blockStart; i<normals.blockCount(); i++)
        Parallel::broadCast(n.at(i), 0, normalEquationInfo.comm);
      if(computeResiduals)
      {
        for(UInt i=0; i<blockStart; i++)
          monteCarlo.at(i) = Matrix(normals.blockSize(i), 100);
        for(UInt i=blockStart; i<blockStart+blockCount; i++)
          Parallel::broadCast(monteCarlo.at(i), 0, normalEquationInfo.comm);
      }

      // reconstruct N12
      // ---------------
      Gnss::DesignMatrix A(normalEquationInfo);
      UInt idLoop = 0;
      logTimerStart;
      for(UInt idEpoch : normalEquationInfo.idEpochs)
      {
        logTimerLoop(idLoop++, normalEquationInfo.idEpochs.size());

        // loop over all receivers
        Gnss::ObservationEquation eqn;
        for(UInt idRecv=0; idRecv<gnss.receiver.size(); idRecv++)
          if(normalEquationInfo.estimateReceiver.at(idRecv) && gnss.receiver.at(idRecv)->isEpochEstimable(normalEquationInfo.analysisType, idEpoch))
            for(UInt idTrans=0; idTrans<gnss.receiver.at(idRecv)->idTransmitterSize(idEpoch); idTrans++)
              if(gnss.basicObservationEquations(normalEquationInfo, idRecv, idTrans, idEpoch, eqn))
              {
                gnssIonosphere->eliminateTecParameter(normalEquationInfo, eqn);
                A.init(eqn.l);
                gnss.designMatrix(normalEquationInfo, eqn, A);
                A.transMult(eqn.l-A.mult(n, blockStart, normals.blockCount()-blockStart), n, 0, blockStart);
                A.transMult(-A.mult(monteCarlo, blockStart, blockCount), monteCarlo, 0, blockStart);
              }
      } // for(idEpoch)
      Parallel::barrier(normalEquationInfo.comm);
      logTimerLoopEnd(normalEquationInfo.idEpochs.size());

      // collect
      // -------
      for(UInt i=0; i<blockStart; i++)
      {
        Parallel::reduceSum(n.at(i), 0, normalEquationInfo.comm);
        if(!Parallel::isMaster(normalEquationInfo.comm))
          n.at(i).setNull();
      }

      if(computeResiduals)
        for(UInt i=0; i<blockStart; i++)
        {
          Parallel::reduceSum(monteCarlo.at(i), 0, normalEquationInfo.comm);
          if(!Parallel::isMaster(normalEquationInfo.comm))
            monteCarlo.at(i).setNull();
        }

      // solve epoch (Forward step)
      // --------------------------
      normals.triangularTransSolve(n, 0, blockStart, FALSE);
      if(computeResiduals)
      {
        normals.triangularTransSolve(monteCarlo, 0, blockStart, FALSE);
        if(Parallel::isMaster(normalEquationInfo.comm))
          for(UInt i=0; i<blockStart; i++)
            axpy(1, Vce::monteCarlo(monteCarlo.at(i).rows(), monteCarlo.at(i).columns()), monteCarlo.at(i));
      }

      // solve epoch (Backward step)
      // ---------------------------
      normals.triangularSolve(n, 0, blockStart);
      if(computeResiduals)
        normals.triangularSolve(monteCarlo, 0, blockStart);
    }

    // ========================================================================================

    // free memory
    normals.initEmpty(normals.blockIndex(), normals.communicator());

    // copy solution to one vector
    // ---------------------------
    Vector x;
    if(Parallel::isMaster(normalEquationInfo.comm))
    {
      x = Vector(normals.parameterCount());
      for(UInt i=0; i<normals.blockCount(); i++)
        if(normals.blockSize(i))
          copy(n.at(i), x.row(normals.blockIndex(i), normals.blockSize(i)));
    }
    n.clear();
    Parallel::broadCast(x, 0, normalEquationInfo.comm);
    if(!computeResiduals)
      return gnss.updateParameter(normalEquationInfo, x, Matrix(), TRUE);

    // copy monte carlo to one vector
    // ------------------------------
    Matrix Wz;
    if(Parallel::isMaster(normalEquationInfo.comm))
    {
      Wz = Matrix(normals.parameterCount(), monteCarlo.at(0).columns());
      for(UInt i=0; i<normals.blockCount(); i++)
        if(normals.blockSize(i))
          copy(monteCarlo.at(i), Wz.row(normals.blockIndex(i), normals.blockSize(i)));
    }
    monteCarlo.clear();
    Parallel::broadCast(Wz, 0, normalEquationInfo.comm);

    // Residual tracking
    // -----------------
    Double                   maxChangeTec = 0;
    std::string              infoMaxChangeTec;
    std::vector<GnssType>    maxResidualType = gnss.types(~GnssType::FREQ_NO);
    std::vector<Double>      maxResidual(maxResidualType.size(), 0);
    std::vector<std::string> maxResidualInfo(maxResidualType.size());

    Double minSTEC   =  1e+99;
    Double maxSTEC   = -1e+99;
    Double meanSTEC  = 0;
    Double stdSTEC   = 0;
    UInt   countSTEC = 0;

    Parallel::barrier(normalEquationInfo.comm);
    logStatus<<"Compute residuals"<<Log::endl;
    Gnss::DesignMatrix A(normalEquationInfo);
    UInt idLoop = 0;
    logTimerStart;
    for(UInt idEpoch : normalEquationInfo.idEpochs)
    {
      logTimerLoop(idLoop++, normalEquationInfo.idEpochs.size());

      // loop over all receivers
      for(UInt idRecv=0; idRecv<gnss.receiver.size(); idRecv++)
        if(normalEquationInfo.estimateReceiver.at(idRecv) && gnss.receiver.at(idRecv)->isEpochEstimable(normalEquationInfo.analysisType, idEpoch))
        {
          Gnss::ObservationEquation eqn;
          for(UInt idTrans=0; idTrans<gnss.receiver.at(idRecv)->idTransmitterSize(idEpoch); idTrans++)
            if(gnss.basicObservationEquations(normalEquationInfo, idRecv, idTrans, idEpoch, eqn))
            {
              // setup observation equations
              A.init(eqn.l);
              gnss.designMatrix(normalEquationInfo, eqn, A);
              Vector We  = eqn.l - A.mult(x); // decorrelated residuals
              Matrix AWz = A.mult(Wz);        // redundancies

              // estimate & reduce B parameters (ionosphere)
              // -------------------------------------------
              gnssIonosphere->updateAndEliminateTecParameter(normalEquationInfo, eqn, We, AWz, maxChangeTec, infoMaxChangeTec);

              // STEC statistics
              // ---------------
              const Double dSTEC = gnss.receiver.at(idRecv)->observation(idTrans, idEpoch)->dSTEC;
              minSTEC   = std::min(dSTEC, minSTEC);
              maxSTEC   = std::max(dSTEC, maxSTEC);
              meanSTEC += dSTEC;
              stdSTEC  += dSTEC*dSTEC;
              countSTEC++;

              // redundancies
              // ------------
              Vector r(We.rows());
              for(UInt i=0; i<We.rows(); i++)
                r(i) = 1. - quadsum(AWz.row(i));

              // find max. residual (for statistics)
              // -----------------------------------
              if(norm(eqn.sigma-eqn.sigma0) < 1e-8) // without outlier
                for(UInt k=0; k<eqn.types.size(); k++)
                {
                  const Double de = We(k)*eqn.sigma(k) - gnss.receiver.at(idRecv)->observation(idTrans, idEpoch)->at(eqn.types.at(k)).residuals;
                  const UInt idType = GnssType::index(maxResidualType, eqn.types.at(k));
                  if(fabs(de) > maxResidual.at(idType))
                  {
                    maxResidual.at(idType) = fabs(de);
                    maxResidualInfo.at(idType) = "  max. residual change ("+(maxResidualType.at(idType)+eqn.transmitter->PRN()).str()+", "
                                              + eqn.receiver->name()+", "+gnss.times.at(idEpoch).dateTimeStr()+") "+ (1e3*de)%"%7.1f mm"s;
                  }
                } // for(k)

              gnss.receiver.at(idRecv)->observation(idTrans, idEpoch)->setDecorrelatedResiduals(eqn.types, We, r);
            } // for(idTrans)
        } // for(idRecv)
    } // for(idEpoch)
    Parallel::barrier(normalEquationInfo.comm);
    logTimerLoopEnd(normalEquationInfo.idEpochs.size());

    // new weights
    // -----------
    for(UInt idRecv=0; idRecv<gnss.receiver.size(); idRecv++)
      if(normalEquationInfo.estimateReceiver.at(idRecv) && gnss.receiver.at(idRecv)->useable())
        gnss.receiver.at(idRecv)->computeWeightsFromResiduals(normalEquationInfo.idEpochs, computeWeights, adjustSigma0);

    // residual analysis
    // -----------------
    std::vector<GnssType> type = maxResidualType;
    std::vector<Double>   ePe(type.size(), 0), redundancy(type.size(), 0);
    std::vector<UInt>     obsCount(type.size(), 0), outlierCount(type.size(), 0);
    for(UInt idRecv=0; idRecv<gnss.receiver.size(); idRecv++)
      if(normalEquationInfo.estimateReceiver.at(idRecv) && gnss.receiver.at(idRecv)->useable())
        gnss.receiver.at(idRecv)->residualsStatistics(NULLINDEX, normalEquationInfo.idEpochs, type, ePe, redundancy, obsCount, outlierCount);
    for(UInt i=0; i<type.size(); i++)
    {
      Parallel::reduceSum(ePe.at(i),          0, normalEquationInfo.comm);
      Parallel::reduceSum(redundancy.at(i),   0, normalEquationInfo.comm);
      Parallel::reduceSum(obsCount.at(i),     0, normalEquationInfo.comm);
      Parallel::reduceSum(outlierCount.at(i), 0, normalEquationInfo.comm);

      if(Parallel::isMaster(normalEquationInfo.comm) && obsCount.at(i))
      {
        logInfo<<"  "<<type.at(i).str()
              <<": sigma0 = "    <<Vce::standardDeviation(ePe.at(i), redundancy.at(i), 2.5/*huber*/, 1.5/*huberPower*/)%"%4.2f"s
              <<", redundancy = "<<(redundancy.at(i)/(obsCount.at(i)-outlierCount.at(i)))%"%4.2f"s
              <<", count = "     <<obsCount.at(i)%"%7i"s
              <<", outliers = "  <<outlierCount.at(i)%"%6i"s<<" ("<<(100.*outlierCount.at(i)/obsCount.at(i))%"%4.2f"s<<" %)"
              <<Log::endl;
      }
    }

    for(UInt i=0; i<maxResidual.size(); i++)
      Gnss::checkMaxChange(maxResidual.at(i), maxResidualInfo.at(i), TRUE/*printStatistics*/, normalEquationInfo.comm);

    // STEC analysis
    // -------------
    Parallel::reduceMin(minSTEC,   0, normalEquationInfo.comm);
    Parallel::reduceMax(maxSTEC,   0, normalEquationInfo.comm);
    Parallel::reduceSum(meanSTEC,  0, normalEquationInfo.comm);
    Parallel::reduceSum(stdSTEC,   0, normalEquationInfo.comm);
    Parallel::reduceSum(countSTEC, 0, normalEquationInfo.comm);
    stdSTEC   = sqrt((stdSTEC-meanSTEC*meanSTEC/countSTEC)/(countSTEC-1));
    meanSTEC /= countSTEC;
    logInfo<<"  STEC: mean = "<<meanSTEC%"%.2f +- "s<<stdSTEC%"%.2f (range "s<<minSTEC%"%.2f -- "s<<maxSTEC%"%.2f) TECU"s<<Log::endl;

    Gnss::checkMaxChange(maxChangeTec, infoMaxChangeTec, TRUE/*printStatistics*/, normalEquationInfo.comm);
    return gnss.updateParameter(normalEquationInfo, x, Wz, TRUE);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
