/***********************************************/
/**
* @file slrProcessingStep.h
*
* @brief Processing steps for SLR normal equations.
*
* @author Torsten Mayer-Guerr
* @date 2022-04-28
*
*/
/***********************************************/

#ifndef __GROOPS_SLRPROCESSINGSTEP__
#define __GROOPS_SLRPROCESSINGSTEP__

// Latex documentation
#ifdef DOCSTRING_SlrProcessingStep
static const char *docstringSlrProcessingStep = R"(
\section{SlrProcessingStep}\label{slrProcessingStepType}
Processing step in \program{SlrProcessing}.

Processing steps enable a dynamic definition of the consecutive steps performed during any kind of SLR processing.
The most common steps are \configClass{estimate}{slrProcessingStepType:estimate}, which performs an iterative least
squares adjustment, and \configClass{writeResults}{slrProcessingStepType:writeResults}, which writes all output files
defined in \program{SlrProcessing} and is usually the last step.
Some steps such as \configClass{selectParametrizations}{slrProcessingStepType:selectParametrizations}
and \configClass{selectStations}{slrProcessingStepType:selectStations} affect all subsequent steps.
In case these steps are used within a \configClass{group}{slrProcessingStepType:group} step,
they only affect the steps within this level.
)";
#endif

/***********************************************/

#include "base/import.h"
#include "parallel/matrixDistributed.h"
#include "slr/slr.h"

/**
* @defgroup slrProcessingStepGroup SlrProcessingStep
* @brief Processing steps in @ref SlrProcessing.
* @ingroup slrGroup
* The interface is given by @ref SlrProcessingStep. */
/// @{

/***** TYPES ***********************************/

class Slr;
class SlrProcessingStep;
class SlrProcessingStepBase;
typedef std::shared_ptr<SlrProcessingStep> SlrProcessingStepPtr;

/***** CLASS ***********************************/

/** @brief Provides a list of stations.
* An Instance of this class can be created by @ref readConfig. */
class SlrProcessingStep
{
  std::vector<SlrProcessingStepBase*> bases;

public:
  class State
  {
  public:
    SlrPtr                 slr;
    SlrNormalEquationInfo  normalEquationInfo;
    Bool                   changedNormalEquationInfo;
    MatrixDistributed      normals;
    std::vector<Matrix>    n;           // at master (after solve)
    Vector                 lPl;         // at master (after solve)
    UInt                   obsCount;    // at master (after solve)
    Vector                 sigmaFactor; // for each station

    /** @brief Constructor. */
    State(SlrPtr slr);

    void   regularizeNotUsedParameters(UInt blockStart, UInt blockCount);
    void   buildNormals(Bool constraintsOnly);
    Double estimateSolution(Bool computeResiduals,  Bool computeWeights, Bool adjustSigma0, Double huber, Double huberPower);
    void   residualsStatistics(UInt idStat, UInt idSat, Double &ePe, Double &redundancy, UInt &obsCount, UInt &outlierCount);
  };

  /** @brief Constructor from config. */
  SlrProcessingStep(Config &config, const std::string &name);

  /// Destructor.
 ~SlrProcessingStep();

  /** @brief Perform the processing steps. */
  void process(State &state);

  /** @brief creates an derived instance of this class. */
  static SlrProcessingStepPtr create(Config &config, const std::string &name) {return SlrProcessingStepPtr(new SlrProcessingStep(config, name));}
};

/***** FUNCTIONS *******************************/

/** @brief Creates an instance of the class SlrProcessingStep.
* Search for a node with @a name in the Config node.
* if @a name is not found the function returns FALSE and @a var is untouched.
* @param config The config node which includes the node with the options for this class
* @param name Tag name in the config.
* @param[out] var Created class.
* @param mustSet If is MUSTSET and @a name is not found, this function throws an exception instead of returning with FALSE.
* @param defaultValue Ignored at the moment.
* @param annotation Description of the function of this class.
* @relates SlrProcessingStep */
template<> Bool readConfig(Config &config, const std::string &name, SlrProcessingStepPtr &var, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation);

/// @}

/***** CLASS ***********************************/

// Internal class
class SlrProcessingStepBase
{
public:
  virtual ~SlrProcessingStepBase() {}

  /** @brief Execute the processing step. */
  virtual void process(SlrProcessingStep::State &state) = 0;
  virtual Bool expectInitializedParameters() const {return TRUE;}
};

/***********************************************/

#endif /* __GROOPS___ */
