/***********************************************/
/**
* @file gnssProcessingStep.h
*
* @brief Processing steps for GNSS normal equations.
*
* @author Torsten Mayer-Guerr
* @date 2021-09-05
*
*/
/***********************************************/

#ifndef __GROOPS_GNSSPROCESSINGSTEP__
#define __GROOPS_GNSSPROCESSINGSTEP__

// Latex documentation
#ifdef DOCSTRING_GnssProcessingStep
static const char *docstringGnssProcessingStep = R"(
\section{GnssProcessingStep}\label{gnssProcessingStepType}
Processing step in \program{GnssProcessing}.

Processing steps enable a dynamic definition of the consecutive steps performed during any kind of GNSS processing.
The most common steps are \configClass{estimate}{gnssProcessingStepType:estimate}, which performs an iterative least
squares adjustment, and \configClass{writeResults}{gnssProcessingStepType:writeResults}, which writes all output files
defined in \program{GnssProcessing} and is usually the last step.
Some steps such as \configClass{selectParametrizations}{gnssProcessingStepType:selectParametrizations},
\configClass{selectEpochs}{gnssProcessingStepType:selectEpochs},
\configClass{selectNormalsBlockStructure}{gnssProcessingStepType:selectNormalsBlockStructure}, and
\configClass{selectReceivers}{gnssProcessingStepType:selectReceivers} affect all subsequent steps.
In case these steps are used within a \configClass{group}{gnssProcessingStepType:group} or
\configClass{forEachReceiverSeparately}{gnssProcessingStepType:forEachReceiverSeparately} step,
they only affect the steps within this level.

For usage examples see cookbooks on \reference{GNSS satellite orbit determination and network analysis}{cookbook.gnssNetwork:processing}
or \reference{Kinematic orbit determination of LEO satellites}{cookbook.kinematicOrbit}.
)";
#endif

/***********************************************/

#include "base/import.h"
#include "parallel/matrixDistributed.h"
#include "gnss/gnss.h"

/**
* @defgroup gnssProcessingStepGroup GnssProcessingStep
* @brief Processing steps in @ref GnssProcessing.
* @ingroup gnssGroup
* The interface is given by @ref GnssProcessingStep.
* An Instance can be created by @ref readConfig. */
/// @{

/***** TYPES ***********************************/

class Gnss;
class GnssProcessingStep;
class GnssProcessingStepBase;
typedef std::shared_ptr<GnssProcessingStep> GnssProcessingStepPtr;

/***** CLASS ***********************************/

/** @brief Provides a list of receivers.
* An Instance of this class can be created by @ref readConfig. */
class GnssProcessingStep
{
  std::vector<GnssProcessingStepBase*> bases;

public:
  class State
  {
  public:
    GnssPtr                            gnss;
    GnssNormalEquationInfo             normalEquationInfo;
    Bool                               changedNormalEquationInfo;
    MatrixDistributed                  normals;
    std::vector<Matrix>                n;        // at master (after solve)
    Vector                             lPl;      // at master (after solve)
    UInt                               obsCount; // at master (after solve)
    std::vector<std::vector<GnssType>> sigmaType;
    std::vector<std::vector<Double>>   sigmaFactor; // for each receiver and type

    /** @brief Constructor. */
    State(GnssPtr gnss, Parallel::CommunicatorPtr comm);

    void regularizeNotUsedParameters(UInt blockStart, UInt blockCount);
    void collectNormalsBlocks       (UInt blockStart, UInt blockCount);
    void buildNormals               (Bool constraintsOnly, Bool solveEpochParameters);
    Double estimateSolution         (const std::function<Vector(const_MatrixSliceRef xFloat, MatrixSliceRef W, const_MatrixSliceRef d, Vector &xInt, Double &sigma)> &searchInteger,
                                     Bool computeResiduals,  Bool computeWeights, Bool adjustSigma0, Double huber, Double huberPower);
    void residualsStatistics        (UInt idRecv, UInt idTrans,
                                     std::vector<GnssType> &types, std::vector<Double> &ePe, std::vector<Double> &redundancy,
                                     std::vector<UInt> &obsCount, std::vector<UInt> &outlierCount);
  };

  /** @brief Constructor from config. */
  GnssProcessingStep(Config &config, const std::string &name);

  /// Destructor.
 ~GnssProcessingStep();

  /** @brief Perform the processing steps. */
  void process(State &state);

  /** @brief creates an derived instance of this class. */
  static GnssProcessingStepPtr create(Config &config, const std::string &name) {return GnssProcessingStepPtr(new GnssProcessingStep(config, name));}
};

/***** FUNCTIONS *******************************/

/** @brief Creates an instance of the class GnssProcessingStep.
* Search for a node with @a name in the Config node.
* if @a name is not found the function returns FALSE and @a var is untouched.
* @param config The config node which includes the node with the options for this class
* @param name Tag name in the config.
* @param[out] var Created class.
* @param mustSet If is MUSTSET and @a name is not found, this function throws an exception instead of returning with FALSE.
* @param defaultValue Ignored at the moment.
* @param annotation Description of the function of this class.
* @relates GnssProcessingStep */
template<> Bool readConfig(Config &config, const std::string &name, GnssProcessingStepPtr &var, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation);

/// @}

/***** CLASS ***********************************/

// Internal class
class GnssProcessingStepBase
{
public:
  virtual ~GnssProcessingStepBase() {}

  /** @brief Execute the processing step. */
  virtual void process(GnssProcessingStep::State &state) = 0;
  virtual Bool expectInitializedParameters() const {return TRUE;}
};

/***********************************************/

#endif /* __GROOPS___ */
