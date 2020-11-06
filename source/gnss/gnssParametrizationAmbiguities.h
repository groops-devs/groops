/***********************************************/
/**
* @file gnssParametrizationAmbiguities.h
*
* @brief integer and float ambiguities.
*
* @author Torsten Mayer-Guerr
* @author Sebastian Strasser
* @date 2013-06-24
*
*/
/***********************************************/

#ifndef __GROOPS_GNSSPARAMETRIZATIONAMBIGUITIES__
#define __GROOPS_GNSSPARAMETRIZATIONAMBIGUITIES__

// Latex documentation
#ifdef DOCSTRING_GnssParametrizationAmbiguities
static const char *docstringGnssParametrizationAmbiguities = R"(
\section{GnssParametrizationAmbiguities}\label{gnssParametrizationAmbiguitiesType}

This class provides functionality for GNSS phase ambiguities and integer ambiguity resolution.

Integer ambiguity resolution is performed based on the least squares ambiguity decorrelation adjustment
(LAMBDA) method (Teunissen 1995, DOI \href{https://doi.org/10.1007/BF00863419}{10.1007/BF00863419}), specifically
the modified algorithm (MLAMBDA) by Chang et al. (2005, DOI \href{https://doi.org/10.1007/s00190-005-0004-x}{10.1007/s00190-005-0004-x}).
First the covariance matrix of the integer ambiguity parameters is computed by eliminating all but those parameters
from the full normal equation matrix and inverting it. Then, a Z-transformation is performed as described by
Chang et al. (2005) to decorrelate the ambiguity parameters without losing their integer nature.

The search process follows MLAMBDA and uses integer minimization of the weighted sum of squared residuals.
It is computationally infeasible to search a hyper-ellipsoid with a dimension of ten thousand or more.
Instead, a blocked search algorithm is performed by moving a window with a length of, for example,
\config{searchBlockSize}=\verb|200| parameters over the decorrelated ambiguities, starting from the most accurate.
In each step, the window is moved by half of its length and the overlapping parts are compared to each other.
If all fixed ambiguities in the overlap agree, the algorithm continues.
Otherwise, both windows are combined and the search is repeated using the combined window, again comparing with the overlapping
part of the preceding window. If not all solutions could be checked for a block after \config{maxSearchSteps},
the selected \config{incompleteAction} is performed.
If the algorithm reaches ambiguities with a standard deviation higher than \config{sigmaMaxResolve},
ambiguity resolution stops and the remaining ambiguities are left as float values.
Otherwise, all ambiguity parameters are fixed to integer values.

In contrast to an integer least squares solution over the full ambiguity vector, it is not guaranteed that the resulting solution
is optimal in the sense of minimal variance with given covariance.
This trade-off is necessary to cope with large numbers of ambiguities.

See also \program{GnssProcessing}.
)";
#endif

/***********************************************/

#include "config/config.h"
#include "gnss/gnss.h"
#include "gnss/gnssReceiver.h"

/** @addtogroup gnssGroup */
/// @{

/***** TYPES ***********************************/

class GnssParametrizationAmbiguities;
typedef std::shared_ptr<GnssParametrizationAmbiguities> GnssParametrizationAmbiguitiesPtr;

/***** CLASS ***********************************/

/** @brief Parametrization of ambiguities.
* An Instance of this class can be created by @ref readConfig. */
class GnssParametrizationAmbiguities : public Gnss::Parametrization
{
  class Ambiguity : public Gnss::Ambiguity
  {
    public:
    Gnss::Track           *track;
    std::vector<GnssType>  types;
    Matrix                 T;              // Matrix to transform ambiguities to observations [cycles -> m]
    Vector                 value;          // ambiguities in cycles
    Gnss::ParameterIndex   parameterIndex; // index in the parameter vector
    Vector                 resolved;
    Bool                   isInteger;

    Vector ambiguities(const std::vector<GnssType> &types) const override;
  };

  typedef std::shared_ptr<Ambiguity> AmbiguityPtr;

  class AmbiguityInfo; // to share information between processes

  // ======================

  // sparse transformation matrix (in both directions: direct and back)
  // created by elementary transformations (swap and reduce)
  class Transformation
  {
    std::vector<std::map<UInt, Double>> columnForward, columnBackward; // for each row a column with a factor

  public:
    Transformation(UInt dim=0);
   ~Transformation() {}

    void reduce(Double alpha, UInt i, UInt k); // row(k) -= alpha * row(i)
    void swap(UInt i, UInt k);

    Matrix transform(const_MatrixSliceRef x) const;
    Matrix transformBack(const_MatrixSliceRef x)  const;
    Matrix distributeBack(const_MatrixSliceRef x) const;
  };

  // ======================

  enum class IncompleteAction {STOP, SHRINKBLOCKSIZE, IGNORE, EXCEPTION};

  std::vector<AmbiguityPtr> ambiguities;
  Double                    sigmaMaxResolve;
  UInt                      searchBlockSize;
  UInt                      maxSearchSteps;
  IncompleteAction          incompleteAction;

  std::string              receiverName;
  FileName                 fileNameAmbiguities;
  std::vector<Matrix>      solutionSteps;  /// resolved ambiguities (for output); key = empty string or station name, value = iterations in searchIntegerBlocked()
  std::vector<std::string> solutionNames;

  /** @brief Decorrelate ambiguities (Melbourne Wuebbena like linear combinations).
  * @param types list of phase observations.
  * @return Transformation matrix from decorrelated ambiguities to phase observations [cycles]->[m]. */
  static Matrix phaseDecorrelation(const std::vector<GnssType> &types, Double wavelengthFactor);

  // LAMBDA method
  static void   choleskyReversePivot(Matrix &N, Transformation &Z);
  static void   choleskySwap(UInt i, Double delta, MatrixSliceRef W, std::vector<Double> &d, Transformation &transformation);
  static void   choleskyReduce(UInt i, UInt k, MatrixSliceRef W, Transformation &transformation);
  static Vector choleskyTransform(MatrixSliceRef W, Transformation &transformation);
  static Bool   searchInteger(const_MatrixSliceRef xFloat, MatrixSliceRef W, const_MatrixSliceRef d, UInt maxSearchSteps, Vector &solution, Double &minNorm);
  void   searchIntegerBlocked(const_MatrixSliceRef xFloat, MatrixSliceRef W, const_MatrixSliceRef d,
                              Vector &xInt, Vector &isNotFixed, Matrix &solutionSteps, Double &sigma) const;

friend class Gnss;
friend class GnssReceiver;

public:
  GnssParametrizationAmbiguities(Config &config, const std::string &name); /// Constructor.
  virtual ~GnssParametrizationAmbiguities() {}    /// Destructor.

  void   initIntervalLate(Gnss::AnalysisType analysisType, const std::vector<Time> &times, const Time &timeMargin, Parallel::CommunicatorPtr comm) override;
  void   initParameter(Gnss::NormalEquationInfo &normalEquationInfo) override;
  void   aprioriParameter(const Gnss::NormalEquationInfo &normalEquationInfo, MatrixSliceRef x0) const override;
  Bool   isDesignMatrix(const Gnss::NormalEquationInfo &normalEquationInfo, UInt idRecv, UInt idTrans, UInt idEpoch) const override;
  void   designMatrix(const Gnss::NormalEquationInfo &normalEquationInfo, const Gnss::ObservationEquation &eqn, Gnss::DesignMatrix &A) const override;
  Double updateParameter(const Gnss::NormalEquationInfo &normalEquationInfo, const_MatrixSliceRef x, const_MatrixSliceRef Wz, Bool printStatistics) override;
  Double ambiguityResolve(const Gnss::NormalEquationInfo &normalEquationInfo, MatrixDistributed &normals, std::vector<Matrix> &n, Double &lPl, UInt &obsCount, Bool dryRun=FALSE);
  void   writeResults(const Gnss::NormalEquationInfo &normalEquationInfo, const std::string &suffix) override;

  /** @brief creates an derived instance of this class. */
  static GnssParametrizationAmbiguitiesPtr create(Config &config, const std::string &name) {return std::make_shared<GnssParametrizationAmbiguities>(config, name);}
};

/***** FUNCTIONS *******************************/

/** @brief Creates an instance of the class GnssParametrizationAmbiguities.
* Search for a node with @a name in the Config node.
* if @a name is not found the function returns FALSE and @a var is untouched.
* @param config The config node which includes the node with the options for this class
* @param name Tag name in the config.
* @param[out] var Created class.
* @param mustSet If is MUSTSET and @a name is not found, this function throws an exception instead of returning with FALSE.
* @param defaultValue Ignored at the moment.
* @param annotation Description of the function of this class.
* @relates GnssParametrizationAmbiguities */
template<> Bool readConfig(Config &config, const std::string &name, GnssParametrizationAmbiguitiesPtr &var, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation);

/// @}

/***********************************************/

#endif /* __GROOPS___ */
