/***********************************************/
/**
* @file normalEquation.h
*
* @brief Systems of normal equations.
* Creation, combination and solving of systems of normal equations.
*
* @author Torsten Mayer-Guerr
* @date 2004-12-10
*
*/
/***********************************************/

#ifndef __GROOPS_NORMALEQUATION__
#define __GROOPS_NORMALEQUATION__

// Latex documentation
#ifdef DOCSTRING_NormalEquation
static const char *docstringNormalEquation = R"(
\section{NormalEquation}\label{normalEquationType}
This class provides a system of normal equations.
This total system is the weighted sum of individual normals.
\begin{equation}
 \M N_{total} =  \sum_{k=1} \frac{1}{\sigma_k^2}\M N_k
 \qquad\text{and}\qquad
\M n_{total} = \sum_{k=1} \frac{1}{\sigma_k^2} \M n_k.
\end{equation}
The normals do not need to have the same dimension. The dimension
of the total combined system is chosen to cover all individual systems.
For each normal a \config{startIndex} is required which indicates
the position of the first unknown of the individual normal within the
combined parameter vector.

The $\sigma_k$ of the relative weights are defined by \config{aprioriSigma}
in a first step. If an apriori solution \configFile{inputfileApproxSolution}{matrix} is
given or the normals are solved iteratively the weights are determined by means
of variance compoment estimation (VCE), see \program{NormalsSolverVCE}:
\begin{equation}
\sigma_k^2 =
\frac{\M e_k^T\M P\M e_k}
{n_k-\frac{1}{\sigma_k^2}\text{trace}\left(\M N_k\M N_{total}^{-1}\right)},
\end{equation}
where $n_k$ is the number of observations. The square sum of the residuals
is calculated by
\begin{equation}
\M e_k^T\M P\M e_k = \M x^T\M N_k\M x - 2\M n_k^T\M x + \M l_k^T\M P_k\M l_k.
\end{equation}
The system of normal equations can be solved with several right hand sides at once. But
only one right hand side, which can be selected with the index \config{rightHandSide}
(counting from zero), can be used to compute the variance factors.
The combined normal $\M N_{total}$ and the solution $\M x$ are taken from the previous
iteration step. In case of \configClass{DesignVCE}{normalEquationType:designVCE} the algorithm
is a little bit different as described below.
)";
#endif

/***********************************************/

#include "base/import.h"
#include "base/parameterName.h"
#include "inputOutput/fileName.h"
#include "config/config.h"
#include "parallel/matrixDistributed.h"

/**
* @defgroup normalEquationGroup NormalEquation
* @brief Systems of normal equations.
* @ingroup classesGroup
* The interface is given by @ref NormalEquation.
* An Instance can be created by @ref readConfig. */
/// @{

/***** TYPES ***********************************/

class NormalEquation;
class NormalEquationBase;
typedef std::shared_ptr<NormalEquation> NormalEquationPtr;

/***** CLASS ***********************************/

/** @brief Systems of normal equations.
* Creation, combination and solving of systems of normal equations.
* An Instance of this class can be created by @ref readConfig. */
class NormalEquation
{
private:
  enum Status {UNKNOWN, INIT, NORMAL, CHOLESKY, INVERSECHOLESKY, INVERSE};
  Status            status;
  UInt              rhsNo;
  MatrixDistributed normals;
  Matrix            n;        // right hand sides (at master)
  Vector            lPl;      // Norm of the observations
  UInt              obsCount;
  Matrix            Wz;       // Monte-Carlo-vector
  Matrix            x;        // current solution

  std::vector<NormalEquationBase*> normalsComponent;

public:
  /// Constructor.
  NormalEquation(Config &config, const std::string &name);

  /// Destructor.
  virtual ~NormalEquation();

  /** @brief Init systems of normal equations.
  * @param blockSize normal matrix is divided into blocks, (0: only one block).
  * @param comm normal matrix is distributed over processes.  */
  void init(UInt blockSize, Parallel::CommunicatorPtr comm);

  /** @brief Number of unknown parameters.
  * Dimension of the normal matrix N. */
  UInt parameterCount() const {return normals.parameterCount();}

  /** @brief Count of the observation vectors on the right hand side.
  * Columns of the right hand side n. */
  UInt rightHandSideCount() const {return n.columns();}

  /** @brief Total count of observations.
  * This includes pseudo observations by the regularization. */
  UInt observationCount() const {return obsCount;}

  /** @brief Names of the unknown parameters. */
  std::vector<ParameterName> parameterNames() const;

  /** @brief Variance factors of the different normal equations systems.
  * Estimated with Variance Component Estimation (VCE). */
  Vector varianceComponentFactors() const;

  /** @brief Set approximate solution.
  * To speed up the iterative Variance Component Estimation (VCE). */
  void setApproximateSolution(const const_MatrixSlice &x0);

  /** @name Status
  * @brief The following functions changes the state of this class.
  * Optimal performance is only given if the functions are called in a specific state.
  * The expected state and the state after the call are given in brackets.
  * (expected state -> state after call). */
  //@{

  /** @brief Accumulate normal equations.
  * The weight of the different normal equation systems are computed
  * in terms of Variance Component Estimation (VCE).
  * This is an iterative process as it depends on residuals which depend on the solution.
  * If convergence is reached TRUE is returned.
  * Change of state of this class: (->NORMAL). */
  Bool build(UInt rightHandSide=0);

  /** @brief Write systems of normal equations to file.
  * Change of state of this class: (NORMAL -> NORMAL). */
  void write(const FileName &name);

  /** @brief Solve the system of normal equations.
  * Change of state of this class: (NORMAL -> CHOLESKY).
  * @return Solution vector as columns of the matrix. */
  Matrix solve();

  /** @brief A posteriori sigma.
  * Change of state of this class: (CHOLESKY -> CHOLESKY). */
  Double aposterioriSigma();

  /** @brief Inverse of the combined normal matrix.
  * The posteriori sigma is applied if possible.
  * Change of state of this class: (CHOLESKY -> CHOLESKYINVERSE).
  * @return Accuracy of the solution (Square root of the diagonals of the inverse). */
  Vector sigmaParameter();

  /** @brief Write variance/covariance matrix.
  * Change of state of this class: (CHOLESKYINVERSE -> INVERSE). */
  void writeCovariance(const FileName &name);

  /** @brief Contribution of normal equation components.
  * Each column of the returned matrix contains the contribution
  * of the normal equation component to the solution.
  * The sum of the elements in each row is one.
  * Change of state of this class: (INVERSE -> INVERSE).
  * @return contribution vectors as coulumns of a matrix */
  Matrix contribution();
  //@}

  /** @brief creates an derived instance of this class. */
  static NormalEquationPtr create(Config &config, const std::string &name) {return NormalEquationPtr(new NormalEquation(config, name));}
};

/***** FUNCTIONS *******************************/

/** @brief Creates an instance of the class NormalEquation.
* Search for a node with @a name in the Config node.
* if @a name is not found the function returns FALSE and class with empty normals is created.
* @param config The config node which includes the node with the options for this class
* @param name Tag name in the config.
* @param[out] normalEquation Created class.
* @param mustSet If is MUSTSET and @a name is not found, this function throws an exception instead of returning with FALSE.
* @param defaultValue Ignored at the moment.
* @param annotation Description of the function of this class.
* @relates NormalEquation */
template<> Bool readConfig(Config &config, const std::string &name, NormalEquationPtr &normalEquation, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation);

/// @}

/***** CLASS ***********************************/

// Internal class
class NormalEquationBase
{
public:
  virtual ~NormalEquationBase() {}

protected:
  virtual UInt   rightHandSideCount() const = 0;
  virtual UInt   parameterCount()     const = 0;
  virtual void   parameterNames(std::vector<ParameterName> &name) const = 0;
  virtual void   init(MatrixDistributed &normals, UInt rhsCount) = 0;
  virtual Bool   addNormalEquation(UInt rightHandSide, const const_MatrixSlice &x, const const_MatrixSlice &Wz,
                                   MatrixDistributed &normals, Matrix &n, Vector &lPl, UInt &obsCount) = 0;
  virtual Vector contribution(MatrixDistributed &Cov) = 0;
  virtual std::vector<Double> varianceComponentFactors() const = 0;

  friend class NormalEquation;
};

/***********************************************/

#endif /* __GROOPS_NORMALEQUATION__ */
