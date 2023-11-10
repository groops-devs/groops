/***********************************************/
/**
* @file matrixGenerator.h
*
* @brief Matrix calculation.
*
* @author Torsten Mayer-Guerr
* @date 2014-03-18
*
*/
/***********************************************/

#ifndef __GROOPS_MATRIXGENERATOR__
#define __GROOPS_MATRIXGENERATOR__

// Latex documentation
#ifdef DOCSTRING_MatrixGenerator
static const char *docstringMatrixGenerator = R"(
\section{MatrixGenerator}\label{matrixGeneratorType}
This class provides a matrix used e.g. by \program{MatrixCalculate}.
If multiple matrices are given the resulting matrix is the sum all
and the size is exandeded to fit all matrices. Before the computation of each submatrix
the variables \verb|rowsBefore| and \verb|columnsBefore| with current size of the overall matrix
are set. As all matrices can be manipulated before, complex matrix operations are possible.
)";
#endif

/***********************************************/

#include "base/import.h"
#include "config/config.h"

/**
* @defgroup matrixGenerator MatrixGenerator
* @brief Matrix calculation.
* @ingroup classesGroup
* The interface is given by @ref MatrixGenerator.
* An Instance can be created by @ref readConfig. */
/// @{

/***** TYPES ***********************************/

class MatrixGenerator;
class MatrixGeneratorBase;
typedef std::shared_ptr<MatrixGenerator> MatrixGeneratorPtr;

/***** CLASS ***********************************/

/** @brief Matrix calculation.
* This class provides a matrix.
* An Instance of this class can be created by @ref readConfig. */
class MatrixGenerator
{
  std::vector<MatrixGeneratorBase*> matrix;

public:
  /// Constructor.
  MatrixGenerator(Config &config, const std::string &name);

  /// Destructor.
  ~MatrixGenerator();

  /// provides a matrix.
  Matrix compute();

  /** @brief creates an derived instance of this class. */
  static MatrixGeneratorPtr create(Config &config, const std::string &name) {return MatrixGeneratorPtr(new MatrixGenerator(config, name));}
};


/***** FUNCTIONS *******************************/

/** @brief Creates an instance of the class MatrixGenerator.
* Search for a node with @a name in the Config node.
* if @a name is not found the function returns FALSE and @a matrixGenerator with zero-size matrix is created.
* @param config The config node which includes the node with the options for this class
* @param name Tag name in the config.
* @param[out] matrixGenerator Created class.
* @param mustSet If is MUSTSET and @a name is not found, this function throws an exception instead of returning with FALSE.
* @param defaultValue Ignored at the moment.
* @param annotation Description of the function of this class.
* @relates MatrixGenerator */
template<> Bool readConfig(Config &config, const std::string &name, MatrixGeneratorPtr &matrixGenerator, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation);

/// @}

/***** CLASS ***********************************/

// Internal class
class MatrixGeneratorBase
{
public:
  virtual ~MatrixGeneratorBase() {}
  virtual void compute(Matrix &A, UInt rowsBefore, UInt columnsBefore, UInt &startRow, UInt &startCol) = 0;
};

/***********************************************/

#endif
