/***********************************************/
/**
* @file normalsCreate.cpp
*
* @brief Create normal equations from calulated matrices
*
* @author Torsten Mayer-Guerr
* @date 2017-09-04
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
Create \file{normal equations}{normalEquation}
from calulated matrices (\configClass{matrixGenerator}{matrixGeneratorType}).

The \configFile{inputfileParameterNames}{parameterName} can be created with \program{ParameterNamesCreate}.

The \configClass{normalMatrix}{matrixGeneratorType} must be symmetric.
The \configClass{rightHandSide}{matrixGeneratorType} must have the same number of rows
and can contain multiple columns for multiple solutions.

The Vector $\M l^T\M P\M l$ is the quadratic sum of observations for each column of the right hand side.
It is used to determine the aposteriori accuracy
\begin{equation}
\hat{\sigma}^2 = \frac{\hat{\M e}^T\M P\hat{\M e}}{n-m} = \frac{\M l^T\M P\M l - \M n^T\hat{\M x}}{n-m}.
\end{equation}
If the vector is not given, it is automatically determined by assuming $\hat{\sigma}^2=1$.

The number of observations~$n$ is given by the expression \config{observationCount}.
The variable \verb|observationCount| can be used, if it is set by a normal equation file
\configFile{inputfileNormalEquationObsCount}{normalEquation}.
)";

/***********************************************/

#include "programs/program.h"
#include "files/fileNormalEquation.h"
#include "files/fileParameterName.h"
#include "classes/matrixGenerator/matrixGenerator.h"

/***** CLASS ***********************************/

/** @brief Create normal equations from calulated matrices.
* @ingroup programsGroup */
class NormalsCreate
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(NormalsCreate, SINGLEPROCESS, "create normal equations from calulated matrices", NormalEquation)

/***********************************************/

void NormalsCreate::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName              fileNameNormals;
    FileName              fileNameParameterName;
    FileName              fileNameNormalsObsCount;
    MatrixGeneratorPtr    generate_N, generate_n, generate_lPl;
    ExpressionVariablePtr exprObsCount;

    renameDeprecatedConfig(config, "outputfileNormalequation",        "outputfileNormalEquation",        date2time(2020, 6, 3));
    renameDeprecatedConfig(config, "inputfileNormalequationObsCount", "inputfileNormalEquationObsCount", date2time(2020, 6, 3));

    readConfig(config, "outputfileNormalEquation",        fileNameNormals,         Config::MUSTSET,  "", "");
    readConfig(config, "inputfileParameterNames",         fileNameParameterName,   Config::OPTIONAL, "", "");
    readConfig(config, "normalMatrix",                    generate_N,              Config::MUSTSET,  "", "");
    readConfig(config, "rightHandSide",                   generate_n,              Config::MUSTSET,  "", "");
    readConfig(config, "lPl",                             generate_lPl,            Config::DEFAULT,  "", "vector with size of rhs columns");
    readConfig(config, "inputfileNormalEquationObsCount", fileNameNormalsObsCount, Config::OPTIONAL, "", "sets the variable observationCount");
    readConfig(config, "observationCount",                exprObsCount,            Config::MUSTSET,  "observationCount", "(variables: rows, columns (rhs), observationCount)");
    if(isCreateSchema(config)) return;

    // create matrices
    // ---------------
    logStatus<<"create matrices"<<Log::endl;
    Matrix N = generate_N->compute();
    Matrix n = generate_n->compute();
    NormalEquationInfo info(std::vector<ParameterName>(N.rows()), generate_lPl->compute());

    // checks
    if((N.rows() != N.columns()) || (N.getType() != Matrix::SYMMETRIC))
      throw(Exception("Normal matrix must quadratix and symmetric"));
    if((n.rows() != N.rows()) || (n.getType() != Matrix::GENERAL))
      throw(Exception("Right hand side: dimension error"));
    if(info.lPl.size() && (info.lPl.rows() != n.columns()))
      throw(Exception("lPl: dimension error"));

    // compute observation count
    // -------------------------
    auto varList = config.getVarList();
    addVariable("rows",    n.rows(),    varList);
    addVariable("columns", n.columns(), varList);
    if(!fileNameNormalsObsCount.empty())
    {
      Matrix n;
      NormalEquationInfo info;
      readFileNormalEquation(fileNameNormalsObsCount, info, n);
      addVariable("observationCount", info.observationCount, varList);
    }
    info.observationCount = static_cast<UInt>(exprObsCount->evaluate(varList));

    // simulate lPl
    // ------------
    if(info.lPl.size() == 0)
    {
      logStatus<<"Simulate lPl"<<Log::endl;
      Matrix x = n;
      Matrix N2 = N;
      for(UInt i=0; i<N2.rows(); i++)
        if(N2(i,i) == 0)
          N2(i,i) = 1.;
      solveInPlace(N2, x);
      // sigma^2 = (lPl-n'x)/(m-n)
      info.lPl = Vector(n.columns());
      for(UInt i=0; i<info.lPl.rows(); i++)
        info.lPl(i) = (info.observationCount-x.rows()) + inner(x.column(i), n.column(i));
    }

    // parameter names
    // ---------------
    if(!fileNameParameterName.empty())
      readFileParameterName(fileNameParameterName, info.parameterName);

    // write normal equations
    // ----------------------
    logStatus<<"write normal equations to <"<<fileNameNormals<<">"<<Log::endl;
    writeFileNormalEquation(fileNameNormals, info, N, n);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
