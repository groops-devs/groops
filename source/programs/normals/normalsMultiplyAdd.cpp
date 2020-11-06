/***********************************************/
/**
* @file normalsMultiplyAdd.cpp
*
* @brief n += alpha*N*x.
*
* @author Torsten Mayer-Guerr
* @date 2011-02-19
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program modifies \configFile{inputfileNormalEquation}{normalEquation} in a way
that $\bar{\M x}$ is estimated instead of $\M x$.
\begin{equation}
 \bar{\M x} := \M x + \alpha\, \M x_0,
\end{equation}
where $\M x_0$ is \configFile{inputfileParameter}{matrix} and $\alpha$ is \config{factor}.
This can be used to re-add reduced reference fields before a combined estimation
at normal equation level.
Therefore the right hand side of the normal equations is modified by
\begin{equation}
 \bar{\M n} := \M n + \alpha\,\M N\M x_0,
\end{equation}
and the quadratic sum of observations by
\begin{equation}
 \bar{\M l^T\M P\M l} := \M l^T\M P\M l + \alpha^2\,\M x_0^T\M N\M x_0 + 2\alpha\,\M x_0^T\M n
\end{equation}

As the normal matrix itself is not modified, rewriting of the matrix can be disabled by setting
\config{writeNormalMatrix} to false.
)";

/***********************************************/

#include "programs/program.h"
#include "parallel/matrixDistributed.h"
#include "files/fileMatrix.h"
#include "files/fileNormalEquation.h"

/***** CLASS ***********************************/

/** @brief n += alpha*N*x.
* @ingroup programsGroup */
class NormalsMultiplyAdd
{
public:
  void run(Config &config);
};

GROOPS_REGISTER_PROGRAM(NormalsMultiplyAdd, PARALLEL, "n += alpha*N*x", NormalEquation)

/***********************************************/

void NormalsMultiplyAdd::run(Config &config)
{
  try
  {
    FileName outName, normalsName, xName;
    Double   factor;
    Bool     writeMatrix;

    renameDeprecatedConfig(config, "outputfileNormalequation", "outputfileNormalEquation", date2time(2020, 6, 3));
    renameDeprecatedConfig(config, "inputfileNormalequation",  "inputfileNormalEquation",  date2time(2020, 6, 3));

    readConfig(config, "outputfileNormalEquation", outName,     Config::MUSTSET,  "",    "");
    readConfig(config, "inputfileNormalEquation",  normalsName, Config::MUSTSET,  "",    "");
    readConfig(config, "inputfileParameter",       xName,       Config::MUSTSET,  "",    "x");
    readConfig(config, "factor",                   factor,      Config::DEFAULT,  "1.0", "alpha");
    readConfig(config, "writeNormalMatrix",        writeMatrix, Config::DEFAULT,  "1",   "write full coefficient matrix, right hand sides and info files");
    if(isCreateSchema(config)) return;

    // ==================================

    logStatus<<"init normal equations"<<Log::endl;
    MatrixDistributed normal;
    Matrix n;
    NormalEquationInfo info;
    readFileNormalEquation(normalsName, info, normal, n);
    logInfo<<"  number of parameters:       "<<normal.parameterCount()<<Log::endl;
    logInfo<<"  number of right hand sides: "<<info.lPl.rows()<<Log::endl;

    Matrix x;
    if(Parallel::isMaster())
    {
      readFileMatrix(xName, x);
      if((x.rows() != n.rows()) || (x.columns() != n.columns()))
      {
        std::stringstream ss;
        ss<<"Dimension error N("<<n.rows()<<" x "<<n.rows()<<")*x("<<x.rows()<<" x "<<x.columns()<<") = n("<<n.rows()<<" x "<<n.columns()<<")";
        throw(Exception(ss.str()));
      }
    }

    // ==================================

    logStatus<<"multiply normal matrix"<<Log::endl;
    Parallel::broadCast(x);
    x *= factor;
    Matrix Nx(x.rows(), x.columns());
    for(UInt i=0; i<normal.blockCount(); i++)
    {
      if(normal.isMyRank(i,i))
        matMult(1., normal.N(i,i), x.row(normal.blockIndex(i), normal.blockSize(i)), Nx.row(normal.blockIndex(i), normal.blockSize(i)));
      for(UInt k=i+1; k<normal.blockCount(); k++)
        if(normal.isMyRank(i,k))
        {
          matMult(1., normal.N(i,k),         x.row(normal.blockIndex(k), normal.blockSize(k)), Nx.row(normal.blockIndex(i), normal.blockSize(i)));
          matMult(1., normal.N(i,k).trans(), x.row(normal.blockIndex(i), normal.blockSize(i)), Nx.row(normal.blockIndex(k), normal.blockSize(k)));
        }
    }
    Parallel::reduceSum(Nx);

    if(Parallel::isMaster())
    {
      for(UInt i=0; i<info.lPl.rows(); i++)
        info.lPl(i) += 2.*inner(x.column(i), n.column(i)) + inner(x.column(i), Nx.column(i));
      n += Nx;
    }

    // ==================================

    // write normal equations
    // ----------------------
    logStatus<<"write normal equations to <"<<outName<<">"<<Log::endl;
    if( (normalsName==outName) || !writeMatrix )
    {
      logStatus<<"(only the information file and the right hand sides)"<<Log::endl;
      if(Parallel::isMaster())
        writeFileNormalEquation(outName, info, n);
    }
    else
      writeFileNormalEquation(outName, info, normal, n);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
