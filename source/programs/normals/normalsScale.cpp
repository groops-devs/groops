/***********************************************/
/**
* @file normalsScale.cpp
*
* @brief Scales rows and columns of a normal equation system.
*
* @author Torsten Mayer-Guerr
* @date 2011-02-23
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
Scales rows and columns of a system of \configFile{inputfileNormalEquation}{normalEquation}
given by a diagonal matrix \configFile{inputfileFactorVector}{matrix} $\M S$
\begin{equation}
  \bar{\M N} := \M S \M N \M S \qquad\text{and}\qquad \bar{\M n} := \M S \M n.
\end{equation}
The estimated solution is now
\begin{equation}
  \bar{\M x} := \M S^{-1} \M x.
\end{equation}
This is effectively the same as rescaling columns of the design matrix.
This program is useful when combining normal equations from different sources,
for example in case the units of certain parameters don't match.
)";

/***********************************************/

#include "programs/program.h"
#include "parallel/matrixDistributed.h"
#include "files/fileMatrix.h"
#include "files/fileNormalEquation.h"

/***** CLASS ***********************************/

/** @brief Scales rows and columns of a normal equation system.
* @ingroup programsGroup */
class NormalsScale
{
public:
  void run(Config &config);
};

GROOPS_REGISTER_PROGRAM(NormalsScale, PARALLEL, "Scales rows and columns of a normal equation system", NormalEquation)

/***********************************************/

void NormalsScale::run(Config &config)
{
  try
  {
    FileName outName, normalsName, factorName;

    renameDeprecatedConfig(config, "outputfileNormalequation", "outputfileNormalEquation", date2time(2020, 6, 3));
    renameDeprecatedConfig(config, "inputfileNormalequation",  "inputfileNormalEquation",  date2time(2020, 6, 3));

    readConfig(config, "outputfileNormalEquation", outName,     Config::MUSTSET, "", "");
    readConfig(config, "inputfileNormalEquation",  normalsName, Config::MUSTSET, "", "");
    readConfig(config, "inputfileFactorVector",    factorName,  Config::MUSTSET, "", "Vector containing the factors");
    if(isCreateSchema(config)) return;

    // ==================================

    logStatus<<"read normal equations <"<<outName<<">"<<Log::endl;
    MatrixDistributed normal;
    Matrix n;
    NormalEquationInfo info;
    readFileNormalEquation(normalsName, info, normal, n);
    logInfo<<"  number of parameters:       "<<normal.parameterCount()<<Log::endl;
    logInfo<<"  number of right hand sides: "<<info.lPl.rows()<<Log::endl;

    logStatus<<"read factor vector <"<<factorName<<">"<<Log::endl;
    Vector factor;
    readFileMatrix(factorName, factor);

    if(Parallel::isMaster() && (factor.rows() != n.rows()))
      throw(Exception("Dimension error factor("+factor.rows()%"%i) != n("s+n.rows()%"%i)"s));

    // ==================================

    logStatus<<"scale normal matrix"<<Log::endl;
    for(UInt i=0; i<normal.blockCount(); i++)
      for(UInt k=i; k<normal.blockCount(); k++)
        if(normal.isMyRank(i,k))
        {
          Matrix &N = normal.N(i,k);
          for(UInt z=0; z<N.rows(); z++)
            N.row(z) *= factor(z+normal.blockIndex(i));
          for(UInt s=0; s<N.columns(); s++)
            N.column(s) *= factor(s+normal.blockIndex(k));
        }

    if(Parallel::isMaster())
      for(UInt z=0; z<n.rows(); z++)
        n.row(z) *= factor(z)*factor(z);

    // ==================================

    // write normal equations
    // ----------------------
    logStatus<<"write normal equations to <"<<outName<<">"<<Log::endl;
    writeFileNormalEquation(outName, info, normal, n);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
