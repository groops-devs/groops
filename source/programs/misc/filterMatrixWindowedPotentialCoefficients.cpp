/***********************************************/
/**
* @file filterMatrixWindowedPotentialCoefficients.cpp
*
* @brief Create a spherical harmonic window matrix.
*
* @author Torsten Mayer-Guerr
* @author Andreas Kvas
* @date 2020-07-21
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
Create a spherical harmonic window matrix. The window matrix $\mathbf{W}$ is generated in space domain through
spherical harmonic synthesis and analysis matrices.
The resulting linear operator can be written as
\begin{equation}
\mathbf{W} = \mathbf{K} \mathbf{A} \mathbf{\Omega} \mathbf{S} \mathbf{K}^{-1}.
\end{equation}
Here, $\mathbf{K}$ is a diagonal matrix with the \configClass{kernel}{kernelType} coefficients on the main diagonal,
$\mathbf{S}$ is the spherical harmonic synthesis matrix, $\mathbf{\Omega}$ is defined by the values in
\file{inputfileGriddedData}{griddedData} and the
expression \config{value}, $\mathbf{A}$ is the spherical harmonic analysis matrix.
The resulting window matrix is written to a \file{matrix}{matrix} file.

The spherical harmonic degree range, and coefficient numbering are defined by
\config{minDegree}, \config{maxDegree}, and \configClass{numbering}{sphericalHarmonicsNumberingType}.

Note that a proper window function $\mathbf{\Omega}$ should contain values in the range [0, 1].
The window function $\mathbf{\Omega}$ can feature a smooth transition between 0 and 1 to avoid ringing effects.
)";

/***********************************************/

#include "programs/program.h"
#include "base/sphericalHarmonics.h"
#include "parser/dataVariables.h"
#include "files/fileMatrix.h"
#include "files/fileGriddedData.h"
#include "classes/grid/grid.h"
#include "classes/kernel/kernel.h"
#include "classes/sphericalHarmonicsNumbering/sphericalHarmonicsNumbering.h"

/***** CLASS ***********************************/

/** @brief Create a spherical harmonic window matrix.
* @ingroup programsGroup */
class FilterMatrixWindowedPotentialCoefficients
{
public:
  void run(Config &config);
};

GROOPS_REGISTER_PROGRAM(FilterMatrixWindowedPotentialCoefficients, PARALLEL, "create a spherical harmonic window matrix.", Misc, PotentialCoefficients, Matrix)

/***********************************************/

void FilterMatrixWindowedPotentialCoefficients::run(Config &config)
{
  try
  {
    FileName  fileNameOut, fileNameIn;
    KernelPtr kernel;
    UInt      minDegree, maxDegree;
    Double    GM, R;
    SphericalHarmonicsNumberingPtr numbering;
    ExpressionVariablePtr exprValue;

    readConfig(config, "outputfileWindowMatrix", fileNameOut, Config::OPTIONAL, "",  "");
    readConfig(config, "inputfileGriddedData",   fileNameIn,  Config::MUSTSET,  "",  "gridded data which defines the window function in space domain");
    readConfig(config, "value",                  exprValue,   Config::MUSTSET,  "data0", "expression to compute the window function (input columns are named data0, data1, ...)");
    readConfig(config, "kernel",                 kernel,      Config::MUSTSET,  "",  "kernel for windowing");
    readConfig(config, "minDegree",              minDegree,   Config::DEFAULT,  "0", "");
    readConfig(config, "maxDegree",              maxDegree,   Config::MUSTSET,  "",  "");
    readConfig(config, "GM",                     GM,          Config::DEFAULT,  STRING_DEFAULT_GM, "Geocentric gravitational constant");
    readConfig(config, "R",                      R,           Config::DEFAULT,  STRING_DEFAULT_R,  "reference radius");
    readConfig(config, "numbering",              numbering,   Config::MUSTSET,  "",  "numbering scheme for solution vector");
    if(isCreateSchema(config)) return;

    std::vector<UInt> n, m, cs;
    numbering->numbering(maxDegree, minDegree, n, m, cs);
    const UInt coefficientsCount = numbering->parameterCount(maxDegree, minDegree);

    GriddedData grid;
    readFileGriddedData(fileNameIn, grid);
    const std::vector<Vector3d> points = grid.points;
    const std::vector<Double>   areas  = grid.areas;

    // evaluate expression
    // -------------------
    auto varList = config.getVarList();
    std::set<std::string> usedVariables;
    exprValue->usedVariables(varList, usedVariables);
    addDataVariables(grid, varList, usedVariables);
    exprValue->simplify(varList);

    Vector windowFunction(points.size());
    for(UInt i=0; i<points.size(); i++)
    {
      evaluateDataVariables(grid, i, varList);
      windowFunction(i) = exprValue->evaluate(varList);
    }

    logStatus<<"create window matrix"<<Log::endl;
    Matrix A(coefficientsCount, coefficientsCount);
    const UInt blockSize = 256; // compute block of points to speed up computation
    Parallel::forEach((points.size()+blockSize-1)/blockSize, [&](UInt idx)
    {
      const UInt pointCount = std::min(points.size(), (idx+1)*blockSize) - idx*blockSize;
      Matrix A1(pointCount, coefficientsCount);
      Matrix A2(pointCount, coefficientsCount);

      for(UInt i=0; i<pointCount; i++)
      {
        const UInt   idPoint = idx*blockSize + i;
        const Double r       = points.at(idPoint).r();

        Matrix Cnm, Snm;
        SphericalHarmonics::CnmSnm(normalize(points.at(idPoint)), maxDegree, Cnm, Snm);

        Vector k1 = std::sqrt(areas.at(idPoint)/(4*PI)) * kernel->inverseCoefficients(points.at(idPoint), maxDegree);
        Vector k2 = std::sqrt(areas.at(idPoint)/(4*PI)) * kernel->coefficients(points.at(idPoint), maxDegree);
        for(UInt n=0; n<=maxDegree; n++)
        {
          k1(n) *= std::pow(R/r, n+1);
          k2(n) *= std::pow(r/R, n+1);
        }

        for(UInt k=0; k<A.columns(); k++)
        {
          A1(i, k) = k1(n[k]) * (cs[k] ? Cnm(n[k], m[k]) : Snm(n[k], m[k])) * windowFunction(idPoint);
          A2(i, k) = k2(n[k]) * (cs[k] ? Cnm(n[k], m[k]) : Snm(n[k], m[k]));
        }
      }

      matMult(1., A2.trans(), A1, A);
    });
    Parallel::reduceSum(A);

    // write
    // -----
    if(Parallel::isMaster())
    {
      logStatus<<"writing window matrix to file <"<<fileNameOut<<">"<<Log::endl;
      writeFileMatrix(fileNameOut, A);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
