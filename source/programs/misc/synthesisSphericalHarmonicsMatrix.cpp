/***********************************************/
/**
* @file synthesisSphericalHarmonicsMatrix.cpp
*
* @brief Matrix for synthesizing spherical harmonics on a grid.
*
* @author Torsten Mayer-Guerr
* @date 2025-01-23
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program builds a linear operator matrix for spherical harmonic analysis or synthesis based on
the points defined in \configClass{grid}{gridType}. Depending on the chosen
\config{type} (synthesis, quadrature, or leastSquares), the resulting matrix can be used to:
\begin{itemize}
  \item \textbf{synthesis}: Map spherical harmonic coefficients to values on a grid,
  \item \textbf{quadrature}: Integrate grid-based functionals into spherical harmonic coefficients by
        a simple quadrature formula,
  \item \textbf{leastSquares}: Estimate coefficients from grid data via a least squares approach.
\end{itemize}

he spherical harmonic degree range is constrained by
\config{minDegree} and \config{maxDegree}, and the ordering of the coefficients is given by
\configClass{numbering}{sphericalHarmonicsNumberingType}. The reference gravitational
constant is \config{GM}, and the reference radius is \config{R}.

The computed matrix is written to \configFile{outputfileMatrix}{matrix} with dimensions
(number of grid points)~$\times$~(number of spherical harmonic coefficients). For
\config{type} = \emph{leastSquares}, the program applies a QR-based pseudo-inverse so that the
output matrix can directly form the normal-equation building blocks for a blockwise
least-squares solution in spherical harmonic space.

See also \program{Gravityfield2GriddedData}, \program{GriddedData2PotentialCoefficients},
\program{Gravityfield2SphericalHarmonicsVector}, and \program{MatrixCalculate} for additional
tools to convert between grids and spherical harmonics.
)";

/***********************************************/

#include "programs/program.h"
#include "base/sphericalHarmonics.h"
#include "files/fileMatrix.h"
#include "classes/grid/grid.h"
#include "classes/kernel/kernel.h"
#include "classes/sphericalHarmonicsNumbering/sphericalHarmonicsNumbering.h"

/***** CLASS ***********************************/

/** @brief Matrix for synthesizing spherical harmonics on a grid.
* @ingroup programsGroup */
class SynthesisSphericalHarmonicsMatrix
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(SynthesisSphericalHarmonicsMatrix, PARALLEL, "matrix for synthesizing spherical harmonics on a grid", Misc, PotentialCoefficients, Matrix)

/***********************************************/

void SynthesisSphericalHarmonicsMatrix::run(Config &config, Parallel::CommunicatorPtr comm)
{
  try
  {
    FileName    fileNameOut;
    GridPtr     grid;
    KernelPtr   kernel;
    UInt        minDegree, maxDegree;
    Double      GM, R;
    SphericalHarmonicsNumberingPtr numbering;
    enum Type {SYNTHESIS, QUADRATURE, LEASTSQUARES};
    Type        type;
    std::string choice;

    readConfig(config, "outputfileMatrix", fileNameOut, Config::OPTIONAL, "",  "");
    readConfig(config, "grid",             grid,        Config::MUSTSET,  "",  "");
    readConfig(config, "kernel",           kernel,      Config::MUSTSET,  "",  "");
    readConfig(config, "minDegree",        minDegree,   Config::DEFAULT,  "0", "");
    readConfig(config, "maxDegree",        maxDegree,   Config::MUSTSET,  "",  "");
    readConfig(config, "GM",               GM,          Config::DEFAULT,  STRING_DEFAULT_GM, "Geocentric gravitational constant");
    readConfig(config, "R",                R,           Config::DEFAULT,  STRING_DEFAULT_R,  "reference radius");
    readConfig(config, "numbering",        numbering,   Config::MUSTSET,  "",  "numbering scheme of sh coefficients");
    if(readConfigChoice(config, "type", choice, Config::MUSTSET, "", ""))
    {
      if(readConfigChoiceElement(config, "synthesis",    choice, "synthesize spherical harmonics on a grid")) type = SYNTHESIS;
      if(readConfigChoiceElement(config, "quadrature",   choice, "calculate spherical harmonics from grid"))  type = QUADRATURE;
      if(readConfigChoiceElement(config, "leastSquares", choice, "estimated spherical harmonics from grid"))  type = LEASTSQUARES;
      endChoice(config);
    }
    if(isCreateSchema(config)) return;

    std::vector<std::vector<UInt>> idxC, idxS;
    numbering->numbering(maxDegree, minDegree, idxC, idxS);

    logStatus<<"calcalute matrix for "<<numbering->parameterCount(maxDegree, minDegree)<<" coefficients"<<Log::endl;
    Matrix A = Matrix(grid->points().size(), numbering->parameterCount(maxDegree, minDegree));
    Parallel::forEach(grid->points().size(), [&](UInt k)
    {
      Matrix Cnm, Snm;
      Vector factors;
      SphericalHarmonics::CnmSnm(1./R * grid->points().at(k), maxDegree, Cnm, Snm, (type == QUADRATURE)/*isInterior*/);
      if(type == QUADRATURE)
        factors = grid->points().at(k).r()/(4*PI*GM) * grid->areas().at(k) * kernel->coefficients(grid->points().at(k), maxDegree);
      else
        factors = GM/R * kernel->inverseCoefficients(grid->points().at(k), maxDegree);

      for(UInt n=minDegree; n<=maxDegree; n++)
      {
        if(idxC[n][0]!=NULLINDEX) A(k, idxC[n][0]) = factors(n) * Cnm(n,0);
        for(UInt m=1; m<=n; m++)
        {
          if(idxC[n][m]!=NULLINDEX) A(k, idxC[n][m]) = factors(n) * Cnm(n,m);
          if(idxS[n][m]!=NULLINDEX) A(k, idxS[n][m]) = factors(n) * Snm(n,m);
        }
      }
    }, comm);
    Parallel::reduceSum(A, 0, comm);
    if(!Parallel::isMaster(comm))
      return;

    // A' := (A'PA)^(-1) A'P
    if(type == LEASTSQUARES)
    {
      logStatus<<"compute pseudo inverse"<<Log::endl;
      for(UInt k=0; k<grid->points().size(); k++)
        A.row(k) *= std::sqrt(grid->areas().at(k));
      const Vector tau = QR_decomposition(A);
      const Matrix R   = A.row(0, A.columns());
      generateQ(A, tau);
      triangularSolve(1., R, A.trans());
      for(UInt k=0; k<grid->points().size(); k++)
        A.row(k) *= std::sqrt(grid->areas().at(k));
    }

    // write
    // -----
    logStatus<<"write matrix to file <"<<fileNameOut<<">"<<Log::endl;
    writeFileMatrix(fileNameOut, (type == SYNTHESIS) ? A : A.trans());
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
