/***********************************************/
/**
* @file griddedTopography2AtmospherePotentialCoefficients.cpp
*
* @brief Estimate interior and exterior potential coefficients for atmosphere above digital terrain models.
*
* @author Christian Pock
* @author Daniel Rieser
* @date 2013-10-16
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
Estimate interior and exterior potential coefficients for atmosphere above digital terrain models.
Coefficients for interior $(1/r)^{n+1}$ and exterior ($r^n$) are computed.
The density of the atmosphere is assumed to be (Sjöberg, 1998)
\begin{equation}
\rho_0\left(\frac{R}{R+h}\right)^\nu,
\end{equation}
where $R$ is the radial distance of the ellipsoid at each point, $h$ the radial height above the ellipsoid,
$\rho_0$ is \config{densitySeaLevel} and \config{nu} $\nu$ is a constant factor. The density is integrated
from \config{radialLowerBound} and \config{upperAtmosphericBoundary} above the ellipsoid.
The \config{radialLowerBound} is typically the topography and can be computed as expression at every point
from \configFile{inputfileGriddedData}{griddedData}.
)";

/***********************************************/

#include "programs/program.h"
#include "base/legendreFunction.h"
#include "parser/dataVariables.h"
#include "files/fileGriddedData.h"
#include "files/fileSphericalHarmonics.h"
#include "misc/miscGriddedData.h"

/***** CLASS ***********************************/

/** @brief Estimate interior and exterior atmospheric potential coefficients above digital terrain models.
* @ingroup programsGroup */
class GriddedTopography2AtmospherePotentialCoefficients
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(GriddedTopography2AtmospherePotentialCoefficients, PARALLEL, "Estimate interoir and exterior atmospheric potential coefficients above digital terrain models", Grid, PotentialCoefficients)

/***********************************************/

void GriddedTopography2AtmospherePotentialCoefficients::run(Config &config, Parallel::CommunicatorPtr comm)
{
  try
  {
    FileName              fileNameOutExterior, fileNameOutInterior, fileNameInGrid;
    ExpressionVariablePtr expressionLower;
    Double                rho, ny, upperBoundary, factor;
    UInt                  minDegree, maxDegree;
    Double                GM, R;

    readConfig(config, "outputfilePotentialCoefficientsExterior", fileNameOutExterior, Config::OPTIONAL, "", "");
    readConfig(config, "outputfilePotentialCoefficientsInterior", fileNameOutInterior, Config::OPTIONAL, "", "");
    readConfig(config, "inputfileGriddedData",     fileNameInGrid,  Config::MUSTSET,  "",       "Digital Terrain Model");
    readConfig(config, "densitySeaLevel",          rho,             Config::DEFAULT,  "1.225",  "[kg/m**3]");
    readConfig(config, "ny",                       ny,              Config::DEFAULT,  "680",    "Constant for Atmosphere");
    readConfig(config, "radialLowerBound",         expressionLower, Config::DEFAULT,  "data0",  "expression (variables 'L', 'B', 'height', 'data', and 'area' are taken from the gridded data");
    readConfig(config, "upperAtmosphericBoundary", upperBoundary,   Config::DEFAULT,  "11000",  "constant upper bound [m]");
    readConfig(config, "factor",                   factor,          Config::DEFAULT,  "1.0",    "the result is multiplied by this factor, set -1 to subtract the field");
    readConfig(config, "minDegree",                minDegree,       Config::DEFAULT,  "0",      "");
    readConfig(config, "maxDegree",                maxDegree,       Config::MUSTSET,  "",       "");
    readConfig(config, "GM",                       GM,              Config::DEFAULT,  STRING_DEFAULT_GM, "Geocentric gravitational constant");
    readConfig(config, "R",                        R,               Config::DEFAULT,  STRING_DEFAULT_R,  "reference radius");
    if(isCreateSchema(config)) return;

    // read rectangular grid
    // ---------------------
    std::vector<Double> lambda, phi, radius;
    Matrix              topo;
    if(Parallel::isMaster(comm))
    {
      // read rectangular grid
      // ---------------------
      logStatus<<"read grid from file <"<<fileNameInGrid<<">"<<Log::endl;
      GriddedDataRectangular grid;
      readFileGriddedData(fileNameInGrid, grid);
      MiscGriddedData::printStatistics(grid);

      // evaluate upper and lower height
      // -------------------------------
      logStatus<<"evaluate expression for lower boundary"<<Log::endl;
      VariableList varList;
      addDataVariables(grid, varList);
      expressionLower->simplify(varList);

      radius.resize(grid.heights.size());
      for(UInt i=0; i<radius.size(); i++)
        radius.at(i) = grid.ellipsoid(Angle(0), grid.latitudes.at(i), grid.heights.at(i)).r();

      topo = Matrix(grid.latitudes.size(), grid.longitudes.size());
      Single::forEach(grid.latitudes.size(), [&](UInt i)
      {
        for(UInt k=0; k<grid.longitudes.size(); k++)
        {
          evaluateDataVariables(grid, i, k, varList);
          topo(i, k) = expressionLower->evaluate(varList);   //  Topography
        }
      });

      grid.cellBorders(lambda, phi);
      for(UInt i=0; i<phi.size(); i++)
        phi.at(i) = grid.ellipsoid(Angle(0), Angle(phi.at(i)), 0.).phi(); // geocentric
    } // if(Parallel::isMaster(comm))

    Parallel::broadCast(topo,    0, comm);
    Parallel::broadCast(radius,  0, comm);
    Parallel::broadCast(lambda,  0, comm);
    Parallel::broadCast(phi,     0, comm);

    // precompute integral_sin, integral_cos
    // -------------------------------------
    Matrix cosm(lambda.size()-1, maxDegree+1);
    Matrix sinm(lambda.size()-1, maxDegree+1);
    for(UInt i=0; i<lambda.size()-1; i++)
    {
      cosm(i,0) = lambda.at(i+1)-lambda.at(i);
      for(UInt m=1; m<=maxDegree; m++)
      {
        cosm(i,m) =  (std::sin(m*lambda.at(i+1)) - std::sin(m*lambda.at(i)))/m;
        sinm(i,m) = -(std::cos(m*lambda.at(i+1)) - std::cos(m*lambda.at(i)))/m;
      }
    }

    // computing quadrature formular
    // -----------------------------
    logStatus<<"computing quadrature formular"<<Log::endl;
    Matrix cnmExt(maxDegree+1, Matrix::TRIANGULAR, Matrix::LOWER);
    Matrix snmExt(maxDegree+1, Matrix::TRIANGULAR, Matrix::LOWER);
    Matrix cnmInt(maxDegree+1, Matrix::TRIANGULAR, Matrix::LOWER);
    Matrix snmInt(maxDegree+1, Matrix::TRIANGULAR, Matrix::LOWER);
    Parallel::forEach(phi.size(), [&](UInt i)
    {
      //     const Matrix Pnm = LegendreFunction::integral(std::sin(phi.at(i+1)), std::sin(phi.at(i)), maxDegree);
      const Matrix Pnm = LegendreFunction::compute(std::cos(0.5*(PI-phi.at(i+1)-phi.at(i))), maxDegree)
                       * std::fabs(std::sin(phi.at(i+1))-std::sin(phi.at(i))); // integral cos(phi) dPhi

      const Double H0   = radius.at(i);
      const Double term = factor * rho * GRAVITATIONALCONSTANT/GM *R*R*R;

      Matrix fExt(lambda.size(), maxDegree+1);
      Matrix fInt(lambda.size(), maxDegree+1);
      for(UInt k=0; k<lambda.size(); k++)
      {
        const Double r1R  = (1.+upperBoundary/H0);
        const Double r2R  = (1.+topo(i,k)/H0);

        Double r1RnExt = term * std::pow(1.+upperBoundary/H0, 3.+minDegree-ny);
        Double r2RnExt = term * std::pow(1.+topo(i,k)/H0,     3.+minDegree-ny);
        Double r1RnInt = term * std::pow(1.+upperBoundary/H0, 2.-minDegree-ny);
        Double r2RnInt = term * std::pow(1.+topo(i,k)/H0,     2.-minDegree-ny);

        for(UInt n=minDegree; n<=maxDegree; n++)
        {
          if (n == ny-3.)
            fExt(k,n) = term/(R*R*R) * std::pow(1./R, n) * std::pow(H0, ny) * std::log((H0+upperBoundary)/(H0+topo(i,k)));
          else
            fExt(k,n) = (r1RnExt-r2RnExt) * std::pow(H0/R, n+3.)/((2.*n+1.)*(3.+n-ny));

          fInt(k,n) = (r1RnInt-r2RnInt) * std::pow(H0/R, 2.-n)/((2.*n+1.)*(2.-n-ny));

          r1RnExt *=r1R;
          r2RnExt *=r2R;
          r1RnInt /=r1R;
          r2RnInt /=r2R;
        }
      }

      for(UInt n=minDegree; n<=maxDegree; n++)
        for(UInt m=0; m<=n; m++)
        {
          cnmExt(n,m) += Pnm(n,m) * inner(cosm.column(m), fExt.column(n));
          snmExt(n,m) += Pnm(n,m) * inner(sinm.column(m), fExt.column(n));
          cnmInt(n,m) += Pnm(n,m) * inner(cosm.column(m), fInt.column(n));
          snmInt(n,m) += Pnm(n,m) * inner(sinm.column(m), fInt.column(n));
        }
    }, comm);
    Parallel::reduceSum(cnmExt, 0, comm);
    Parallel::reduceSum(snmExt, 0, comm);
    Parallel::reduceSum(cnmInt, 0, comm);
    Parallel::reduceSum(snmInt, 0, comm);

    // save potential coefficients
    // ---------------------------
    if(Parallel::isMaster(comm) && !fileNameOutExterior.empty())
    {
      logStatus<<"write potential coefficients to file <"<<fileNameOutExterior<<">"<<Log::endl;
      writeFileSphericalHarmonics(fileNameOutExterior, SphericalHarmonics(GM, R, cnmExt, snmExt));
    }

    if(Parallel::isMaster(comm) && !fileNameOutInterior.empty())
    {
      logStatus<<"write potential coefficients (interior) to file <"<<fileNameOutInterior<<">"<<Log::endl;
      writeFileSphericalHarmonics(fileNameOutInterior, SphericalHarmonics(GM, R, cnmInt, snmInt, TRUE));
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
