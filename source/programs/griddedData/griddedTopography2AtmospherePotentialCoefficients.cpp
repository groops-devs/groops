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
#include "parser/dataVariables.h"
#include "files/fileGriddedData.h"
#include "files/fileSphericalHarmonics.h"
#include "misc/miscGriddedData.h"

/***** CLASS ***********************************/

/** @brief Estimate interior and exterior atmospheric potential coefficients above digital terrain models.
* @ingroup programsGroup */
class GriddedTopography2AtmospherePotentialCoefficients
{
  Double                rho;
  Double                ny;
  Double                factor;
  Double                GM, R;
  Double                upperBoundary;
  UInt                  minDegree, maxDegree;
  Matrix                cnm, snm, cnmInt, snmInt;
  Matrix                cosm, sinm;

  UInt                  rows, cols;
  std::vector<Double>   dLambda, dPhi;
  std::vector<Angle>    lambda;
  std::vector<Angle>    phi;
  std::vector<Double>   radius;
  Matrix                rLower, rUpper, topo;

  void computeCoefficientsRow(UInt row);

public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(GriddedTopography2AtmospherePotentialCoefficients, PARALLEL, "Estimate interoir and exterior atmospheric potential coefficients above digital terrain models", Grid, PotentialCoefficients)

/***********************************************/

void GriddedTopography2AtmospherePotentialCoefficients::run(Config &config, Parallel::CommunicatorPtr comm)
{
  try
  {
    FileName outName, outNameInt, gridName;
    ExpressionVariablePtr expressionLower;

    readConfig(config, "outputfilePotentialCoefficientsExterior", outName,    Config::MUSTSET,  "", "");
    readConfig(config, "outputfilePotentialCoefficientsInterior", outNameInt, Config::MUSTSET,  "", "");
    readConfig(config, "inputfileGriddedData",     gridName,        Config::MUSTSET,  "",       "Digital Terrain Model");
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

    if(Parallel::isMaster(comm))
    {
      // read rectangular grid
      // ---------------------
      logStatus<<"read grid from file <"<<gridName<<">"<<Log::endl;
      GriddedDataRectangular grid;
      readFileGriddedData(gridName, grid);
      MiscGriddedData::printStatistics(grid);

      grid.geocentric(lambda, phi, radius, dLambda, dPhi);
      rows = phi.size();
      cols = lambda.size();

      // evaluate expression for lower boundary
      // -------------------------------
      logStatus<<"evaluate expression for lower boundary"<<Log::endl;
      VariableList varList;
      addDataVariables(grid, varList);
      varList.undefineVariable("area");
      expressionLower->simplify(varList);

      topo = Matrix(rows,cols);
      for(UInt z=0; z<rows; z++)
        for(UInt s=0; s<cols; s++)
        {
          evaluateDataVariables(grid, z, s, varList);
          varList.setVariable("area",  dLambda.at(s)*dPhi.at(z)*cos(phi.at(z)) ); // area
          topo(z,s) = expressionLower->evaluate(varList);   //  Topography
        }
    } // if(Parallel::isMaster(comm))

    Parallel::broadCast(topo,    0, comm);
    Parallel::broadCast(radius,  0, comm);
    Parallel::broadCast(lambda,  0, comm);
    Parallel::broadCast(phi,     0, comm);
    Parallel::broadCast(dLambda, 0, comm);
    Parallel::broadCast(dPhi,    0, comm);
    rows = phi.size();
    cols = lambda.size();

    // precompute sin, cos
    // -------------------
    cosm = Matrix(lambda.size(), maxDegree+1);
    sinm = Matrix(lambda.size(), maxDegree+1);
    for(UInt i=0; i<lambda.size(); i++)
      for(UInt m=0; m<=maxDegree; m++)
      {
        cosm(i,m) = cos(m*static_cast<Double>(lambda.at(i)));
        sinm(i,m) = sin(m*static_cast<Double>(lambda.at(i)));
      }

    // computing quadrature formular
    // -----------------------------
    logStatus<<"computing quadrature formular"<<Log::endl;
    cnm = Matrix(maxDegree+1, Matrix::TRIANGULAR, Matrix::LOWER);
    snm = Matrix(maxDegree+1, Matrix::TRIANGULAR, Matrix::LOWER);
    cnmInt = Matrix(maxDegree+1, Matrix::TRIANGULAR, Matrix::LOWER);
    snmInt = Matrix(maxDegree+1, Matrix::TRIANGULAR, Matrix::LOWER);
    Parallel::forEach(phi.size(), [this](UInt i){computeCoefficientsRow(i);}, comm);
    Parallel::reduceSum(cnm,    0, comm);
    Parallel::reduceSum(snm,    0, comm);
    Parallel::reduceSum(cnmInt, 0, comm);
    Parallel::reduceSum(snmInt, 0, comm);

    // save potential coefficients
    // ---------------------------
    if(Parallel::isMaster(comm))
    {
      logStatus<<"write exterior potential coefficients to file <"<<outName<<">"<<Log::endl;
      writeFileSphericalHarmonics(outName, SphericalHarmonics(GM, R, cnm, snm));
      logStatus<<"write interior potential coefficients to file <"<<outNameInt<<">"<<Log::endl;
      writeFileSphericalHarmonics(outNameInt, SphericalHarmonics(GM, R, cnmInt, snmInt));
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GriddedTopography2AtmospherePotentialCoefficients::computeCoefficientsRow(UInt row)
{
  try
  {
    const Double dB    = this->dPhi.at(row);
    const Double cosB0 = cos(phi.at(row));
    const Double H0 = radius.at(row);
    Matrix f(lambda.size(), maxDegree+1);
    Matrix g(lambda.size(), maxDegree+1);

    for(UInt k=0; k<lambda.size(); k++)
    {
      const Double term = factor * rho * GRAVITATIONALCONSTANT/GM * cosB0*dLambda.at(k)*dB *R*R*R;
      const Double r1R  = (1.+upperBoundary/H0);
      const Double r2R  = (1.+topo(row,k)/H0);

      Double r1Rn = term *(pow(1.+upperBoundary/H0,3.+minDegree-ny));
      Double r2Rn = term *(pow(1.+topo(row,k)/H0,3.+minDegree-ny));
      Double r1RnInt = term *(pow(1.+upperBoundary/H0,2.-minDegree-ny));
      Double r2RnInt = term *(pow(1.+topo(row,k)/H0,2.-minDegree-ny));

      for(UInt n=minDegree; n<=maxDegree; n++)
      {
        if (n == ny-3.)
          f(k,n) = factor * rho * GRAVITATIONALCONSTANT/GM * cosB0*dLambda.at(k)*dB  * pow(1./R,n) * pow(H0,ny) * log((H0+upperBoundary)/(H0+topo(row,k)));
        else
          f(k,n) = (r1Rn-r2Rn)*pow(H0/R,n+3.)/((2.*n+1.)*(3.+n-ny));

        g(k,n) = (r1RnInt-r2RnInt)*pow(H0/R,2.-n)/((2.*n+1.)*(2.-n-ny));

        r1Rn *=r1R;
        r2Rn *=r2R;
        r1RnInt /=r1R;
        r2RnInt /=r2R;
      }
    }

    const Matrix Pnm = SphericalHarmonics::Pnm(Angle(PI/2-phi.at(row)), 1.0, maxDegree);
    for(UInt n=minDegree; n<=maxDegree; n++)
      for(UInt m=0; m<=n; m++)
      {
        cnm(n,m) += Pnm(n,m) * inner(cosm.column(m), f.column(n));
        snm(n,m) += Pnm(n,m) * inner(sinm.column(m), f.column(n));
        cnmInt(n,m) += Pnm(n,m) * inner(cosm.column(m), g.column(n));
        snmInt(n,m) += Pnm(n,m) * inner(sinm.column(m), g.column(n));
      }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
