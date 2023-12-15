/***********************************************/
/**
* @file griddedTopography2PotentialCoefficients.cpp
*
* @brief Estimate potential coefficients from digital terrain models.
*
* @author Torsten Mayer-Guerr
* @date 2011-11-01
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
Estimate potential coefficients from digital terrain models.
Coefficients for interior $(1/r)^{n+1}$ and exterior ($r^n$) are computed.
)";

/***********************************************/

#include "programs/program.h"
#include "base/legendreFunction.h"
#include "parser/dataVariables.h"
#include "files/fileGriddedData.h"
#include "files/fileSphericalHarmonics.h"
#include "misc/miscGriddedData.h"

/***** CLASS ***********************************/

/** @brief Estimate potential coefficients from digital terrain models.
* @ingroup programsGroup */
class GriddedTopography2PotentialCoefficients
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(GriddedTopography2PotentialCoefficients, PARALLEL, "Estimate potential coefficients from digital terrain models", Grid, PotentialCoefficients)

/***********************************************/

void GriddedTopography2PotentialCoefficients::run(Config &config, Parallel::CommunicatorPtr comm)
{
  try
  {
    FileName              fileNameOutExterior, fileNameOutInterior;
    FileName              fileNameInGrid;
    ExpressionVariablePtr expressionUpper, expressionLower, expressionRho;
    Double                factor;
    UInt                  minDegree, maxDegree;
    Double                GM, R;

    readConfig(config, "outputfilePotentialCoefficients",         fileNameOutExterior, Config::OPTIONAL, "", "");
    readConfig(config, "outputfilePotentialCoefficientsInterior", fileNameOutInterior, Config::OPTIONAL, "", "");
    readConfig(config, "inputfileGriddedData", fileNameInGrid,  Config::MUSTSET, "",      "Digital Terrain Model");
    readConfig(config, "density",              expressionRho,   Config::DEFAULT, "2670",  "expression [kg/m^3]");
    readConfig(config, "radialUpperBound",     expressionUpper, Config::DEFAULT, "data0", "expression (variables 'L', 'B', 'height', 'data', and 'area' are taken from the gridded data");
    readConfig(config, "radialLowerBound",     expressionLower, Config::DEFAULT, "0",     "expression (variables 'L', 'B', 'height', 'data', and 'area' are taken from the gridded data");
    readConfig(config, "factor",               factor,          Config::DEFAULT, "1.0",   "the result is multiplied by this factor");
    readConfig(config, "minDegree",            minDegree,       Config::DEFAULT, "0",     "");
    readConfig(config, "maxDegree",            maxDegree,       Config::MUSTSET, "",      "");
    readConfig(config, "GM",                   GM,              Config::DEFAULT, STRING_DEFAULT_GM, "Geocentric gravitational constant");
    readConfig(config, "R",                    R,               Config::DEFAULT, STRING_DEFAULT_R,  "reference radius");
    if(isCreateSchema(config)) return;

    std::vector<Double> lambda, phi;
    Matrix              rLower, rUpper, rho;
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
      logStatus<<"evaluate upper and lower height"<<Log::endl;
      VariableList varList;
      addDataVariables(grid, varList);
      expressionUpper->simplify(varList);
      expressionLower->simplify(varList);
      expressionRho  ->simplify(varList);

      std::vector<Double> radius(grid.heights.size());
      for(UInt i=0; i<radius.size(); i++)
        radius.at(i) = grid.ellipsoid(Angle(0), grid.latitudes.at(i), grid.heights.at(i)).r();

      rLower = rUpper = rho = Matrix(grid.latitudes.size(), grid.longitudes.size());
      Single::forEach(grid.latitudes.size(), [&](UInt i)
      {
        for(UInt k=0; k<grid.longitudes.size(); k++)
        {
          evaluateDataVariables(grid, i, k, varList);
          rUpper(i,k) = radius.at(i) + expressionUpper->evaluate(varList);
          rLower(i,k) = radius.at(i) + expressionLower->evaluate(varList);
          rho(i,k)    = expressionRho->evaluate(varList);
        }
      });

      grid.cellBorders(lambda, phi);
      for(UInt i=0; i<phi.size(); i++)
        phi.at(i) = grid.ellipsoid(Angle(0), Angle(phi.at(i)), 0.).phi(); // geocentric
    } // if(Parallel::isMaster(comm))

    if(Parallel::size(comm) > 1)
      logStatus<<"broadcast data"<<Log::endl;
    Parallel::broadCast(rUpper,  0, comm);
    Parallel::broadCast(rLower,  0, comm);
    Parallel::broadCast(rho,     0, comm);
    Parallel::broadCast(lambda,  0, comm);
    Parallel::broadCast(phi,     0, comm);

    // precompute integral_sin, integral_cos
    // -------------------------------------
    Matrix cosm(lambda.size()-1, maxDegree+1);
    Matrix sinm(lambda.size()-1, maxDegree+1);
    for(UInt i=0; i<lambda.size()-1; i++)
    {
      cosm(i,0) = std::remainder(lambda.at(i+1)-lambda.at(i), 2*PI);
      for(UInt m=1; m<=maxDegree; m++)
      {
        cosm(i,m) =  (std::sin(m*lambda.at(i+1)) - std::sin(m*lambda.at(i)))/m;
        sinm(i,m) = -(std::cos(m*lambda.at(i+1)) - std::cos(m*lambda.at(i)))/m;
      }
    }

    // computing quadrature formular
    // -----------------------------
    logStatus<<"computing quadrature formular"<<Log::endl;
    Matrix cnmExt, snmExt, cnmInt, snmInt;
    if(!fileNameOutExterior.empty()) cnmExt = Matrix(maxDegree+1, Matrix::TRIANGULAR, Matrix::LOWER);
    if(!fileNameOutExterior.empty()) snmExt = Matrix(maxDegree+1, Matrix::TRIANGULAR, Matrix::LOWER);
    if(!fileNameOutInterior.empty()) cnmInt = Matrix(maxDegree+1, Matrix::TRIANGULAR, Matrix::LOWER);
    if(!fileNameOutInterior.empty()) snmInt = Matrix(maxDegree+1, Matrix::TRIANGULAR, Matrix::LOWER);
    Parallel::forEach(phi.size()-1, [&](UInt i)
    {
      const Matrix Pnm = LegendreFunction::integral(std::sin(std::min(phi.at(i), phi.at(i+1))),
                                                    std::sin(std::max(phi.at(i), phi.at(i+1))), maxDegree);
      // const Matrix Pnm = LegendreFunction::compute(std::sin(0.5*(phi.at(i+1)+phi.at(i))), maxDegree)
      //                  * std::fabs(std::sin(phi.at(i+1))-std::sin(phi.at(i))); // integral cos(phi) dPhi

      if(!fileNameOutExterior.empty())
      {
        Matrix fExt(lambda.size()-1, maxDegree+1);
        for(UInt k=0; k<lambda.size()-1; k++)
          if(rho(i,k) && (std::fabs(rUpper(i,k)-rLower(i,k)) > 0.001))
          {
            const Double term = factor * rho(i,k) * GRAVITATIONALCONSTANT/GM * R*R*R;
            const Double r1R  = rLower(i,k)/R;
            const Double r2R  = rUpper(i,k)/R;
            Double r1RnExt = term * std::pow(r1R, minDegree+3);
            Double r2RnExt = term * std::pow(r2R, minDegree+3);
            for(UInt n=minDegree; n<=maxDegree; n++)
            {
              fExt(k,n) = (r2RnExt-r1RnExt)/((2.*n+1)*(n+3.));
              r1RnExt *= r1R;
              r2RnExt *= r2R;
            } // for(n)
          } // for(i)

        for(UInt n=minDegree; n<=maxDegree; n++)
          for(UInt m=0; m<=n; m++)
          {
            cnmExt(n,m) += Pnm(n,m) * inner(cosm.column(m), fExt.column(n));
            snmExt(n,m) += Pnm(n,m) * inner(sinm.column(m), fExt.column(n));
          }
      } // if(!fileNameOutExterior.empty())

      if(!fileNameOutInterior.empty())
      {
        Matrix fInt(lambda.size()-1, maxDegree+1);
        for(UInt k=0; k<lambda.size()-1; k++)
         if(rho(i,k) && (std::fabs(rUpper(i,k)-rLower(i,k)) > 0.001))
         {
            const Double term = factor * rho(i,k) * GRAVITATIONALCONSTANT/GM * R*R*R;
            const Double Rr1  = R/rLower(i,k);
            const Double Rr2  = R/rUpper(i,k);
            Double r1RnInt = term * std::pow(Rr1, minDegree-2.);
            Double r2RnInt = term * std::pow(Rr2, minDegree-2.);
            for(UInt n=minDegree; n<=maxDegree; n++)
            {
              if(n != 2)
                fInt(k,n) = (r2RnInt-r1RnInt)/((2*n+1)*(2.-n));
              else
                fInt(k,n) = term * std::log(rUpper(i,k)/rLower(i,k))/(2*n+1);
              r1RnInt *= Rr1;
              r2RnInt *= Rr2;
            } // for(n)
          } // for(i)

        for(UInt n=minDegree; n<=maxDegree; n++)
          for(UInt m=0; m<=n; m++)
          {
            cnmInt(n,m) += Pnm(n,m) * inner(cosm.column(m), fInt.column(n));
            snmInt(n,m) += Pnm(n,m) * inner(sinm.column(m), fInt.column(n));
          }
      } // if(!fileNameOutInterior.empty())
    }, comm);
    if(!fileNameOutExterior.empty()) Parallel::reduceSum(cnmExt, 0, comm);
    if(!fileNameOutExterior.empty()) Parallel::reduceSum(snmExt, 0, comm);
    if(!fileNameOutInterior.empty()) Parallel::reduceSum(cnmInt, 0, comm);
    if(!fileNameOutInterior.empty()) Parallel::reduceSum(snmInt, 0, comm);

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
