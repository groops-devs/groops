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
  FileName              fileNameOutExterior, fileNameOutInterior;
  Bool                  isExterior, isInterior;
  Double                factor;
  Double                GM, R;
  UInt                  minDegree, maxDegree;
  Matrix                cnmExt, snmExt;
  Matrix                cnmInt, snmInt;
  Matrix                cosm, sinm;

  UInt                  rows, cols;
  std::vector<Double>   dLambda, dPhi;
  std::vector<Angle>    lambda;
  std::vector<Angle>    phi;
  Matrix                rLower, rUpper, rho;

  void computeCoefficientsRow(UInt row);

public:
  void run(Config &config);
};

GROOPS_REGISTER_PROGRAM(GriddedTopography2PotentialCoefficients, PARALLEL, "Estimate potential coefficients from digital terrain models", Grid, PotentialCoefficients)

/***********************************************/

void GriddedTopography2PotentialCoefficients::run(Config &config)
{
  try
  {
    FileName              fileNameOutExterior, fileNameOutInterior;
    FileName              fileNameInGrid;
    ExpressionVariablePtr expressionUpper, expressionLower, expressionRho;

    isExterior = readConfig(config, "outputfilePotentialCoefficients",         fileNameOutExterior, Config::OPTIONAL, "", "");
    isInterior = readConfig(config, "outputfilePotentialCoefficientsInterior", fileNameOutInterior, Config::OPTIONAL, "", "");
    readConfig(config, "inputfileGriddedData",     fileNameInGrid,  Config::MUSTSET,  "",      "Digital Terrain Model");
    readConfig(config, "density",                  expressionRho,   Config::DEFAULT,  "2670",  "expression [kg/m**3]");
    readConfig(config, "radialUpperBound",         expressionUpper, Config::DEFAULT,  "data0", "expression (variables 'L', 'B', 'height', 'data', and 'area' are taken from the gridded data");
    readConfig(config, "radialLowerBound",         expressionLower, Config::DEFAULT,  "0",     "expression (variables 'L', 'B', 'height', 'data', and 'area' are taken from the gridded data");
    readConfig(config, "factor",                   factor,          Config::DEFAULT,  "1.0",   "the result is multplied by this factor");
    readConfig(config, "minDegree",                minDegree,       Config::DEFAULT,  "0",     "");
    readConfig(config, "maxDegree",                maxDegree,       Config::MUSTSET,  "",      "");
    readConfig(config, "GM",                       GM,              Config::DEFAULT,  STRING_DEFAULT_GM, "Geocentric gravitational constant");
    readConfig(config, "R",                        R,               Config::DEFAULT,  STRING_DEFAULT_R,  "reference radius");
    if(isCreateSchema(config)) return;

    if(Parallel::isMaster())
    {
      // read rectangular grid
      // ---------------------
      logStatus<<"read grid from file <"<<fileNameInGrid<<">"<<Log::endl;
      GriddedDataRectangular grid;
      readFileGriddedData(fileNameInGrid, grid);
      MiscGriddedData::printStatistics(grid);

      std::vector<Double> radius;
      grid.geocentric(lambda, phi, radius, dLambda, dPhi);
      rows = phi.size();
      cols = lambda.size();

      // evaluate upper and lower height
      // -------------------------------
      logStatus<<"evaluate upper and lower height"<<Log::endl;
      auto varList = config.getVarList();
      std::set<std::string> usedVariables;
      expressionUpper->usedVariables(varList, usedVariables);
      expressionLower->usedVariables(varList, usedVariables);
      expressionRho  ->usedVariables(varList, usedVariables);
      addDataVariables(grid, varList, usedVariables);
      addVariable("area", varList);
      expressionUpper->simplify(varList);
      expressionLower->simplify(varList);
      expressionRho  ->simplify(varList);

      rLower = rUpper =  rho = Matrix(rows,cols);
      for(UInt z=0; z<rows; z++)
        for(UInt s=0; s<cols; s++)
        {
          evaluateDataVariables(grid, z, s, varList);
          varList["area"]->setValue( dLambda.at(s)*dPhi.at(z)*cos(phi.at(z)) ); // area
          rUpper(z,s) = radius.at(z) + expressionUpper->evaluate(varList);
          rLower(z,s) = radius.at(z) + expressionLower->evaluate(varList);
          rho(z,s)    = expressionRho->evaluate(varList);
        }
    } // if(Parallel::isMaster())

    Parallel::broadCast(rUpper);
    Parallel::broadCast(rLower);
    Parallel::broadCast(rho);
    Parallel::broadCast(lambda);
    Parallel::broadCast(phi);
    Parallel::broadCast(dLambda);
    Parallel::broadCast(dPhi);
    rows = phi.size();
    cols = lambda.size();

    // precompute integral_sin, integral_cos
    // -------------------------------------
    cosm = Matrix(lambda.size(), maxDegree+1);
    sinm = Matrix(lambda.size(), maxDegree+1);
    for(UInt i=0; i<lambda.size(); i++)
    {
      cosm(i,0) = dLambda.at(i);
      for(UInt m=1; m<=maxDegree; m++)
      {
        cosm(i,m) = cos(m*static_cast<Double>(lambda.at(i))) * 2*sin(m*dLambda.at(i)/2)/m;
        sinm(i,m) = sin(m*static_cast<Double>(lambda.at(i))) * 2*sin(m*dLambda.at(i)/2)/m;
      }
    }

    // computing quadrature formular
    // -----------------------------
    logStatus<<"computing quadrature formular"<<Log::endl;
    if(isExterior) cnmExt = Matrix(maxDegree+1, Matrix::TRIANGULAR, Matrix::LOWER);
    if(isExterior) snmExt = Matrix(maxDegree+1, Matrix::TRIANGULAR, Matrix::LOWER);
    if(isInterior) cnmInt = Matrix(maxDegree+1, Matrix::TRIANGULAR, Matrix::LOWER);
    if(isInterior) snmInt = Matrix(maxDegree+1, Matrix::TRIANGULAR, Matrix::LOWER);
    Parallel::forEach(phi.size(), [this](UInt i){computeCoefficientsRow(i);});
    if(isExterior) Parallel::reduceSum(cnmExt);
    if(isExterior) Parallel::reduceSum(snmExt);
    if(isInterior) Parallel::reduceSum(cnmInt);
    if(isInterior) Parallel::reduceSum(snmInt);

    // save potential coefficients
    // ---------------------------
    if(Parallel::isMaster())
    {
      if(!fileNameOutExterior.empty())
      {
        logStatus<<"write potential coefficients to file <"<<fileNameOutExterior<<">"<<Log::endl;
        writeFileSphericalHarmonics(fileNameOutExterior, SphericalHarmonics(GM, R, cnmExt, snmExt));
      }

      if(!fileNameOutInterior.empty())
      {
        logStatus<<"write potential coefficients (interior) to file <"<<fileNameOutInterior<<">"<<Log::endl;
        writeFileSphericalHarmonics(fileNameOutInterior, SphericalHarmonics(GM, R, cnmInt, snmInt, TRUE));
      }
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GriddedTopography2PotentialCoefficients::computeCoefficientsRow(UInt row)
{
  try
  {
    Matrix fExt, fInt;

    if(isExterior)
    {
      fExt = Matrix(lambda.size(), maxDegree+1);
      for(UInt k=0; k<lambda.size(); k++)
      {
        if(fabs(rUpper(row,k)-rLower(row,k))<0.001)
          continue;

        const Double dense= rho(row,k);
        const Double term = factor * dense * GRAVITATIONALCONSTANT/GM * R*R*R;
        const Double r1R  = rLower(row,k)/R;
        const Double r2R  = rUpper(row,k)/R;
        Double r1RnExt = term * pow(r1R, minDegree+3);
        Double r2RnExt = term * pow(r2R, minDegree+3);

        for(UInt n=minDegree; n<=maxDegree; n++)
        {
          fExt(k,n) = (r2RnExt-r1RnExt)/((2.*n+1)*(n+3.));
          r1RnExt *= r1R;
          r2RnExt *= r2R;
        } // for(n)
      } // for(i)
    } // if(isExterior)

    if(isInterior)
    {
      fInt = Matrix(lambda.size(), maxDegree+1);
      for(UInt k=0; k<lambda.size(); k++)
      {
        if(fabs(rUpper(row,k)-rLower(row,k))<0.001)
          continue;

        const Double dense= rho(row,k);
        const Double term = factor * dense * GRAVITATIONALCONSTANT/GM * R*R*R;
        const Double Rr1  = R/rLower(row,k);
        const Double Rr2  = R/rUpper(row,k);
        Double r1RnInt = term * pow(Rr1, minDegree-2.);
        Double r2RnInt = term * pow(Rr2, minDegree-2.);

        for(UInt n=minDegree; n<=maxDegree; n++)
        {
          if(n!=2.)
            fInt(k,n) = (r2RnInt-r1RnInt)/((2.*n+1)*(2.-n));
          else
            fInt(k,n) = term * log(rUpper(row,k)/rLower(row,k))/(2.*n+1);
          r1RnInt *= Rr1;
          r2RnInt *= Rr2;
        } // for(n)
      } // for(i)
    } // if(isInterior)

    // --------------------------------

//     const Matrix Pnm = LegendreFunction::integral(cos(PI/2-phi.at(row)+fabs(dPhi.at(row))/2),
//                                                   cos(PI/2-phi.at(row)-fabs(dPhi.at(row))/2), maxDegree);
    const Matrix Pnm = LegendreFunction::compute(cos(PI/2-phi.at(row)), maxDegree)
                     * (cos(PI/2-phi.at(row)-fabs(dPhi.at(row))/2) - cos(PI/2-phi.at(row)+fabs(dPhi.at(row))/2));

    if(isExterior)
    {
      Matrix cnm(maxDegree+1, Matrix::TRIANGULAR, Matrix::LOWER);
      Matrix snm(maxDegree+1, Matrix::TRIANGULAR, Matrix::LOWER);
      for(UInt n=minDegree; n<=maxDegree; n++)
        for(UInt m=0; m<=n; m++)
        {
          cnm(n,m) += Pnm(n,m) * inner(cosm.column(m), fExt.column(n));
          snm(n,m) += Pnm(n,m) * inner(sinm.column(m), fExt.column(n));
        }
      cnmExt += cnm;
      snmExt += snm;
    } // if(isExterior)

    if(isInterior)
    {
      Matrix cnm(maxDegree+1, Matrix::TRIANGULAR, Matrix::LOWER);
      Matrix snm(maxDegree+1, Matrix::TRIANGULAR, Matrix::LOWER);
      for(UInt n=minDegree; n<=maxDegree; n++)
        for(UInt m=0; m<=n; m++)
        {
          cnm(n,m) += Pnm(n,m) * inner(cosm.column(m), fInt.column(n));
          snm(n,m) += Pnm(n,m) * inner(sinm.column(m), fInt.column(n));
        }
      cnmInt += cnm;
      snmInt += snm;
    } // if(isInterior)
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
