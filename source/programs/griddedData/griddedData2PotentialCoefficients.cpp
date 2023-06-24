/***********************************************/
/**
* @file griddedData2PotentialCoefficients.cpp
*
* @brief Estimate potential coefficients from gridded gravity field functionals.
* With quadrature formular or
* least squares adjustment with block diagonal normals (order by order).
*
* @author Annette Eicker
* @author Torsten Mayer-Guerr
* @date 2005-05-01
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program estimate potential coefficients from \configFile{inputfileGriddedData}{griddedData}
gravity field functionals. It used a simple quadrature formular
\begin{equation}
  c_{nm} = \frac{1}{4\pi}\frac{R}{GM} \sum_i f_i \left(\frac{r_i}{R}\right)^{n+1} k_n C_{nm}(\lambda_i,\vartheta_i)\,\Delta\Phi_i
\end{equation}
or a \config{leastSquares} adjustment with block diagonal normal matrix (order by order).
For the latter one the data must be regular distributed.

The \config{value}s $f_i$ and the \config{weight}s $\Delta\Phi_i$ are expressions
using the common data variables for grids, see \reference{dataVariables}{general.parser:dataVariables}.
The type of the gridded data (e.g gravity anomalies or geoid heights)
must be set with \configClass{kernel}{kernelType} $k_n$.

The expansion is limited in the range between \config{minDegree}
and \config{maxDegree} inclusively. The coefficients are related
to the reference radius~\config{R} and the Earth gravitational constant \config{GM}.

For irregular distributed data and using the full variance covariance matrix use
\program{NormalsSolverVCE} together with \configClass{oberservation:terrestrial}{observationType:terrestrial}
and \configClass{parametrizationGravity:sphericalHarmonics}{parametrizationGravityType:sphericalHarmonics}.
)";

/***********************************************/

#include "programs/program.h"
#include "parser/dataVariables.h"
#include "files/fileSphericalHarmonics.h"
#include "files/fileGriddedData.h"
#include "classes/kernel/kernel.h"
#include "misc/miscGriddedData.h"

/***** CLASS ***********************************/

/** @brief Estimate potential coefficients from gridded gravity field functionals.
* With quadrature formular or
* least squares adjustment with block diagonal normals (order by order).
* @ingroup programsGroup */
class GriddedData2PotentialCoefficients
{
  KernelPtr             kernel;
  Double                GM, R;
  UInt                  minDegree, maxDegree;

  GriddedData           grid;
  std::vector<Angle>    lambda;
  std::vector<Angle>    phi;
  std::vector<Double>   radius;

  UInt                  blocking;
  Matrix                cossinm;
  std::vector<Matrix>   N;
  std::vector<Vector>   n;
  Double                lPl;


  SphericalHarmonics computeQuadrature(Bool isRectangle, Parallel::CommunicatorPtr comm);
  SphericalHarmonics computeLeastSquares(Bool isRectangle, Parallel::CommunicatorPtr comm);
  void               buildNormals(UInt i);
  void               buildNormalsFast(UInt i);

public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(GriddedData2PotentialCoefficients, PARALLEL, "Estimate potential coefficients from gridded gravity field functionals", Grid, PotentialCoefficients)

/***********************************************/

void GriddedData2PotentialCoefficients::run(Config &config, Parallel::CommunicatorPtr comm)
{
  try
  {
    FileName      fileNameOut, fileNameGrid;
    ExpressionVariablePtr exprValue, exprArea;
    Bool          useLeastSquares;

    readConfig(config, "outputfilePotentialCoefficients", fileNameOut,     Config::MUSTSET,  "",      "");
    readConfig(config, "inputfileGriddedData",            fileNameGrid,    Config::MUSTSET,  "",      "");
    readConfig(config, "value",                           exprValue,       Config::MUSTSET,  "data0", "expression to compute values (input columns are named data0, data1, ...)");
    readConfig(config, "weight",                          exprArea,        Config::MUSTSET,  "area",  "expression to compute values (input columns are named data0, data1, ...)");
    readConfig(config, "kernel",                          kernel,          Config::MUSTSET,  "",      "data type of input values");
    readConfig(config, "minDegree",                       minDegree,       Config::DEFAULT,  "0",     "");
    readConfig(config, "maxDegree",                       maxDegree,       Config::MUSTSET,  "",      "");
    readConfig(config, "GM",                              GM,              Config::DEFAULT,  STRING_DEFAULT_GM, "Geocentric gravitational constant");
    readConfig(config, "R",                               R,               Config::DEFAULT,  STRING_DEFAULT_R,  "reference radius");
    readConfig(config, "leastSquares",                    useLeastSquares, Config::DEFAULT,  "0",     "false: quadrature formular, true: least squares adjustment order by order");
    if(isCreateSchema(config)) return;

    // reading grids
    // -------------
    logStatus<<"read grid from file <"<<fileNameGrid<<">"<<Log::endl;
    readFileGriddedData(fileNameGrid, grid);
    if(!grid.areas.size())
      logWarningOnce<<"  no areas"<<Log::endl;
    if(!grid.values.size())
      logWarningOnce<<"  no values"<<Log::endl;

    // evaluate expression
    // -------------------
    VariableList varList;
    addDataVariables(grid, varList);
    exprValue->simplify(varList);
    exprArea ->simplify(varList);

    for(UInt i=0; i<grid.points.size(); i++)
    {
      evaluateDataVariables(grid, i, varList);
      grid.values.at(0).at(i) = exprValue->evaluate(varList);
      grid.areas.at(i)        = exprArea->evaluate(varList);
    }
    grid.values.resize(1);

    MiscGriddedData::printStatistics(grid);

    // Are fast algorithms possible?
    // -----------------------------
    Bool isRectangle = grid.isRectangle(lambda, phi, radius);
    if(isRectangle)
      logInfo<<"  regular grid -> use of fast algorithm possible"<<Log::endl;

    SphericalHarmonics harm;
    if(useLeastSquares)
      harm = computeLeastSquares(isRectangle, comm);
    else
      harm = computeQuadrature(isRectangle, comm);

    // write potential coefficients
    // ----------------------------
    if(Parallel::isMaster(comm))
    {
      logStatus<<"write potential coefficients to file <"<<fileNameOut<<">"<<Log::endl;
      writeFileSphericalHarmonics(fileNameOut, harm);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

SphericalHarmonics GriddedData2PotentialCoefficients::computeQuadrature(Bool /*isRectangle*/, Parallel::CommunicatorPtr comm)
{
  try
  {
    Matrix cnm(maxDegree+1, Matrix::TRIANGULAR, Matrix::LOWER);
    Matrix snm(maxDegree+1, Matrix::TRIANGULAR, Matrix::LOWER);
    logStatus<<"computing quadrature formular"<<Log::endl;
    Parallel::forEach(grid.points.size(), [&](UInt i)
    {
      const Vector kn = kernel->coefficients(grid.points.at(i), maxDegree);
      Matrix Cnm, Snm;
      SphericalHarmonics::CnmSnm(grid.points.at(i)/grid.points.at(i).r(), maxDegree, Cnm, Snm);
      for(UInt n=minDegree; n<=maxDegree; n++)
      {
        axpy((kn(n)* R/(4*PI*GM) * std::pow(grid.points.at(i).r()/R, n+1) * grid.values.at(0).at(i) * grid.areas.at(i)), Cnm.row(n), cnm.row(n));
        axpy((kn(n)* R/(4*PI*GM) * std::pow(grid.points.at(i).r()/R, n+1) * grid.values.at(0).at(i) * grid.areas.at(i)), Snm.row(n), snm.row(n));
      }
    }, comm);
    Parallel::reduceSum(cnm, 0, comm);
    Parallel::reduceSum(snm, 0, comm);
    return SphericalHarmonics(GM, R, cnm, snm);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

SphericalHarmonics GriddedData2PotentialCoefficients::computeLeastSquares(Bool isRectangle, Parallel::CommunicatorPtr comm)
{
  try
  {
    logStatus<<"least squares adjustment (order by order)"<<Log::endl;

    lPl = 0;
    N.resize(2*maxDegree+1);
    n.resize(2*maxDegree+1);
    // order m=0
    N.at(0) = Matrix(maxDegree+1, Matrix::SYMMETRIC, Matrix::UPPER);
    n.at(0) = Vector(maxDegree+1);
    // remaining orders
    for(UInt m=1; m<=maxDegree; m++)
    {
      N.at(2*m-1) = Matrix(maxDegree+1-m, Matrix::SYMMETRIC, Matrix::UPPER);
      N.at(2*m+0) = Matrix(maxDegree+1-m, Matrix::SYMMETRIC, Matrix::UPPER);
      n.at(2*m-1) = Vector(maxDegree+1-m);
      n.at(2*m+0) = Vector(maxDegree+1-m);
    }

    // accumulate normals
    // ------------------
    if(!isRectangle)
    {
      blocking = 100;
      Parallel::forEach(grid.points.size()/blocking+1, [this](UInt i){buildNormals(i);}, comm);
      Parallel::reduceSum(lPl, 0, comm);
    }
    else
    {
      // fast algorithm
      lPl = 0;
      for(UInt i=0; i<grid.points.size(); i++)
        lPl += grid.values.at(0).at(i) * grid.areas.at(i)/(4*PI) * grid.values.at(0).at(i);

      cossinm = Matrix(lambda.size(), 2*maxDegree+1);
      for(UInt k=0; k<lambda.size(); k++)
      {
        cossinm(k,0) = 1.;
        for(UInt m=1; m<=maxDegree; m++)
        {
          cossinm(k,2*m-1) = cos(m*static_cast<Double>(lambda.at(k)));
          cossinm(k,2*m+0) = sin(m*static_cast<Double>(lambda.at(k)));
        }
      }

      Parallel::forEach(phi.size(), [this](UInt i){buildNormalsFast(i);}, comm);
    } // if(isRectangle)

    for(UInt i=0; i<N.size(); i++)
    {
      Parallel::reduceSum(N.at(i), 0, comm);
      Parallel::reduceSum(n.at(i), 0, comm);
    }

    // solve normals
    // -------------
    if(Parallel::isMaster(comm))
    {
      logStatus<<"solve the system of equations"<<Log::endl;
      std::vector<Vector> x(2*maxDegree+1);
      std::vector<Vector> sigma2x(2*maxDegree+1);
      Double ePe = lPl;
      UInt   parameterCount = 0;
      for(UInt i=0; i<x.size(); i++)
      {
        parameterCount += n.at(i).rows();
        for(UInt k=0; k<N.at(i).rows(); k++)
          if(N.at(i)(k,k)==0)
          {
            N.at(i)(k,k) += 1;
            parameterCount--;
            logWarning<<k<<". parameter has zero diagonal element -> set to one"<<Log::endl;
          }
        x.at(i) = solve(N.at(i), n.at(i));
        ePe -= inner(x.at(i), n.at(i)); // quad sum of residuals
        inverse(N.at(i)); // inverse of the cholesky matrix
        sigma2x.at(i) = Vector(x.at(i).rows());
        for(UInt k=0; k<N.at(i).rows(); k++)
          sigma2x.at(i)(k) = quadsum(N.at(i).slice(k,k,1,N.at(i).columns()-k));
      }

      // aposteriori sigma
      Double sigma2 = ePe/(grid.points.size()-parameterCount);
      if((sigma2<=0)||(sigma2!=sigma2))
      {
        logWarning<<"sigma^2 = "<<sigma2<<" not applied to covariance matrix"<<Log::endl;
        sigma2 = 1.;
      }
      logInfo<<"  aposteriori sigma = "<<sqrt(sigma2)<<Log::endl;
      for(UInt i=0; i<x.size(); i++)
        sigma2x.at(i) *= sigma2;

      // potential coefficients
      Matrix cnm      (maxDegree+1, Matrix::TRIANGULAR, Matrix::LOWER);
      Matrix snm      (maxDegree+1, Matrix::TRIANGULAR, Matrix::LOWER);
      Matrix sigma2cnm(maxDegree+1, Matrix::TRIANGULAR, Matrix::LOWER);
      Matrix sigma2snm(maxDegree+1, Matrix::TRIANGULAR, Matrix::LOWER);

      // order m=0
      copy(x.at(0),       cnm.slice(0,0,maxDegree+1,1));
      copy(sigma2x.at(0), sigma2cnm.slice(0,0,maxDegree+1,1));

      // remaining orders
      for(UInt m=1; m<=maxDegree; m++)
      {
        copy(x.at(2*m-1),       cnm.slice(m,m,maxDegree+1-m,1));
        copy(x.at(2*m+0),       snm.slice(m,m,maxDegree+1-m,1));
        copy(sigma2x.at(2*m-1), sigma2cnm.slice(m,m,maxDegree+1-m,1));
        copy(sigma2x.at(2*m+0), sigma2snm.slice(m,m,maxDegree+1-m,1));
      }

      return SphericalHarmonics(GM, R, cnm, snm, sigma2cnm, sigma2snm).get(maxDegree, minDegree);
    }
    return SphericalHarmonics();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GriddedData2PotentialCoefficients::buildNormals(UInt i)
{
  try
  {
    UInt start = blocking * i;
    if(start>=grid.points.size())
      return;
    UInt count = std::min(grid.points.size(), start+blocking)-start;

    // design matrix for each order
    std::vector<Matrix> At(2*maxDegree+1);
    At.at(0) = Matrix(maxDegree+1, count);
    for(UInt m=1; m<=maxDegree; m++)
    {
      At.at(2*m-1) = Matrix(maxDegree+1-m, count);
      At.at(2*m+0) = Matrix(maxDegree+1-m, count);
    }

    Vector l(count);
    for(UInt i=0; i<count; i++)
    {
      Double weight = sqrt(grid.areas.at(i+start)/(4*PI));

      Matrix Cnm, Snm;
      SphericalHarmonics::CnmSnm((1/R)*grid.points.at(i+start), maxDegree, Cnm, Snm);

      l(i) = weight * grid.values.at(0).at(i+start);

      Vector kn = kernel->inverseCoefficients(grid.points.at(i+start), maxDegree);
      kn *= weight * GM/R;
      Cnm(0,0) *= kn(0);
      for(UInt n=1; n<=maxDegree; n++)
      {
        Cnm.slice(n,0,1,n+1) *= kn(n);
        Snm.slice(n,1,1,n)   *= kn(n);
      }

      // zero order (only cnm)
      // ------------------------
      copy(Cnm.column(0), At.at(0).column(i));

      // other orders (m>0)
      // ------------------
      for(UInt m=1; m<=maxDegree; m++)
      {
        copy(Cnm.slice(m,m,maxDegree+1-m,1), At.at(2*m-1).column(i));
        copy(Snm.slice(m,m,maxDegree+1-m,1), At.at(2*m+0).column(i));
      }
    }

    // accumulate normals
    // ------------------
    for(UInt i=0; i<At.size(); i++)
    {
      rankKUpdate(1., At.at(i).trans(), N.at(i));
      matMult    (1., At.at(i), l,      n.at(i));
    }
    lPl += quadsum(l);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GriddedData2PotentialCoefficients::buildNormalsFast(UInt i)
{
  try
  {
    // legendre functions with kernel coefficients
    Vector3d p   = polar(lambda.at(0), phi.at(i), radius.at(i));
    Vector   kn  = kernel->inverseCoefficients(p, maxDegree);
    Matrix   Pnm = SphericalHarmonics::Pnm(Angle(PI/2-phi.at(i)), radius.at(i)/R, maxDegree);
    for(UInt n=0; n<=maxDegree; n++)
      Pnm.slice(n,0,1,n+1) *= GM/R*kn(n);
    Pnm.setType(Matrix::GENERAL);

    // right hand side
    // ---------------
    // order m=0
    Double c = 0;
    for(UInt k=0; k<lambda.size(); k++)
      c += grid.areas.at(k+i*lambda.size())/(4*PI) * grid.values.at(0).at(k+i*lambda.size());
    axpy(c, Pnm.column(0), this->n.at(0));

    // remaining orders
    for(UInt m=1; m<=maxDegree; m++)
    {
      Double c = 0, s = 0;
      for(UInt k=0; k<lambda.size(); k++)
      {
        c += grid.areas.at(k+i*lambda.size())/(4*PI) * grid.values.at(0).at(k+i*lambda.size()) * cossinm(k,2*m-1);
        s += grid.areas.at(k+i*lambda.size())/(4*PI) * grid.values.at(0).at(k+i*lambda.size()) * cossinm(k,2*m+0);
      }
      axpy(c, Pnm.slice(m,m,maxDegree+1-m,1), this->n.at(2*m-1));
      axpy(s, Pnm.slice(m,m,maxDegree+1-m,1), this->n.at(2*m+0));
    }

    // normals matrix
    // --------------
    // order m=0
    c = 0;
    for(UInt k=0; k<lambda.size(); k++)
      c += grid.areas.at(k+i*lambda.size())/(4*PI);
    rankKUpdate(c, Pnm.column(0).trans(), N.at(0));

    // remaining orders
    for(UInt m=1; m<=maxDegree; m++)
    {
      Double cc = 0, ss = 0;
      for(UInt k=0; k<lambda.size(); k++)
      {
        cc += grid.areas.at(k+i*lambda.size())/(4*PI) * cossinm(k,2*m-1) * cossinm(k,2*m-1);
        ss += grid.areas.at(k+i*lambda.size())/(4*PI) * cossinm(k,2*m+0) * cossinm(k,2*m+0);
      }
      rankKUpdate(cc, Pnm.slice(m,m,maxDegree+1-m,1).trans(), N.at(2*m-1));
      rankKUpdate(ss, Pnm.slice(m,m,maxDegree+1-m,1).trans(), N.at(2*m+0));
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

// Statt Mittelpunkte Integrale ueber die Bloecke
// ==============================================
// {
//   const Double deltaTheta  = 0.125*DEG2RAD;
//   const Double deltaLambda = 0.125*DEG2RAD;
//   Cnm = Matrix(maxDegree+1, Matrix::SYMMETRIC, Matrix::LOWER);
//   Snm = Matrix(maxDegree+1, Matrix::SYMMETRIC, Matrix::LOWER);
//
//   Double lambda = grid.points.at(i+start).lambda();
//   Double theta  = grid.points.at(i+start).theta();
//   Double Rr     = R/grid.points.at(i+start).r();
//   Double factor = 1/grid.areas.at(i+start);
//
//   if((theta >= (deltaTheta/2)) && (theta <= (PI-deltaTheta/2)))
//   {
//     Matrix Pnm = LegendreFunction::integral(cos(theta+deltaTheta/2), cos(theta-deltaTheta/2), maxDegree);
//     for(UInt n=0; n<=maxDegree; n++)
//     {
//       factor *= Rr;
//       Cnm(n,0) = factor * Pnm(n,0) * deltaLambda;
//       for(UInt m=1; m<=n; m++)
//       {
//         Cnm(n,m) = factor * Pnm(n,m) * 2.*cos(m*lambda)*sin(m*deltaLambda/2)/m;
//         Snm(n,m) = factor * Pnm(n,m) * 2.*sin(m*lambda)*sin(m*deltaLambda/2)/m;
//       }
//     }
//   }
//   else
//   {
//     Vector Pn = LegendrePolynomial::integral(cos(deltaTheta/2), maxDegree);
//     if(theta>PI-deltaTheta)
//       Rr *= -1;
//     for(UInt n=0; n<=maxDegree; n++)
//     {
//       factor  *=  Rr;
//       Cnm(n,0) = -factor * Pn(n) * 2*PI;
//     }
//   }
// }
// ======================================
