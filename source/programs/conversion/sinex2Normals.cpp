/***********************************************/
/**
* @file sinex2Normals.cpp
*
* @brief Convert SINEX to GROOPS normal equations.
*
* @author Sebastian Strasser
* @date 2017-05-16
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
Convert normal equations from \href{http://www.iers.org/IERS/EN/Organization/AnalysisCoordinator/SinexFormat/sinex.html}{SINEX format}
to \file{normal equations}{normalEquation}.

See also \program{GnssNormals2Sinex} and \program{NormalsSphericalHarmonics2Sinex}.
)";

/***********************************************/

#include "programs/program.h"
#include "inputOutput/fileSinex.h"
#include "files/fileMatrix.h"
#include "files/fileNormalEquation.h"

/***** CLASS ***********************************/

/** @brief Convert SINEX to GROOPS normal equations.
* @ingroup programsConversionGroup */
class Sinex2Normals
{
  void correlation2covariance(MatrixSliceRef correlationMatrix) const;
  Matrix normalsMatrix(Sinex::SinexSolutionMatrixPtr sinexSolutionMatrix) const;

public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(Sinex2Normals, SINGLEPROCESS, "Convert SINEX to GROOPS normal equations.", Conversion, NormalEquation)

/***********************************************/

void Sinex2Normals::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName outNameNormals, outNameNormalsConstraint, outNameSolutionApriori, outNameSolution, inNameSinex;

    readConfig(config, "outputfileNormals",           outNameNormals,           Config::OPTIONAL, "", "N, n: unconstrained normal equations");
    readConfig(config, "outputfileNormalsConstraint", outNameNormalsConstraint, Config::OPTIONAL, "", "N0, n0: normal equations of applied constraints");
    readConfig(config, "outputfileSolution",          outNameSolution,          Config::OPTIONAL, "", "x: parameter vector");
    readConfig(config, "outputfileSolutionApriori",   outNameSolutionApriori,   Config::OPTIONAL, "", "x0: a priori parameter vector");
    readConfig(config, "inputFileSinex",              inNameSinex,              Config::MUSTSET,  "", "");
    if(isCreateSchema(config)) return;

    logInfo << "read SINEX file" << Log::endl;
    Sinex sinex(inNameSinex);

    // convert SINEX to normal equations
    // ---------------------------------
    NormalEquationInfo info, infoConstraint;
    Vector n, n0;
    Matrix N, N0;
    Vector x0 = sinex.getBlock<Sinex::SinexSolutionVector>("SOLUTION/APRIORI")->vector();
    Vector x;
    if(sinex.hasBlock("SOLUTION/ESTIMATE"))
      x = sinex.getBlock<Sinex::SinexSolutionVector>("SOLUTION/ESTIMATE")->vector();

    // get lPl and obsCount from solution statistics if available
    if(sinex.hasBlock("SOLUTION/STATISTICS"))
    {
      info.lPl    = Vector(1);
      info.lPl(0) = sinex.getBlock<Sinex::SinexSolutionStatistics>("SOLUTION/STATISTICS")->value("WEIGHTED SQUARE SUM OF O-C");
      info.observationCount = static_cast<UInt>(sinex.getBlock<Sinex::SinexSolutionStatistics>("SOLUTION/STATISTICS")->value("NUMBER OF OBSERVATIONS"));
    }

    // normal equations storage method 6b and 6c from SINEX documentation
    if(sinex.hasBlock("SOLUTION/NORMAL_EQUATION_MATRIX") && sinex.hasBlock("SOLUTION/NORMAL_EQUATION_VECTOR") && sinex.hasBlock("SOLUTION/STATISTICS"))
    {
      // unconstrained normal equations
      n = sinex.getBlock<Sinex::SinexSolutionVector>("SOLUTION/NORMAL_EQUATION_VECTOR")->vector();
      N = sinex.getBlock<Sinex::SinexSolutionMatrix>("SOLUTION/NORMAL_EQUATION_MATRIX")->matrix();

      // normal equations of applied constraints, if available
      if(sinex.hasBlock("SOLUTION/MATRIX_APRIORI") && sinex.hasBlock("SOLUTION/ESTIMATE"))
      {
        N0 = normalsMatrix(sinex.getBlock<Sinex::SinexSolutionMatrix>("SOLUTION/MATRIX_APRIORI"));
        n0 = Vector(N0.rows());
        infoConstraint.lPl = Vector(1, 0.0);
        infoConstraint.observationCount = 0;
        infoConstraint.blockIndex = {0, n0.rows()};
        infoConstraint.parameterName = sinex.getBlock<Sinex::SinexSolutionVector>("SOLUTION/NORMAL_EQUATION_VECTOR")->parameterNames();
        for(UInt i = 0; i < N0.rows(); i++)
          if(N0(i,i)!=0)
            infoConstraint.observationCount++;
      }

      info.parameterName = sinex.getBlock<Sinex::SinexSolutionVector>("SOLUTION/NORMAL_EQUATION_VECTOR")->parameterNames();
    }
    // normal equations storage method 6a from SINEX documentation
    else if(sinex.hasBlock("SOLUTION/MATRIX_ESTIMATE") && sinex.hasBlock("SOLUTION/ESTIMATE") && sinex.hasBlock("SOLUTION/APRIORI"))
    {
      // constrained normal equations
      N = normalsMatrix(sinex.getBlock<Sinex::SinexSolutionMatrix>("SOLUTION/MATRIX_ESTIMATE"));
      n = N * (x-x0);

      // if solution statistics are not provided, "restore" lPl and obsCount
      if(!sinex.hasBlock("SOLUTION/STATISTICS"))
      {
        info.lPl = x.trans()*(N*x);
        info.observationCount = x.size();
      }

      // normal equations of applied constraints, if available
      if(sinex.hasBlock("SOLUTION/MATRIX_APRIORI"))
      {
        N0 = normalsMatrix(sinex.getBlock<Sinex::SinexSolutionMatrix>("SOLUTION/MATRIX_APRIORI"));
        n0 = N0 * (x-x0);
        n -= n0;
        N -= N0;
      }
      else
        logWarning << "Normal equations of applied constraints not available. Normal equations may be constrained." << Log::endl;

      info.parameterName = sinex.getBlock<Sinex::SinexSolutionVector>("SOLUTION/ESTIMATE")->parameterNames();
    }
    else
      throw(Exception("SINEX file does not contain normal equations in any recognized storage method"));

    // write output files
    // ------------------
    if(!outNameNormals.empty() && N.size() && n.size() && info.lPl.size() && info.observationCount > 0)
    {
      logStatus<<"write unconstrained normal equations to <"<<outNameNormals<<">"<<Log::endl;
      logInfo<<"  unknown parameters: "<<N.columns()<<Log::endl;
      logInfo<<"  right hand sides:   "<<n.columns()<<Log::endl;
      logInfo<<"  observations:       "<<info.observationCount<<Log::endl;
      writeFileNormalEquation(outNameNormals, info, N.isUpper() ? N : N.trans(), n);
    }
    if(!outNameNormalsConstraint.empty() && N0.size() && n0.size() && infoConstraint.lPl.size() && infoConstraint.observationCount > 0)
    {
      logStatus<<"write normal equations of applied constraints to <"<<outNameNormalsConstraint<<">"<<Log::endl;
      logInfo<<"  unknown parameters: "<<N0.columns()<<Log::endl;
      logInfo<<"  right hand sides:   "<<n0.columns()<<Log::endl;
      logInfo<<"  observations:       "<<infoConstraint.observationCount<<Log::endl;
      writeFileNormalEquation(outNameNormalsConstraint, infoConstraint, N0.isUpper() ? N0 : N0.trans(), n0);
    }
    if(!outNameSolution.empty())
    {
      logStatus<<"write solution vector to <"<<outNameSolution<<">"<<Log::endl;
      writeFileMatrix(outNameSolution, x);
    }
    if(!outNameSolutionApriori.empty())
    {
      logStatus<<"write a priori solution vector to <"<<outNameSolutionApriori<<">"<<Log::endl;
      writeFileMatrix(outNameSolutionApriori, x0);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void Sinex2Normals::correlation2covariance(MatrixSliceRef A) const
{
  try
  {
    const Bool isUpper = A.isUpper();
    const UInt size    = A.rows();
    for(UInt i = 0; i < size; i++)
    {
      // off diagonal elements; r_ij -> sigma_ij
      for(UInt j = i+1; j < size; j++)
      {
        if(isUpper)
          A(i,j) *= A(i,i) * A(j,j);
        else
          A(j,i) *= A(i,i) * A(j,j);
      }
      // main diagonal elements: sigma_ii -> sigma_ii^2
      A(i,i) *= A(i,i);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Matrix Sinex2Normals::normalsMatrix(Sinex::SinexSolutionMatrixPtr sinexSolutionMatrix) const
{
  try
  {
    Matrix N = sinexSolutionMatrix->matrix();
    if(sinexSolutionMatrix->type() == Sinex::SinexSolutionMatrix::CORRELATION)
      correlation2covariance(N);
    if(sinexSolutionMatrix->type() == Sinex::SinexSolutionMatrix::CORRELATION || sinexSolutionMatrix->type() == Sinex::SinexSolutionMatrix::COVARIANCE)
      inverse(N);
    return N;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
