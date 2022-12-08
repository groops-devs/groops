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
#include "base/string.h"
#include "inputOutput/fileSinex.h"
#include "files/fileMatrix.h"
#include "files/fileNormalEquation.h"

/***** CLASS ***********************************/

/** @brief Convert SINEX to GROOPS normal equations.
* @ingroup programsConversionGroup */
class Sinex2Normals
{
  Vector readVector(const Sinex &sinex, UInt &dimension, std::vector<ParameterName> &parameterNames, const std::string &label) const;
  Matrix readMatrix(const Sinex &sinex, UInt dimension, const std::string &label) const;

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

    logInfo<<"read SINEX file"<<Log::endl;
    Sinex sinex;
    readFileSinex(inNameSinex, sinex);

    // dimension of system of equations
    // --------------------------------
    NormalEquationInfo info;
    UInt dimension = 0;
    if(sinex.header.size() > 65)
      dimension = static_cast<UInt>(String::toInt(sinex.header.substr(60, 5)));

    // SOLUTION/STATISTICS
    auto iter = std::find_if(sinex.blocks.begin(), sinex.blocks.end(), [&](const auto &b) {return b->label == "SOLUTION/STATISTICS";});
    if(iter != sinex.blocks.end())
      for(auto &line : (*iter)->lines)
      {
        if(String::startsWith(line, " NUMBER OF DEGREES OF FREEDOM")) info.observationCount = String::toDouble(line.substr(32, 22)) + dimension;
        if(String::startsWith(line, " NUMBER OF OBSERVATIONS"))       info.observationCount = String::toDouble(line.substr(32, 22));
        if(String::startsWith(line, " WEIGHTED SQUARE SUM OF O-C"))   info.lPl(0)           = String::toDouble(line.substr(32, 22));
      }

    Vector n  = readVector(sinex, dimension, info.parameterName, "SOLUTION/NORMAL_EQUATION_VECTOR");
    Vector x  = readVector(sinex, dimension, info.parameterName, "SOLUTION/ESTIMATE");
    Vector x0 = readVector(sinex, dimension, info.parameterName, "SOLUTION/APRIORI");
    Matrix N0 = readMatrix(sinex, dimension, "SOLUTION/MATRIX_APRIORI");
    Matrix N  = readMatrix(sinex, dimension, "SOLUTION/NORMAL_EQUATION_MATRIX");

    // try reconstruct missing information
    if(!x0.size())
      x0 = Vector(x.size());
    if(!N.size())
    {
      N = readMatrix(sinex, dimension, "SOLUTION/MATRIX_ESTIMATE");
      if(N0.size())
        N -= N0;
    }
    if(!n.size() && N.size())
      n = N * (x-x0);
    if(!info.observationCount)
      info.observationCount = x.size();
    if(!info.lPl(0))
      info.lPl = (x-x0).trans()*(N*(x-x0));

    // write output files
    // ------------------
    if(!outNameNormals.empty() && N.size())
    {
      logStatus<<"write unconstrained normal equations to <"<<outNameNormals<<">"<<Log::endl;
      logInfo<<"  unknown parameters: "<<N.columns()<<Log::endl;
      logInfo<<"  observations:       "<<info.observationCount<<Log::endl;
      writeFileNormalEquation(outNameNormals, info, N, n);
    }

    if(!outNameNormalsConstraint.empty() && N0.size())
    {
      NormalEquationInfo infoConstraint(info.parameterName);
        for(UInt i=0; i<N0.rows(); i++)
          if(N0(i,i))
            infoConstraint.observationCount++;
      logStatus<<"write normal equations of applied constraints to <"<<outNameNormalsConstraint<<">"<<Log::endl;
      logInfo<<"  unknown parameters: "<<N0.columns()<<Log::endl;
      logInfo<<"  observations:       "<<infoConstraint.observationCount<<Log::endl;
      writeFileNormalEquation(outNameNormalsConstraint, infoConstraint, N0, Vector(N0.rows()));
    }

    if(!outNameSolution.empty() && x.size())
    {
      logStatus<<"write solution vector to <"<<outNameSolution<<">"<<Log::endl;
      writeFileMatrix(outNameSolution, x);
    }

    if(!outNameSolutionApriori.empty() && x0.size())
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

Vector Sinex2Normals::readVector(const Sinex &sinex, UInt &dimension, std::vector<ParameterName> &parameterNames, const std::string &label) const
{
  try
  {
    auto iter = std::find_if(sinex.blocks.begin(), sinex.blocks.end(), [&](const auto &b) {return b->label == label;});
    if(iter == sinex.blocks.end())
      return Vector();

    if(dimension == 0)
      for(auto &line : (*iter)->lines)
        dimension = std::max(dimension, static_cast<UInt>(String::toInt(line.substr(1, 5))));

    parameterNames.resize(dimension);
    Vector x(dimension);
    for(auto &line : (*iter)->lines)
    {
      // *INDEX _TYPE_ CODE PT SOLN _REF_EPOCH__ UNIT S ___ESTIMATED_VALUE___ __STD_DEV__
      const UInt idx = static_cast<UInt>(String::toInt(line.substr(1, 5)))-1;
      const std::string parameterType = String::trim(line.substr(7, 6));
      x(idx) = String::toDouble(line.substr(47, 21));

      // spherical harmonics coefficients
      if(parameterType == "CN" || parameterType == "SN")
      {
        std::string type = "sphericalHarmonics."
                         + std::string(parameterType == "CN" ? "c_" : "s_")
                         + String::trim(line.substr(14, 4)) + "_" // degree (siteCode)
                         + String::trim(line.substr(22, 4));      // order (solutionId)
        parameterNames.at(idx) = ParameterName("", type);
        continue;
      }

      std::string object = String::trim(line.substr(14, 4)); // siteCode
      if(object == "----")
        object = "";
      if(parameterType.substr(0,3) != "SAT")
        object = String::lowerCase(object);

      std::string type = parameterType; // not all types implemented yet, see SINEX documentation;
      if(     parameterType == "STAX")   type = "position.x";
      else if(parameterType == "STAY")   type = "position.y";
      else if(parameterType == "STAZ")   type = "position.z";
      else if(parameterType == "VELX")   type = "velocity.x";
      else if(parameterType == "VELY")   type = "velocity.y";
      else if(parameterType == "VELZ")   type = "velocity.z";
      else if(parameterType == "XGC")    type = "geocenter.x";
      else if(parameterType == "YGC")    type = "geocenter.y";
      else if(parameterType == "ZGC")    type = "geocenter.z";
      else if(parameterType == "LOD")    type = "LOD";
      else if(parameterType == "UT")     type = "UT1";
      else if(parameterType == "XPO")    type = "polarMotion.xp";
      else if(parameterType == "YPO")    type = "polarMotion.yp";
      else if(parameterType == "XPOR")   type = "polarMotionRate.xp";
      else if(parameterType == "YPOR")   type = "polarMotionRate.yp";
      else if(parameterType == "NUT_X")  type = "nutation.X";
      else if(parameterType == "NUT_Y")  type = "nutation.Y";
      else if(parameterType == "NUTR_X") type = "nutationRate.X";
      else if(parameterType == "NUTR_Y") type = "nutationRate.Y";
      else if(parameterType == "SAT__X") type = "position.x";
      else if(parameterType == "SAT__Y") type = "position.y";
      else if(parameterType == "SAT__Z") type = "position.z";
      else if(parameterType == "SAT_VX") type = "velocity.x";
      else if(parameterType == "SAT_VY") type = "velocity.y";
      else if(parameterType == "SAT_VZ") type = "velocity.z";
      else if(parameterType == "SATA_X") type = "antennaCenterVariations.xOffset";
      else if(parameterType == "SATA_Y") type = "antennaCenterVariations.yOffset";
      else if(parameterType == "SATA_Z") type = "antennaCenterVariations.zOffset";

      parameterNames.at(idx) = ParameterName(object, type);
    }

    return x;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Matrix Sinex2Normals::readMatrix(const Sinex &sinex, UInt dimension, const std::string &label) const
{
  try
  {
    auto iter = std::find_if(sinex.blocks.begin(), sinex.blocks.end(), [&](const auto &b) {return String::startsWith(b->label, label);});
    if(iter == sinex.blocks.end())
      return Matrix();

    Matrix N(dimension, dimension);
    Bool isLower = FALSE;
    for(auto &line : (*iter)->lines)
    {
      const UInt i = static_cast<UInt>(String::toInt(line.substr(1, 5))) - 1;
      const UInt k = static_cast<UInt>(String::toInt(line.substr(7, 5))) - 1;
      isLower = isLower || (k < i);
      for(UInt l=0; l<3; l++)
        if(line.length() >= 13+l*22+21)
          N(i,k+l) += String::toDouble(line.substr(13+l*22, 21));
    }
    N.setType(Matrix::SYMMETRIC, isLower ? Matrix::LOWER : Matrix::UPPER);
    fillSymmetric(N);
    N.setType(Matrix::SYMMETRIC, Matrix::UPPER);

    // convert correlation or covariance matrix to normals
    if(String::endsWith((*iter)->label, "CORR"))
    {
      for(UInt i=1; i<dimension; i++)
        N.column(0, i-1) *= N(i,i);
      for(UInt i=0; i<dimension; i++)
        N.row(i, dimension-i) *= N(i,i);
    }
    if(String::endsWith((*iter)->label, "CORR") || String::endsWith((*iter)->label, "COVA"))
      inverse(N);

    return N;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
