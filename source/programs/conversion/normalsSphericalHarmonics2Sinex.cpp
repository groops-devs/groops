/***********************************************/
/**
* @file normalsSphericalHarmonics2Sinex.cpp
*
* @brief Write potential coefficients and normal equations to SINEX format.
*
* @author Saniya Behzadpour
* @author Sebastian Strasser
* @date 2015-03-17
*
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
Write potential coefficients and \file{normal equations}{normalEquation} to
\href{http://www.iers.org/IERS/EN/Organization/AnalysisCoordinator/SinexFormat/sinex.html}{SINEX format}.

See also \program{Sinex2Normals} and \program{GnssNormals2Sinex}.
)";

/***********************************************/

#include "programs/program.h"
#include "inputOutput/fileSinex.h"
#include "files/fileMatrix.h"
#include "files/fileNormalEquation.h"

/***** CLASS ***********************************/

/** @brief Write potential coefficients and normal equations to SINEX format.
* @ingroup programsConversionGroup */
class NormalsSphericalHarmonics2Sinex
{
  static void addVector(Sinex::SinexSolutionVectorPtr vector, const Time &time, const std::vector<ParameterName> &parameterName, const Vector x, const Vector sigma, const std::vector<Bool> &parameterIsContrained);

public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(NormalsSphericalHarmonics2Sinex, SINGLEPROCESS, "Write potential coefficients and normal equations to SINEX format.", Conversion, NormalEquation)

/***********************************************/

void NormalsSphericalHarmonics2Sinex::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName    fileNameSinex;
    FileName    fileNameNormals;
    FileName    fileNameSolution, fileNameSigmax, fileNameApriori, fileNameAprMat;
    Time        time0;
    Sinex       sinex;

    readConfig(config, "outputfileSinex",         fileNameSinex,    Config::MUSTSET,  "", "solutions in SINEX format");
    readConfig(config, "inputfileNormals",        fileNameNormals,  Config::MUSTSET,  "", "normal equation matrix");
    readConfig(config, "inputfileSolution",       fileNameSolution, Config::OPTIONAL, "", "parameter vector");
    readConfig(config, "inputfileSigmax",         fileNameSigmax,   Config::OPTIONAL, "", "standard deviations of the parameters (sqrt of the diagonal of the inverse normal equation)");
    readConfig(config, "inputfileApriori",        fileNameApriori,  Config::MUSTSET,  "", "apriori parameter vector");
    readConfig(config, "inputfileAprioriMatrix",  fileNameAprMat,   Config::OPTIONAL, "", "normal equation matrix of applied constraints");
    readConfig(config, "time",                    time0,            Config::MUSTSET,  "", "reference time for parameters");
    sinex.readConfigHeader(config);
    if(isCreateSchema(config)) return;

    // ==================================================

    // read data from files
    // --------------------
    Matrix x;
    if(!fileNameSolution.empty())
    {
      logStatus<<"reading solution from <"<<fileNameSolution<<">"<<Log::endl;
      readFileMatrix(fileNameSolution, x);
    }

    Matrix sigmax;
    if(!fileNameSigmax.empty())
    {
      logStatus<<"reading standard deviations from <"<<fileNameSigmax<<">"<<Log::endl;
      readFileMatrix(fileNameSigmax, sigmax);
    }

    Matrix x0;
    if(!fileNameApriori.empty())
    {
      logStatus<<"reading apriori solution from <"<<fileNameApriori<<">"<<Log::endl;
      readFileMatrix(fileNameApriori, x0);
    }

    Matrix N, n;
    NormalEquationInfo info;
    logStatus<<"reading normal equation matrix from <"<<fileNameNormals<<">"<<Log::endl;
    readFileNormalEquation(fileNameNormals, info, N, n);
    UInt countParameter = N.rows();

    Matrix dN;
    std::vector<Bool> parameterIsConstrained(countParameter, FALSE);
    if(!fileNameAprMat.empty())
    {
      logStatus<<"reading normal equation matrix of applied constraints <"<<fileNameAprMat<<">"<<Log::endl;
      Vector n;
      NormalEquationInfo info;
      readFileNormalEquation(fileNameAprMat, info, dN, n);
      if(dN.rows() != parameterIsConstrained.size())
        throw(Exception("Parameter count in constraint matrix and normal equation matrix differs (" + dN.rows()%"%i"s + " vs. "+ N.rows()%"%i"s +" )."));
      for(UInt i = 0; i < dN.rows(); i++)
        parameterIsConstrained.at(i) = (dN(i, i) != 0.0);
    }

    // ==================================================

    // add data to SINEX
    // -----------------
    // SOLUTION/STATISTICS
    Sinex::SinexSolutionStatisticsPtr solutionStatistics = sinex.addBlock<Sinex::SinexSolutionStatistics>("SOLUTION/STATISTICS");
    solutionStatistics->addValue("NUMBER OF OBSERVATIONS", info.observationCount);
    solutionStatistics->addValue("NUMBER OF UNKNOWNS", countParameter);
    solutionStatistics->addValue("NUMBER OF DEGREES OF FREEDOM", info.observationCount-countParameter);
    solutionStatistics->addValue("WEIGHTED SQUARE SUM OF O-C", info.lPl(0));

    // SOLUTION/ESTIMATE
    if(x.size())
    {
      Sinex::SinexSolutionVectorPtr solutionEstimate = sinex.addBlock<Sinex::SinexSolutionVector>("SOLUTION/ESTIMATE");
      addVector(solutionEstimate, time0, info.parameterName, x, sigmax.size() ? sigmax : Vector(), parameterIsConstrained);
    }

    // SOLUTION/APRIORI
    if(x0.size())
    {
      Sinex::SinexSolutionVectorPtr solutionApriori = sinex.addBlock<Sinex::SinexSolutionVector>("SOLUTION/APRIORI");
      addVector(solutionApriori, time0, info.parameterName, x0, Vector(), parameterIsConstrained);
    }

    // SOLUTION/NORMAL_EQUATION_VECTOR
    Sinex::SinexSolutionVectorPtr solutionNormalEquationVector = sinex.addBlock<Sinex::SinexSolutionVector>("SOLUTION/NORMAL_EQUATION_VECTOR");
    addVector(solutionNormalEquationVector, time0, info.parameterName, n, Vector(), parameterIsConstrained);

    // SOLUTION/NORMAL_EQUATION_MATRIX
    Sinex::SinexSolutionMatrixPtr solutionNormalEquationMatrix = sinex.addBlock<Sinex::SinexSolutionMatrix>("SOLUTION/NORMAL_EQUATION_MATRIX " + std::string(N.isUpper() ? "U" : "L"));
    solutionNormalEquationMatrix->setMatrix(N);

    // SOLUTION/MATRIX_APRIORI
    if(dN.size())
    {
      Sinex::SinexSolutionMatrixPtr solutionNormalAprioriMatrix = sinex.addBlock<Sinex::SinexSolutionMatrix>("SOLUTION/MATRIX_APRIORI " + std::string(dN.isUpper() ? "U" : "L") + " INFO");
      solutionNormalAprioriMatrix->setMatrix(dN);
    }

    // ==================================================

    // write SINEX file
    // ----------------
    logStatus<<"write SINEX file <"<<fileNameSinex<<">"<<Log::endl;
    sinex.writeFile(fileNameSinex);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void NormalsSphericalHarmonics2Sinex::addVector(Sinex::SinexSolutionVectorPtr vector, const Time &time0, const std::vector<ParameterName> &parameterName, const Vector x, const Vector sigma, const std::vector<Bool> &parameterIsContrained)
{
  try
  {
    for(UInt i = 0; i < x.size(); i++)
    {
      const UInt idxDegree = parameterName.at(i).type.find_first_of('_')+1;
      const UInt idxOrder  = parameterName.at(i).type.find_last_of('_')+1;
      if(parameterName.at(i).type.substr(0,18) != "sphericalHarmonics")
        throw(Exception("non spherical harmonics parameter: " + parameterName.at(i).str()));

      Sinex::Parameter parameter;
      if(parameterName.at(i).type[idxDegree-2]=='c')
        parameter.parameterType = "CN";
      else if(parameterName.at(i).type[idxDegree-2]=='s')
        parameter.parameterType = "SN";
      else
        throw(Exception("unknown parameter type: " + parameterName.at(i).str()));
      parameter.parameterIndex = i+1;
      parameter.siteCode       = std::atoi(parameterName.at(i).type.substr(idxDegree,idxOrder-idxDegree-1).c_str())%"% 4i"s; // degree
      parameter.solutionId     = std::atoi(parameterName.at(i).type.substr(idxOrder).c_str())%"% 4i"s;                       // order
      parameter.pointCode      = "--";
      parameter.unit           = "----";
      parameter.constraintCode = parameterIsContrained.at(i) ? "1" : "2";
      parameter.time           = time0;
      parameter.value          = x(i);
      if(sigma.size())
        parameter.sigma        = sigma(i);

      vector->addParameter(parameter);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
