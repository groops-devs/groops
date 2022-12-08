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

    readConfig(config, "outputfileSinex",        fileNameSinex,    Config::MUSTSET,  "", "solutions in SINEX format");
    readConfig(config, "inputfileNormals",       fileNameNormals,  Config::MUSTSET,  "", "normal equation matrix");
    readConfig(config, "inputfileSolution",      fileNameSolution, Config::OPTIONAL, "", "parameter vector");
    readConfig(config, "inputfileSigmax",        fileNameSigmax,   Config::OPTIONAL, "", "standard deviations of the parameters (sqrt of the diagonal of the inverse normal equation)");
    readConfig(config, "inputfileApriori",       fileNameApriori,  Config::MUSTSET,  "", "apriori parameter vector");
    readConfig(config, "inputfileAprioriMatrix", fileNameAprMat,   Config::OPTIONAL, "", "normal equation matrix of applied constraints");
    readConfig(config, "time",                   time0,            Config::MUSTSET,  "", "reference time for parameters");
    readConfig(config, "sinexHeader",            sinex,            Config::MUSTSET,  "", "");
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
    fillSymmetric(N);
    UInt countParameter = N.rows();

    Matrix dN;
    std::vector<Bool> parameterIsConstrained(countParameter, FALSE);
    if(!fileNameAprMat.empty())
    {
      logStatus<<"reading normal equation matrix of applied constraints <"<<fileNameAprMat<<">"<<Log::endl;
      Vector n;
      NormalEquationInfo info;
      readFileNormalEquation(fileNameAprMat, info, dN, n);
      fillSymmetric(dN);
      if(dN.rows() != parameterIsConstrained.size())
        throw(Exception("Parameter count in constraint matrix and normal equation matrix differs ("+dN.rows()%"%i"s+" vs. "+N.rows()%"%i"s+" )."));
      for(UInt i=0; i<dN.rows(); i++)
        parameterIsConstrained.at(i) = (dN(i, i) != 0.0);
    }

    // ==================================================

    auto writeVector = [&](SinexBlockPtr block, const Vector x, const Vector sigma=Vector())
    {
      for(UInt i=0; i<x.size(); i++)
      {
        if(info.parameterName.at(i).type.substr(0,18) != "sphericalHarmonics")
          throw(Exception("non spherical harmonics parameter: " + info.parameterName.at(i).str()));

        const UInt idxDegree = info.parameterName.at(i).type.find_first_of('_')+1;
        const UInt idxOrder  = info.parameterName.at(i).type.find_last_of('_')+1;

        const std::string degree = std::atoi(info.parameterName.at(i).type.substr(idxDegree,idxOrder-idxDegree-1).c_str())%"%4i"s; // degree
        const std::string order  = std::atoi(info.parameterName.at(i).type.substr(idxOrder).c_str())%"%4i"s;                       // order

        *block<<(i+1)%" %5i "s<<((info.parameterName.at(i).type[idxDegree-2]=='c') ? "CN" : "SN")<<"     ";
        *block<<degree<<" -- "<<order<<" "<<Sinex::time2str(time0)<<" ---- "<<(parameterIsConstrained.at(i) ? "1" : "2")<<x(i)%" %21.14e"s;
        if(sigma.size())
          *block<<sigma(i)%" %11.4e"s;
        *block<<std::endl;
      }
    };

    auto writeMatrix = [&](SinexBlockPtr block, const Matrix &N)
    {
      for(UInt i=0; i<N.rows(); i++)
        for(UInt k=i; k<N.rows(); k++)
          if(N(i,k))
          {
            *block<<(i+1)%" %5i"s<<(k+1)%" %5i"s<<N(i, k)%" %21.14e"s;
            for(UInt l=1; (l<3) && (k+1<N.rows()) && N(i,k+1); l++, k++)
              *block<<N(i, k+1)%" %21.14e"s;
            *block<<std::endl;
          }
    };

    // ==================================================

    // add data to SINEX
    // -----------------
    {
      SinexBlockPtr block = sinex.addBlock("SOLUTION/STATISTICS");
      *block<<"*____STATISTICAL_PARAMETER_____ _______VALUE(S)_______"<<std::endl;
      *block<<" NUMBER OF OBSERVATIONS         "<<info.observationCount%"%22i"s<<std::endl;
      *block<<" NUMBER OF UNKNOWNS             "<<countParameter%"%22i"s<<std::endl;
      *block<<" NUMBER OF DEGREES OF FREEDOM   "<<(info.observationCount-countParameter)%"%22i"s<<std::endl;
      *block<<" WEIGHTED SQUARE SUM OF O-C     "<<info.lPl(0)%"%22.15e"s<<std::endl;
    }

    if(x.size())
    {
      SinexBlockPtr block = sinex.addBlock("SOLUTION/ESTIMATE");
      *block<<"*INDEX _TYPE_ CODE PT SOLN _REF_EPOCH__ UNIT S ___ESTIMATED_VALUE___ __STD_DEV__"<<std::endl;
      writeVector(block, x, sigmax);
    }

    if(x0.size())
    {
      SinexBlockPtr block = sinex.addBlock("SOLUTION/APRIORI");
      *block<<"*INDEX _TYPE_ CODE PT SOLN _REF_EPOCH__ UNIT S ____APRIORI_VALUE____ __STD_DEV__"<<std::endl;
      writeVector(block, x0, Vector(x0.size()));
    }

    if(n.size())
    {
      SinexBlockPtr block = sinex.addBlock("SOLUTION/NORMAL_EQUATION_VECTOR");
      *block<<"*INDEX _TYPE_ CODE PT SOLN _REF_EPOCH__ UNIT S ___RIGHT_HAND_SIDE___"<<std::endl;
      writeVector(block, n);
    }

    if(N.size())
    {
      SinexBlockPtr block = sinex.addBlock("SOLUTION/NORMAL_EQUATION_MATRIX U");
      *block<<"*PARA1 PARA2 _______PARA2+0_______ _______PARA2+1_______ _______PARA2+2_______"<<std::endl;
      writeMatrix(block, N);
    }

    if(dN.size())
    {
      SinexBlockPtr block = sinex.addBlock("SOLUTION/MATRIX_APRIORI U INFO");
      *block<<"*PARA1 PARA2 _______PARA2+0_______ _______PARA2+1_______ _______PARA2+2_______"<<std::endl;
      writeMatrix(block, dN);
    }

    logStatus<<"write SINEX file <"<<fileNameSinex<<">"<<Log::endl;
    if(sinex.header.size() > 65)
      sinex.header.replace(60, 5, countParameter%"%05i"s);
    writeFileSinex(fileNameSinex, sinex);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
