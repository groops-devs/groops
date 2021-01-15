/***********************************************/
/**
* @file iers2OceanPoleTide.cpp
*
* @brief Read Ocean pole tide model.
*
* @author Torsten Mayer-Guerr
* @date 2007-04-08
*
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
Read ocean pole tide model according to IERS conventions
and convert into \file{oceanPoleTide file}{oceanPoleTide}.
)";

/***********************************************/

#include "programs/program.h"
#include "inputOutput/file.h"
#include "files/fileOceanPoleTide.h"
#include "files/fileMatrix.h"

/***** CLASS ***********************************/

/** @brief Read ocean pole tide model.
* @ingroup programsConversionGroup */
class Iers2OceanPoleTide
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(Iers2OceanPoleTide, SINGLEPROCESS, "Read ocean pole tide model", Conversion)

/***********************************************/

void Iers2OceanPoleTide::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName outputName, inputName, loveNumberName;
    Double GM, R, Omega;
    Double rho, G, g;
    UInt   maxDegree;

    readConfig(config, "outputfileOceanPole",        outputName,     Config::MUSTSET, "",                "");
    readConfig(config, "inputfile",                  inputName,      Config::MUSTSET, "",                "");
    readConfig(config, "inputfileLoadingLoveNumber", loveNumberName, Config::MUSTSET, "{groopsDataDir}/loading/loadLoveNumbers_Gegout97.txt", "");
    readConfig(config, "maxDegree",                  maxDegree,      Config::MUSTSET, "",                "");
    readConfig(config, "GM",                         GM,             Config::DEFAULT, STRING_DEFAULT_GM, "Geocentric gravitational constant");
    readConfig(config, "R",                          R,              Config::DEFAULT, STRING_DEFAULT_R,  "Reference radius");
    readConfig(config, "Omega",                      Omega,          Config::DEFAULT, "7.292115e-05",    "[rad/s] earth rotation");
    readConfig(config, "rho",                        rho,            Config::DEFAULT, "1025",            "[kg/m**3] density of sea water");
    readConfig(config, "G",                          G,              Config::DEFAULT, "6.673e-11",       "[m**3/(kg*s**2)] gravitational constant");
    readConfig(config, "g",                          g,              Config::DEFAULT, "9.7803278",       "[m/s**2] gravity");
    if(isCreateSchema(config)) return;

    // Loading Love Numbers
    // --------------------
    Vector kn;
    readFileMatrix(loveNumberName, kn);

    // read potential coefficients
    // ---------------------------
    logStatus<<"read file <"<<inputName<<">"<<Log::endl;
    Matrix cnmReal(maxDegree+1, Matrix::TRIANGULAR, Matrix::LOWER);
    Matrix snmReal(maxDegree+1, Matrix::TRIANGULAR, Matrix::LOWER);
    Matrix cnmImag(maxDegree+1, Matrix::TRIANGULAR, Matrix::LOWER);
    Matrix snmImag(maxDegree+1, Matrix::TRIANGULAR, Matrix::LOWER);
    {
      InFile file(inputName);
      std::string line;
      std::getline(file, line); // read header
      while(file.good())
      {
        std::getline(file, line);
        if(line.empty())
          continue;
        std::stringstream ss(line);
        ss.exceptions(std::ios::badbit | std::ios::failbit);
        UInt n,m;
        ss>>n>>m;
        if(n<=maxDegree)
          ss>>cnmReal(n,m)>>snmReal(n,m)>>cnmImag(n,m)>>snmImag(n,m);
      }
    }

    // apply factors
    // -------------
    for(UInt n=0; n<=maxDegree; n++)
    {
      const Double f = pow(Omega,2)*pow(R,4)*4*PI*G*rho/GM/g * (1.+kn(n))/(2.*n+1.);
      cnmReal.row(n) *= f;
      snmReal.row(n) *= f;
      cnmImag.row(n) *= f;
      snmImag.row(n) *= f;
    }

    // write result
    // ------------
    logStatus<<"write ocean pole tide <"<<outputName<<">"<<Log::endl;
    writeFileOceanPoleTide(outputName, SphericalHarmonics(GM, R, cnmReal, snmReal), SphericalHarmonics(GM, R, cnmImag, snmImag));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
