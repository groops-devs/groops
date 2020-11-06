/***********************************************/
/**
* @file gravityfield2SphericalHarmonicsVector.cpp
*
* @brief Converts a gravity field to a solution vector with potential coeffcients.
*
* @author Torsten Mayer-Guerr
* @date 2002-05-18
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program evaluates a time variable \configClass{gravityfield}{gravityfieldType} at a given \config{time}
and saves a \file{vector}{matrix} with the coefficients of a spherical harmonics expansion in the sequence given by
\configClass{numbering}{sphericalHarmonicsNumberingType}.
If set the expansion is limited in the range between \config{minDegree} and \config{maxDegree} inclusively.
The coefficients are related to the reference radius~\config{R} and the Earth gravitational constant \config{GM}.

This coefficients vector can be used as a approximate solution, see \program{NormalsMultiplyAdd},
or as pseudo oberservations for regularization,
see \configClass{normalEquation:regularization}{normalEquationType:regularization}.

For back transformation use \program{Gravityfield2PotentialCoefficients}
with \configClass{gravityfield:fromParametrization}{gravityfieldType:fromParametrization}.
)";

/***********************************************/

#include "programs/program.h"
#include "files/fileMatrix.h"
#include "classes/gravityfield/gravityfield.h"
#include "classes/sphericalHarmonicsNumbering/sphericalHarmonicsNumbering.h"

/***** CLASS ***********************************/

/** @brief Converts a gravity field to a solution vector with potential coeffcients.
* @ingroup programsGroup */
class Gravityfield2SphericalHarmonicsVector
{
public:
  void run(Config &config);
};

GROOPS_REGISTER_PROGRAM(Gravityfield2SphericalHarmonicsVector, SINGLEPROCESS, "converts a gravity field to a solution vector with potential coeffcients", Gravityfield, Matrix)

/***********************************************/

void Gravityfield2SphericalHarmonicsVector::run(Config &config)
{
  try
  {
    FileName        outName;
    UInt            startIndex;
    Time            time;
    Double          GM, R;
    UInt            minDegree, maxDegree;
    Bool            useSigma;
    GravityfieldPtr gravityfield;
    SphericalHarmonicsNumberingPtr numbering;

    readConfig(config, "outputfileVector", outName,  Config::MUSTSET,  "",  "");
    readConfig(config, "gravityfield", gravityfield, Config::MUSTSET,  "",  "");
    readConfig(config, "startIndex",   startIndex,   Config::DEFAULT,  "0", "start index to put the coefficents in the solution vector");
    readConfig(config, "minDegree",    minDegree,    Config::MUSTSET,  "2", "");
    readConfig(config, "maxDegree",    maxDegree,    Config::MUSTSET,  "",  "");
    readConfig(config, "GM",           GM,           Config::DEFAULT,  STRING_DEFAULT_GM, "Geocentric gravitational constant");
    readConfig(config, "R",            R,            Config::DEFAULT,  STRING_DEFAULT_R,  "reference radius");
    readConfig(config, "numbering",    numbering,    Config::MUSTSET,  "",  "numbering scheme for solution vector");
    readConfig(config, "time",         time,         Config::OPTIONAL, "",  "at this time the gravity field will be evaluated");
    readConfig(config, "useSigma",     useSigma,     Config::DEFAULT,  "0", "use formal errors instead of coefficients");
    if(isCreateSchema(config)) return;

    // create potential coefficients
    // -----------------------------
    logStatus<<"create spherical harmonics"<<Log::endl;
    SphericalHarmonics harm = gravityfield->sphericalHarmonics(time, maxDegree, minDegree, GM, R);

    // Create vector
    // -------------
    logStatus<<"sorting potential coefficients into vector"<<Log::endl;
    UInt rhsCount = 1; // right hand sides
    std::vector< std::vector<UInt> > idxC, idxS;
    numbering->numbering(maxDegree, minDegree, idxC, idxS);
    UInt dim = numbering->parameterCount(maxDegree, minDegree) + startIndex;
    Matrix x(dim,rhsCount);
    logInfo<<"dimension of the vector: "<<dim<<Log::endl;

    Matrix cnm = (useSigma) ? harm.sigma2cnm() : harm.cnm();
    Matrix snm = (useSigma) ? harm.sigma2snm() : harm.snm();

    for(UInt j=0; j<rhsCount; j++)
      for(UInt n=minDegree; n<=maxDegree; n++)
      {
        if(idxC[n][0]!=NULLINDEX) x(idxC[n][0]+startIndex,j) = (useSigma) ? std::sqrt(cnm(n,0)) : cnm(n,0);
        for(UInt m=1; m<=n; m++)
        {
          if(idxC[n][m]!=NULLINDEX) x(idxC[n][m]+startIndex,j) = (useSigma) ? std::sqrt(cnm(n,m)) : cnm(n,m);
          if(idxS[n][m]!=NULLINDEX) x(idxS[n][m]+startIndex,j) = (useSigma) ? std::sqrt(snm(n,m)) : snm(n,m);
        }
      }

    // Write
    // -----
    logStatus<<"write vector to file <"<<outName<<">"<<Log::endl;
    writeFileMatrix(outName, x);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
