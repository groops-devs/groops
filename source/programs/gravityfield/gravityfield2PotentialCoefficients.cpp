/***********************************************/
/**
* @file gravityfield2PotentialCoefficients.cpp
*
* @brief Writes a gravity field to a file with potential coefficients.
*
* @author Torsten Mayer-Guerr
* @date 2002-05-18
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program evaluates a time variable \configClass{gravityfield}{gravityfieldType}
at a given \config{time} and saves it as a \file{spherical harmonics file}{potentialCoefficients}.
If set the expansion is limited in the range between \config{minDegree}
and \config{maxDegree} inclusivly.
The coefficients are related to the reference radius~\config{R}
and the Earth gravitational constant \config{GM}.
)";

/***********************************************/

#include "programs/program.h"
#include "files/fileSphericalHarmonics.h"
#include "classes/gravityfield/gravityfield.h"

/***** CLASS ***********************************/

/** @brief Writes a gravity field to a file with potential coefficients.
* @ingroup programsGroup */
class Gravityfield2PotentialCoefficients
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(Gravityfield2PotentialCoefficients, SINGLEPROCESS, "writes a gravity field to a file with potential coefficients", Gravityfield, PotentialCoefficients)

/***********************************************/

void Gravityfield2PotentialCoefficients::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName        fileNameCoeff;
    UInt            minDegree, maxDegree = INFINITYDEGREE;
    Time            time;
    Double          GM, R;
    GravityfieldPtr gravityfield;

    readConfig(config, "outputfilePotentialCoefficients", fileNameCoeff, Config::MUSTSET,  "",  "");
    readConfig(config, "gravityfield",                    gravityfield,  Config::MUSTSET,  "",  "");
    readConfig(config, "minDegree",                       minDegree,     Config::DEFAULT,  "0", "");
    readConfig(config, "maxDegree",                       maxDegree,     Config::OPTIONAL, "",  "");
    readConfig(config, "GM",                              GM,            Config::DEFAULT,  STRING_DEFAULT_GM, "Geocentric gravitational constant");
    readConfig(config, "R",                               R,             Config::DEFAULT,  STRING_DEFAULT_R,  "reference radius");
    readConfig(config, "time",                            time,          Config::OPTIONAL, "", "at this time the gravity field will be evaluated");
    if(isCreateSchema(config)) return;

    // Create potential coefficients
    // -----------------------------
    logStatus<<"create spherical harmonics"<<Log::endl;
    SphericalHarmonics harm = gravityfield->sphericalHarmonics(time, maxDegree, minDegree, GM, R);

    // write
    // -----
    logStatus<<"writing potential coefficients to file <"<<fileNameCoeff<<">"<<Log::endl;
    writeFileSphericalHarmonics(fileNameCoeff, harm);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
