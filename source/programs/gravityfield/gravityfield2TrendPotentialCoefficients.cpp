/***********************************************/
/**
* @file gravityfield2TrendPotentialCoefficients.cpp
*
* @brief Estimate temporal parametrization (e.g. trend, annual).
*
* @author Torsten Mayer-Guerr
* @date 2023-07-24
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program estimates \configClass{parametrizationTemporal}{parametrizationTemporalType}
(e.g. mean, trend, annual) from a time variable gravity field.

In a first step a time variable \configClass{gravityfield}{gravityfieldType}
is sampled at \configClass{timeSeries}{timeSeriesType}
and converted to coefficients of a spherical harmonics expansion.
The expansion is limited in the range between \config{minDegree}
and \config{maxDegree} inclusively.
The coefficients are related to the reference radius~\config{R}
and the Earth gravitational constant \config{GM}.

These coefficients serves as observations of
a \reference{robust least squares adjustment}{fundamentals.robustLeastSquares} to estimate
\configClass{parametrizationTemporal}{parametrizationTemporalType} parameters.
For each temporal parameter an
\configFile{outputfilePotentialCoefficients}{potentialCoefficients}
is generated.
)";

/***********************************************/

#include "programs/program.h"
#include "files/fileSphericalHarmonics.h"
#include "classes/gravityfield/gravityfield.h"
#include "classes/timeSeries/timeSeries.h"
#include "classes/parametrizationTemporal/parametrizationTemporal.h"
#include "misc/varianceComponentEstimation.h"

/***** CLASS ***********************************/

/** @brief Estimate temporal parametrization (e.g. trend, annual).
* @ingroup programsGroup */
class Gravityfield2TrendPotentialCoefficients
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(Gravityfield2TrendPotentialCoefficients, SINGLEPROCESS, "estimate temporal parametrization (e.g. trend, annual)", Gravityfield)

/***********************************************/

void Gravityfield2TrendPotentialCoefficients::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    std::vector<FileName>      fileNamesOut;
    GravityfieldPtr            gravityfield;
    TimeSeriesPtr              timeSeries;
    ParametrizationTemporalPtr temporal;
    UInt                       minDegree, maxDegree = INFINITYDEGREE;
    Double                     GM, R;
    Double                     huber, huberPower;
    UInt                       maxIter;

    readConfig(config, "outputfilePotentialCoefficients", fileNamesOut, Config::MUSTSET,  "",  "for each temporal parameter");
    readConfig(config, "gravityfield",                    gravityfield, Config::MUSTSET,  "",  "");
    readConfig(config, "timeSeries",                      timeSeries,   Config::MUSTSET,  "",  "");
    readConfig(config, "parametrizationTemporal",         temporal,     Config::MUSTSET,  "",  "");
    readConfig(config, "minDegree",                       minDegree,    Config::MUSTSET,  "0", "");
    readConfig(config, "maxDegree",                       maxDegree,    Config::OPTIONAL, "",  "");
    readConfig(config, "GM",                              GM,           Config::DEFAULT,  STRING_DEFAULT_GM, "Geocentric gravitational constant");
    readConfig(config, "R",                               R,            Config::DEFAULT,  STRING_DEFAULT_R,  "reference radius");
    readConfig(config, "huber",                           huber,        Config::DEFAULT,  "2.5", "for robust least squares");
    readConfig(config, "huberPower",                      huberPower,   Config::DEFAULT,  "1.5", "for robust least squares");
    readConfig(config, "huberMaxIteration",               maxIter,      Config::DEFAULT,  "15",  "(maximum) number of iterations for robust estimation");
    if(isCreateSchema(config)) return;

    const std::vector<Time> times = timeSeries->times();

    if(maxDegree == INFINITYDEGREE)
      maxDegree = gravityfield->sphericalHarmonics(times.front()).maxDegree();

    // sampling of gravity field (observation vectors)
    Matrix L(times.size(), (maxDegree+1)*(maxDegree+1));
    for(UInt idEpoch=0; idEpoch<times.size(); idEpoch++)
      copy(gravityfield->sphericalHarmonics(times.at(idEpoch), maxDegree, minDegree, GM, R).x().trans(), L.row(idEpoch));

    // design matrix
    temporal->setInterval(times.front(), times.back()+medianSampling(times), TRUE);
    Matrix A(times.size(), temporal->parameterCount());
    for(UInt idEpoch=0; idEpoch<times.size(); idEpoch++)
      copy(temporal->factors(times.at(idEpoch)).trans(), A.row(idEpoch));

    logStatus<<"estimate "<<temporal->parameterCount()<<" temporal parameters from "<<times.size()<<" epochs"<<Log::endl;
    Vector sigma;
    Matrix x = Vce::robustLeastSquares(A, L, 1, huber, huberPower, maxIter, sigma);
    for(UInt idEpoch=0; idEpoch<times.size(); idEpoch++)
      if(sigma(idEpoch) > 1)
        logInfo<<"  "<<times.at(idEpoch).dateTimeStr()<<" outlier with sigma = "<<sigma(idEpoch)<<Log::endl;


    for(UInt i=0; i<fileNamesOut.size(); i++)
    {
      logStatus<<"write potential coefficients to file <"<<fileNamesOut.at(i)<<">"<<Log::endl;
      Matrix cnm(maxDegree+1, Matrix::SYMMETRIC, Matrix::LOWER);
      Matrix snm(maxDegree+1, Matrix::SYMMETRIC, Matrix::LOWER);
      UInt idx = 0;
      for(UInt n=0; n<=maxDegree; n++)
      {
        cnm(n, 0) = x(i, idx++);
        for(UInt m=1; m<=n; m++)
        {
          cnm(n, m) = x(i, idx++);
          snm(n, m) = x(i, idx++);
        }
      }
      writeFileSphericalHarmonics(fileNamesOut.at(i), SphericalHarmonics(GM, R, cnm, snm));
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
