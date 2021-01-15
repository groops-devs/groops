/***********************************************/
/**
* @file gravityfield2PotentialCoefficientsTimeSeries.cpp
*
* @brief Time series of potential coefficients.
*
* @author Torsten Mayer-Guerr
* @date 2020-10-02
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program computes a \configClass{timeSeries}{timeSeriesType}
of a time variable \configClass{gravityfield}{gravityfieldType}
and converts to coefficients of a spherical harmonics expansion.
The expansion is limited in the range between \config{minDegree}
and \config{maxDegree} inclusivly.
The coefficients are related to the reference radius~\config{R}
and the Earth gravitational constant \config{GM}.

The \configFile{outputfileTimeSeries}{instrument} contains the potential coefficients
as data columns for each epoch in the sequence given by
\configClass{numbering}{sphericalHarmonicsNumberingType}.
)";

/***********************************************/

#include "programs/program.h"
#include "base/sphericalHarmonics.h"
#include "files/fileInstrument.h"
#include "classes/gravityfield/gravityfield.h"
#include "classes/timeSeries/timeSeries.h"
#include "classes/sphericalHarmonicsNumbering/sphericalHarmonicsNumbering.h"

/***** CLASS ***********************************/

/** @brief Time series of potential coefficients.
* @ingroup programsGroup */
class Gravityfield2PotentialCoefficientsTimeSeries
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(Gravityfield2PotentialCoefficientsTimeSeries, SINGLEPROCESS, "time series of potential coefficients", Gravityfield, TimeSeries)

/***********************************************/

void Gravityfield2PotentialCoefficientsTimeSeries::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName                       fileNameOut;
    GravityfieldPtr                gravityfield;
    TimeSeriesPtr                  timeSeries;
    UInt                           minDegree, maxDegree;
    Double                         GM, R;
    SphericalHarmonicsNumberingPtr numbering;

    readConfig(config, "outputfileTimeSeries", fileNameOut,  Config::MUSTSET, "",  "instrument file (MISCVALUES)");
    readConfig(config, "gravityfield",         gravityfield, Config::MUSTSET, "",  "");
    readConfig(config, "timeSeries",           timeSeries,   Config::MUSTSET, "",  "");
    readConfig(config, "minDegree",            minDegree,    Config::MUSTSET, "0", "");
    readConfig(config, "maxDegree",            maxDegree,    Config::MUSTSET, "",  "");
    readConfig(config, "GM",                   GM,           Config::DEFAULT, STRING_DEFAULT_GM, "Geocentric gravitational constant");
    readConfig(config, "R",                    R,            Config::DEFAULT, STRING_DEFAULT_R,  "reference radius");
    readConfig(config, "numbering",            numbering,    Config::MUSTSET, "",  "numbering scheme");
    if(isCreateSchema(config)) return;

    const std::vector<Time>        times = timeSeries->times();
    std::vector<std::vector<UInt>> idxC, idxS;
    numbering->numbering(maxDegree, minDegree, idxC, idxS);

    Matrix A(times.size(), 1+numbering->parameterCount(maxDegree, minDegree));
    for(UInt idEpoch=0; idEpoch<times.size(); idEpoch++)
    {
      SphericalHarmonics harm = gravityfield->sphericalHarmonics(times.at(idEpoch), maxDegree, minDegree, GM, R);
      for(UInt n=minDegree; n<=maxDegree; n++)
      {
        if(idxC[n][0]!=NULLINDEX) A(idEpoch, 1+idxC[n][0]) = harm.cnm()(n, 0);
        for(UInt m=1; m<=n; m++)
        {
          if(idxC[n][m]!=NULLINDEX) A(idEpoch, 1+idxC[n][m]) = harm.cnm()(n, m);
          if(idxS[n][m]!=NULLINDEX) A(idEpoch, 1+idxS[n][m]) = harm.snm()(n, m);
        }
      }
    }

    logStatus<<"write time series of potential coefficients to file <"<<fileNameOut<<">"<<Log::endl;
    InstrumentFile::write(fileNameOut, Arc(times, A));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
