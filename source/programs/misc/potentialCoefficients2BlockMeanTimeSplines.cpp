/***********************************************/
/**
* @file potentialCoefficients2BlockMeanTimeSplines.cpp
*
* @brief Write monthly potential coeffcients into one time spline file.
*
* @author Torsten Mayer-Guerr
* @date 2008-08-08
*
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program is a simplified version of \program{Gravityfield2TimeSplines}.
It reads a series of potential coefficient files (\configFile{inputfilePotentialCoefficients}{potentialCoefficients})
and creates a time splines file with spline degree 0 (temporal block means) or degree 1 (linear splines).
The time intervals in which the potential coefficients are valid are defined between adjacent
points in time given by \config{splineTimeSeries}. Therefore one more point in time is needed
than the number of potential coefficient files for degree 0.

The coefficients can be filtered with \configClass{filter}{sphericalHarmonicsFilterType}.
If set the expansion is limited in the range between \config{minDegree} and \config{maxDegree} inclusivly.
The coefficients are related to the reference radius~\config{R} and the Earth gravitational constant \config{GM}.

This program is useful e.g. to combine monthly GRACE solutions to one file.
)";

/***********************************************/

#include "programs/program.h"
#include "base/sphericalHarmonics.h"
#include "files/fileSphericalHarmonics.h"
#include "files/fileTimeSplinesGravityfield.h"
#include "classes/timeSeries/timeSeries.h"
#include "classes/sphericalHarmonicsFilter/sphericalHarmonicsFilter.h"

/***** CLASS ***********************************/

/** @brief Write monthly potential coeffcients into one time spline file.
* @ingroup programsGroup */
class PotentialCoefficients2BlockMeanTimeSplines
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(PotentialCoefficients2BlockMeanTimeSplines, SINGLEPROCESS, "write monthly potential coeffcients into one time spline file", Misc, PotentialCoefficients, TimeSplines)

/***********************************************/

void PotentialCoefficients2BlockMeanTimeSplines::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName              outputName, covName;
    std::vector<FileName> inputName;
    TimeSeriesPtr         timeSeries;
    SphericalHarmonicsFilterPtr filter;
    UInt                  minDegree, maxDegree = INFINITYDEGREE;
    Double                GM, R;
    Bool                  removeMean;
    Bool                  interpolate;
    UInt                  splineDegree;

    readConfig(config, "outputfileTimeSplines",           outputName, Config::MUSTSET,  "", "");
    readConfig(config, "outputfileTimeSplinesCovariance", covName,    Config::OPTIONAL, "", "only the variances are saved");
    readConfig(config, "inputfilePotentialCoefficients",  inputName,  Config::MUSTSET,  "{groopsDataDir}/potential/", "");
    readConfig(config, "filter",           filter,        Config::DEFAULT,  "",  "");
    readConfig(config, "minDegree",        minDegree,     Config::DEFAULT,  "0", "");
    readConfig(config, "maxDegree",        maxDegree,     Config::OPTIONAL, "",  "");
    readConfig(config, "GM",               GM,            Config::DEFAULT,  STRING_DEFAULT_GM, "Geocentric gravitational constant");
    readConfig(config, "R",                R,             Config::DEFAULT,  STRING_DEFAULT_R,  "reference radius");
    readConfig(config, "removeMean",       removeMean,    Config::DEFAULT,  "0", "remove the temporal mean of the series before estimating the splines");
    readConfig(config, "interpolate",      interpolate,   Config::DEFAULT,  "0", "interpolate missing files");
    readConfig(config, "splineTimeSeries", timeSeries,    Config::MUSTSET,  "",  "input files must be between points in time");
    readConfig(config, "splineDegree",     splineDegree,  Config::DEFAULT,  "0", "degree of splines");
    if(isCreateSchema(config)) return;

    std::vector<Time> times  = timeSeries->times();
    const UInt timeCount = times.size();
    const UInt fileCount = inputName.size();
    if(timeCount-1+splineDegree != fileCount)
      throw(Exception("fileCount("+fileCount%"%i) != timeCount("s+timeCount%"%i)-1+splineDegree"s));

    // read data
    // ---------
    logStatus<<"read potential coefficients from files"<<Log::endl;
    std::vector<Matrix> cnmList(fileCount), snmList(fileCount);
    std::vector<Matrix> sigma2List(fileCount);
    std::vector<Bool>   isZero(fileCount, FALSE);
    Single::forEach(fileCount, [&](UInt i)
    {
      SphericalHarmonics harm;
      try
      {
        readFileSphericalHarmonics(inputName.at(i), harm);
        harm = filter->filter(harm);
      }
      catch(std::exception &e)
      {
        logError<<e.what()<<": continue..."<<Log::endl;
        harm = SphericalHarmonics();
        isZero.at(i) = TRUE;
      }

      harm      = harm.get(maxDegree, minDegree, GM, R);
      maxDegree = harm.maxDegree();
      GM        = harm.GM();
      R         = harm.R();

      cnmList.at(i)    = harm.cnm();
      snmList.at(i)    = harm.snm();
      sigma2List.at(i) = harm.sigma2x();
    });

    // interpolate missing data
    // ------------------------
    if(interpolate)
    {
      logStatus<<"interpolate missing data"<<Log::endl;
      UInt idxStart = NULLINDEX;
      for(UInt i=0; i<fileCount; i++)
      {
        if(!isZero.at(i))
          idxStart = i;
        else
        {
          UInt idxEnd = NULLINDEX;
          for(UInt k=i+1; k<fileCount; k++)
            if(!isZero.at(k))
            {
              idxEnd = k;
              break;
            }
          if((idxStart==NULLINDEX)||(idxEnd==NULLINDEX))
            throw(Exception("Cannot interpolate data at begin or end"));
          Double tau = static_cast<Double>(i-idxStart)/(idxEnd-idxStart);
          cnmList.at(i)    = (1-tau) * cnmList.at(idxStart)    + tau * cnmList.at(idxEnd);
          snmList.at(i)    = (1-tau) * snmList.at(idxStart)    + tau * snmList.at(idxEnd);
          sigma2List.at(i) = (1-tau) * sigma2List.at(idxStart) + tau * sigma2List.at(idxEnd);
        }
      }
    }

    // remove mean
    // -----------
    if(removeMean)
    {
      logStatus<<"remove mean"<<Log::endl;
      UInt count = 0;
      for(UInt i=0; i<fileCount; i++)
        if(!isZero.at(i))
          count++;

      Matrix cnmMean = cnmList.at(0);
      for(UInt i=1; i<fileCount; i++)
        if(!isZero.at(i))
          cnmMean += cnmList.at(i);
      cnmMean *= 1./count;
      for(UInt i=0; i<fileCount; i++)
        if(interpolate || !isZero.at(i))
          cnmList.at(i) -= cnmMean;

      Matrix snmMean = snmList.at(0);
      for(UInt i=1; i<fileCount; i++)
        if(!isZero.at(i))
          snmMean += snmList.at(i);
      snmMean *= 1./count;
      for(UInt i=0; i<fileCount; i++)
        if(interpolate || !isZero.at(i))
          snmList.at(i) -= snmMean;
    }

    // write timeSplines file
    // ----------------------
    if(!outputName.empty())
    {
      logStatus<<"write time splines to file <"<<outputName<<">"<<Log::endl;
      writeFileTimeSplinesGravityfield(outputName, GM, R, splineDegree, times, cnmList, snmList);
    }

    if(!covName.empty())
    {
      logStatus<<"write covariance time splines to file <"<<covName<<">"<<Log::endl;
      writeFileTimeSplinesCovariance(covName, GM, R, minDegree, maxDegree, splineDegree, times, sigma2List);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
