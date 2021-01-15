/***********************************************/
/**
* @file graceCoefficients2BlockMeanTimeSplines.cpp
*
* @brief read GRACE data.
*
* @author Andreas Kvas
* @date 2019-09-19
*
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program converts potential coefficients from the GRACE SDS RL06 format
into \configFile{outputfileTimeSplines}{timeSplinesGravityField}.

The \configFile{outputfileTimeSeries}{instrument} contains the mid points
of non-empty intervals and \configFile{outputfileTimeIntervals}{instrument}
contains the monthly interval boundaries from first to last solution.

The output will always be monthly block means. If the SDS solutions do vary or overlap,
the nearest solution in terms of reference epoch is used.
)";

/***********************************************/

#include "programs/program.h"
#include "base/string.h"
#include "inputOutput/file.h"
#include "files/fileInstrument.h"
#include "files/fileTimeSplinesGravityfield.h"

/***** CLASS ***********************************/

/** @brief Read GRACE data.
* @ingroup programsConversionGroup */
class GraceCoefficients2BlockMeanTimeSplines
{
private:
 SphericalHarmonics readPotentialCoefficients(const FileName &name, Time &timeStart, Time &timeEnd);

public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(GraceCoefficients2BlockMeanTimeSplines, SINGLEPROCESS, "read GRACE data", Conversion, Grace, PotentialCoefficients)

/***********************************************/

void GraceCoefficients2BlockMeanTimeSplines::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName fileNameOut, fileNameCov, fileNameTimes, fileNameIntervals;
    std::vector<FileName> fileNameIn;

    readConfig(config, "outputfileTimeSplines",           fileNameOut,       Config::MUSTSET,  "", "");
    readConfig(config, "outputfileTimeSplinesCovariance", fileNameCov,       Config::OPTIONAL, "", "only the variances are saved");
    readConfig(config, "outputfileTimeSeries",            fileNameTimes,     Config::OPTIONAL, "", "mid points of non-empty intervals");
    readConfig(config, "outputfileTimeIntervals",         fileNameIntervals, Config::OPTIONAL, "", "monthly interval boundaries from first to last solution");
    readConfig(config, "inputfile",                       fileNameIn,        Config::MUSTSET,  "", "");
    if(isCreateSchema(config)) return;

    std::vector<SphericalHarmonics> harmonics;
    std::vector<Time> solutionStart, solutionEnd;
    Double GM = DEFAULT_GM;
    Double R  = DEFAULT_R;
    UInt   maxDegree = 0;
    for(auto &f : fileNameIn)
    {
      logStatus<<"read file <"<<f<<">"<<Log::endl;
      Time timeStart, timeEnd;
      harmonics.push_back(readPotentialCoefficients(f, timeStart, timeEnd));

      solutionStart.push_back(timeStart);
      solutionEnd.push_back(timeEnd);
      GM = harmonics.back().GM();
      R  = harmonics.back().R();
      maxDegree = std::max(maxDegree, harmonics.back().maxDegree());
    }
    const auto minTime = std::min_element(solutionStart.begin(), solutionStart.end());
    const auto maxTime = std::max_element(solutionEnd.begin(), solutionEnd.end());

    // monthly time intervals
    // ----------------------
    UInt year, month, day, dummy;
    Double second;
    minTime->date(year, month, day, dummy, dummy, second);

    std::vector<Time> intervals = {date2time(year, month, 1)};
    while(intervals.back() < *maxTime)
    {
      UInt year, month, day, dummy;
      Double second;
      intervals.back().date(year, month, day, dummy, dummy, second);
      intervals.push_back(date2time(year, month+1, 1));
    }

    std::vector<Matrix> cnmList, snmList, sigma2List;
    std::vector<Time> nonEmptyIntervals;
    for(UInt idxInterval=0; idxInterval<intervals.size()-1; idxInterval++)
    {
      Time epochMid = 0.5*(intervals.at(idxInterval)+intervals.at(idxInterval+1));

      std::vector<UInt> overlappingSolutions;
      for(UInt k=0; k<solutionStart.size(); k++)
        if(!((intervals.at(idxInterval) >= solutionEnd.at(k)) || (solutionStart.at(k) >= intervals.at(idxInterval+1))))
          overlappingSolutions.push_back(k);

      if(overlappingSolutions.size()>1)
        logWarning<<"Multiple solutions for interval ("<<intervals.at(idxInterval).dateStr()<<", "<<intervals.at(idxInterval+1).dateStr()<<"). Using nearest neighbour."<<Log::endl;

      SphericalHarmonics harm;
      if(overlappingSolutions.size() > 0)
      {
        auto it = std::min_element(overlappingSolutions.begin(), overlappingSolutions.end(), [&](UInt i, UInt j)
        {return std::abs((epochMid - (solutionStart.at(i)*0.5+solutionEnd.at(i)*0.5)).mjd()) < std::abs((epochMid - (solutionStart.at(j)*0.5+solutionEnd.at(j)*0.5)).mjd());});

        harm = harmonics.at(*it).get(maxDegree, 0, GM, R);
        nonEmptyIntervals.push_back(epochMid);
      }
      else
      {
        logWarning<<"Interval ("<<intervals.at(idxInterval).dateStr()<<", "<<intervals.at(idxInterval+1).dateStr()<<") is empty."<<Log::endl;
        harm = harm.get(maxDegree, 0, GM, R);
      }
      cnmList.push_back(harm.cnm());
      snmList.push_back(harm.snm());
      sigma2List.push_back(harm.sigma2x());
    }

    // write timeSplines file
    // ----------------------
    logStatus<<"write time splines to file <"<<fileNameOut<<">"<<Log::endl;
    writeFileTimeSplinesGravityfield(fileNameOut, GM, R, 0, intervals, cnmList, snmList);

    if(!fileNameCov.empty())
    {
      logStatus<<"write covariance time splines to file <"<<fileNameCov<<">"<<Log::endl;
      writeFileTimeSplinesCovariance(fileNameCov, GM, R, 0, maxDegree, 0, intervals, sigma2List);
    }
    if(!fileNameTimes.empty())
    {
      logStatus<<"write midpoints of non-empty intervals to <"<<fileNameTimes<<">"<<Log::endl;
      InstrumentFile::write(fileNameTimes, Arc(nonEmptyIntervals, Matrix(nonEmptyIntervals.size(), 1)));
    }
    if(!fileNameIntervals.empty())
    {
      logStatus<<"write monthly intervals to <"<<fileNameIntervals<<">"<<Log::endl;
      InstrumentFile::write(fileNameIntervals, Arc(intervals, Matrix(intervals.size(), 1)));
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

SphericalHarmonics GraceCoefficients2BlockMeanTimeSplines::readPotentialCoefficients(const FileName &name, Time &timeStart, Time &timeEnd)
{
  try
  {
    InFile file(name);

    Double  GM = 0.3986004415e15;
    Double  R  = 0.6378136460e07;
    Matrix  cnm, snm, sigma2cnm, sigma2snm;
    timeStart = date2time(9999,1,1);
    timeEnd   = date2time(1, 1, 1);

    Bool seekR = FALSE;
    Bool seekGM = FALSE;
    Bool seekDegree = FALSE;

    std::string line;
    while(std::getline(file, line))
    {
      if(line.find("dimension") != std::string::npos)
        seekDegree=TRUE;
      if( (line.find("degree") != std::string::npos) && seekDegree )
      {
        auto start_search = line.find(":");
        std::stringstream ss(line.substr(start_search+1));
        UInt degree = INFINITYDEGREE;
        ss>>degree;
        cnm       = Matrix(degree+1, Matrix::TRIANGULAR, Matrix::LOWER);
        snm       = Matrix(degree+1, Matrix::TRIANGULAR, Matrix::LOWER);
        sigma2cnm = Matrix(degree+1, Matrix::TRIANGULAR, Matrix::LOWER);
        sigma2snm = Matrix(degree+1, Matrix::TRIANGULAR, Matrix::LOWER);
        seekDegree = FALSE;
      }

      if(line.find("SHM ")==0)
      {
        UInt degree = String::toInt(line.substr(6, 5));
        cnm       = Matrix(degree+1, Matrix::TRIANGULAR, Matrix::LOWER);
        snm       = Matrix(degree+1, Matrix::TRIANGULAR, Matrix::LOWER);
        sigma2cnm = Matrix(degree+1, Matrix::TRIANGULAR, Matrix::LOWER);
        sigma2snm = Matrix(degree+1, Matrix::TRIANGULAR, Matrix::LOWER);
      }
      if(line.find("EARTH ")==0)
      {
        GM = String::toDouble(line.substr(6, 16));
        R  = String::toDouble(line.substr(23, 16));
      }

      if(line.find("mean_equator_radius") != std::string::npos)
        seekR = TRUE;
      if( (line.find("value") != std::string::npos) && seekR )
      {
        auto start_search = line.find(":");
        std::stringstream ss(line.substr(start_search+1));
        ss>>R;
        seekR = FALSE;
      }

      if(line.find("earth_gravity_param") != std::string::npos)
        seekGM = TRUE;
      if( (line.find("value") != std::string::npos) && seekGM )
      {
        auto start_search = line.find(":");
        std::stringstream ss(line.substr(start_search+1));
        ss>>GM;
        seekGM = FALSE;
      }

      if((line.find("GRCOEF")==0)||(line.find("GRCOF2")==0))
      {
        const UInt n    = String::toInt(line.substr(6, 5));
        const UInt m    = String::toInt(line.substr(11, 5));
        cnm(n,m)        = String::toDouble(line.substr(17, 18));
        snm(n,m)        = String::toDouble(line.substr(36, 18));
        sigma2cnm(n,m)  = String::toDouble(line.substr(55, 10));
        sigma2snm(n,m)  = String::toDouble(line.substr(66, 10));
        sigma2cnm(n,m) *= sigma2cnm(n,m);
        sigma2snm(n,m) *= sigma2snm(n,m);

        const UInt yearStart  = String::toInt(line.substr(77, 4));
        const UInt monthStart = String::toInt(line.substr(81, 2));
        const UInt dayStart   = String::toInt(line.substr(83, 2));
        timeStart = std::min(timeStart, date2time(yearStart, monthStart, dayStart));

        const UInt yearEnd    = String::toInt(line.substr(91, 4));
        const UInt monthEnd   = String::toInt(line.substr(95, 2));
        const UInt dayEnd     = String::toInt(line.substr(97, 2));
        timeEnd = std::max(timeEnd, date2time(yearEnd, monthEnd, dayEnd));
      }
    }

    return SphericalHarmonics(GM, R, cnm, snm, sigma2cnm, sigma2snm);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
