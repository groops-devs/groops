/***********************************************/
/**
* @file sinex2StationPostSeismicDeformation.cpp
*
* @brief Convert ITRF post-seismic deformation SINEX file to Instrument VECTOR3D file.
*
* @author Sebastian Strasser
* @date 2018-05-23
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
Convert ITRF post-seismic deformation
\href{http://www.iers.org/IERS/EN/Organization/AnalysisCoordinator/SinexFormat/sinex.html}{SINEX file}
to \configFile{outputfileInstrument}{instrument} (VECTOR3D).

See also \program{Sinex2StationPosition} and \program{Sinex2StationDiscontinuities}.
)";

/***********************************************/

#include "programs/program.h"
#include "base/string.h"
#include "classes/timeSeries/timeSeries.h"
#include "files/fileInstrument.h"
#include "inputOutput/fileSinex.h"

/***** CLASS ***********************************/

/** @brief Convert ITRF post-seismic deformation SINEX file to Instrument VECTOR3D file.
* @ingroup programsConversionGroup */
class Sinex2StationPostSeismicDeformation
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(Sinex2StationPostSeismicDeformation, SINGLEPROCESS, "Convert ITRF post-seismic deformation SINEX file to Instrument VECTOR3d file.", Conversion, Gnss, Instrument)

/***********************************************/

void Sinex2StationPostSeismicDeformation::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName fileNameInstrument, fileNameSinex;
    TimeSeriesPtr timeSeries;
    std::string stationName;
    Bool localLevelFrame;

    readConfig(config, "outputfileInstrument", fileNameInstrument, Config::MUSTSET,  "",  "deformation time series");
    readConfig(config, "inputfileSinex",       fileNameSinex,      Config::MUSTSET,  "",  "ITRF post-seismic deformation SINEX file");
    readConfig(config, "timeSeries",           timeSeries,         Config::MUSTSET,  "",  "compute deformation for these epochs");
    readConfig(config, "stationName",          stationName,        Config::MUSTSET,  "",  "");
    readConfig(config, "localLevelFrame",      localLevelFrame,    Config::DEFAULT,  "0", "output in North, East, Up local-level frame");
    if(isCreateSchema(config)) return;

    std::vector<Time> times = timeSeries->times();
    Matrix A(times.size(), 4);
    const std::vector<Char> axes = {'N', 'E', 'H'};
    std::transform(stationName.begin(), stationName.end(), stationName.begin(), ::toupper);

    logStatus << "read SINEX file <" << fileNameSinex << ">" << Log::endl;
    Sinex sinex(fileNameSinex);
    std::vector<Sinex::Parameter> parameters = sinex.getBlock<Sinex::SinexSolutionVector>("SOLUTION/ESTIMATE")->parameters();
    for(UInt i = 0; i < parameters.size(); i+=2)
    {
      if(parameters.at(i).siteCode != stationName)
        continue;

      if(parameters.at(i).parameterType.front() != 'A' || parameters.at(i+1).parameterType.front() != 'T')
        throw(Exception(i%"incorrect parameter pair at index %i: "s + parameters.at(i).parameterType + " and " + parameters.at(i+1).parameterType));

      const Double amplitude      = parameters.at(i).value;
      const Double relaxationTime = parameters.at(i+1).value;
      const Double referenceTime  = parameters.at(i).time.decimalYear();

      const UInt col = 1 + std::distance(axes.begin(), std::find(axes.begin(), axes.end(), parameters.at(i).parameterType.back()));
      for(UInt idEpoch = 0; idEpoch < times.size(); idEpoch++)
      {
        if(times.at(idEpoch).decimalYear() < referenceTime)
          continue;

        if(parameters.at(i).parameterType.substr(1,3) == "EXP")
          A(idEpoch, col) += amplitude * (1. - std::exp(-(times.at(idEpoch).decimalYear() - referenceTime) / relaxationTime));
        else if(parameters.at(i).parameterType.substr(1,3) == "LOG")
          A(idEpoch, col) += amplitude * std::log(1. + (times.at(idEpoch).decimalYear() - referenceTime) / relaxationTime);
        else
          throw(Exception(i%"unknown parameter type at index %i: "s + parameters.at(i).parameterType));
      }
    }

    if(!localLevelFrame)
    {
      std::vector<std::string> stationInfos = sinex.getBlock<Sinex::SinexText>("SITE/ID")->lines();
      auto iter = std::find_if(stationInfos.begin(), stationInfos.end(), [&](const std::string &s){ return s.substr(1,4) == stationName; });
      if(iter != stationInfos.end())
      {
        const Double longitude = String::toDouble(iter->substr(44, 3)) + String::toDouble(iter->substr(48, 2))/60 + String::toDouble(iter->substr(51, 4))/3600;
        const Double latitude  = String::toDouble(iter->substr(56, 3)) + String::toDouble(iter->substr(60, 2))/60 + String::toDouble(iter->substr(63, 4))/3600;
        const Double height    = String::toDouble(iter->substr(68, 7));

        Ellipsoid ellipsoid;
        const Transform3d lnof2trf = localNorthEastUp(ellipsoid(Angle(longitude*DEG2RAD), Angle(latitude*DEG2RAD), height));
        copy(A.column(1,3) * lnof2trf.matrix().trans(), A.column(1,3));
      }
    }

    logStatus << "write Instrument file <" << fileNameInstrument << ">" << Log::endl;
    InstrumentFile::write(fileNameInstrument, Arc(times, A, Epoch::VECTOR3D));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
