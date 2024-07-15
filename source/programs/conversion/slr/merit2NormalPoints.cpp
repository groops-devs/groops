/***********************************************/
/**
* @file merit2NormalPoints.cpp
*
* @brief Read SLR data (normal points) from merit format.
*
* @author Barbara Suesser-Rechberger
* @date 2023-01-17
*
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
Converts \href{https://ilrs.gsfc.nasa.gov/data_and_products/data/npt/npt_format.html}{MERIT II file}
and writes an \file{instrument file (METEOROLOGICAL)}{instrument} including meteorological data like
temperature, air pressure and humidity as well as an \file{instrument file (SATELLITELASERRANGING)}{instrument}
including normal point data like range, accuracy, redundancy, wavelength, window size, azimuth and elevation.
)";

/***********************************************/

#include "programs/program.h"
#include "base/string.h"
#include "inputOutput/file.h"
#include "files/fileInstrument.h"
#include <string_view>

/***** CLASS ***********************************/

/** @brief Read SLR data (normal points) from merit format.
* @ingroup programsConversionGroup */
class Merit2NormalPoints
{
  std::vector<std::string> splitString(const std::string& str);

public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(Merit2NormalPoints, SINGLEPROCESS, "Read SLR data (normal points) from merit format.", Conversion, Slr, Instrument)

/***********************************************/

void Merit2NormalPoints::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName fileNameNormalPoints;
    FileName fileNameMeteorological;
    FileName fileNameIn;

    readConfig(config, "outputfileNormalPoints",    fileNameNormalPoints,   Config::MUSTSET,  "output/normalsPoints_{station}.dat", "variable {station} available");
    readConfig(config, "outputfileMeteorological",  fileNameMeteorological, Config::MUSTSET,  "output/meteorological_{station}.dat", "variable {station} available");
    readConfig(config, "inputfileSlrData",          fileNameIn,             Config::MUSTSET,  "", "SLR MERIT II file");
    if(isCreateSchema(config)) return;

    logStatus<<"read file <"<<fileNameIn<<">"<<Log::endl;
    InFile file(fileNameIn);
    std::string line;

    std::string stationMonumentNumber;
    std::multimap<std::string, MeteorologicalEpoch> metEpochStations;
    std::multimap<std::string, SatelliteLaserRangingEpoch> normalPointEpochStations;
    std::vector<std::string> stationMonumentNumbers;

    while(std::getline(file, line))
    {
      if(line.empty())
        continue;
      SatelliteLaserRangingEpoch nPointEpoch;
      MeteorologicalEpoch mEpoch;

      Int normalPointWindowIndicator = String::toInt(line.substr(114, 1));
      //logStatus<<"normalPointWindowIndicator: <"<< normalPointWindowIndicator<<">"<<Log::endl;


      if(normalPointWindowIndicator > 0)
      {
        // Store different station monument numbers
        std::string stationMonumentNumber = line.substr(24, 4);
        if(!stationMonumentNumbers.empty())
        {
          if(std::find(stationMonumentNumbers.begin(), stationMonumentNumbers.end(), stationMonumentNumber) == stationMonumentNumbers.end())
            stationMonumentNumbers.push_back(stationMonumentNumber);
        }
        else
        {
          stationMonumentNumbers.push_back(stationMonumentNumber);
        }

        Int year = String::toInt(line.substr(7,2));
        //logStatus<<"Year: <"<< year<<">"<<Log::endl;
        year += (year <= 50) ? 2000 : 1900;

        // Session time in UTC
        Time sessionDate = date2time(year,1,1)+mjd2time(String::toInt(line.substr(9,3))-1);

        // In merit II file the seconds are given in the units of 0.1 us
        Time time = sessionDate + seconds2time(std::stod(line.substr(12,12))*1e-7);
        Int epoch = String::toInt(line.substr(119,1));
        Double timeOfFlight = std::stod(line.substr(45, 12))*1e-12;

        switch(epoch)
        {
          case 0: time -= seconds2time(timeOfFlight);     break; // ground receive time (at SRP) (two-way)
          case 1: time -= seconds2time(0.5*timeOfFlight); break; // spacecraft bounce time (two-way)
          case 2: /* standard */                          break; // ground transmit time (at SRP) (two-way)
          default:
            logWarning<<epoch<<" Epoch Event (indicates the time event reference) not implemented"<<Log::endl;
        }

        nPointEpoch.time = timeUTC2GPS(time);
        nPointEpoch.range = 0.5*timeOfFlight*LIGHT_VELOCITY;

        /* ILRS: The user of the data should interpret the value given as follows:
          3000 - 9999: units are 0.1 nanometer
          1000 - 2999: units are 1.0 nanometer
        */
        Double transmitWavelength = std::stod(line.substr(64, 4));
        if(transmitWavelength >= 3000)
          nPointEpoch.wavelength = transmitWavelength*1e-10;
        else
          nPointEpoch.wavelength = transmitWavelength*1e-9;

        // In the merit format there is no information about the bin RMS, we take the pass RMS
        Double passRMS = std::stod(line.substr(58, 7));
        nPointEpoch.redundancy = std::stod(line.substr(115, 4));
        logStatus<<"epoch: <"<< epoch<<">"<<Log::endl;

        if(nPointEpoch.redundancy > 1)
          nPointEpoch.accuracy = 0.5*passRMS*1e-12/std::sqrt(nPointEpoch.redundancy-1)*LIGHT_VELOCITY;
        else
          nPointEpoch.accuracy = 0.5*passRMS*1e-12*LIGHT_VELOCITY;


        if(normalPointWindowIndicator == 1)
          nPointEpoch.window = 5;
        else if(normalPointWindowIndicator == 3)
          nPointEpoch.window = 15;
        else if(normalPointWindowIndicator == 4)
          nPointEpoch.window = 20;
        else if(normalPointWindowIndicator == 5)
          nPointEpoch.window = 30;
        else if(normalPointWindowIndicator == 6)
          nPointEpoch.window = 60;
        else if(normalPointWindowIndicator == 7)
          nPointEpoch.window = 120;
        else if(normalPointWindowIndicator == 8)
          nPointEpoch.window = 180;
        else if(normalPointWindowIndicator == 9)
          nPointEpoch.window = 300;
        else
          nPointEpoch.window = NAN;

        // Azimuth and Elevation given in 0.1 mdegrees
        nPointEpoch.azmiuth = DEG2RAD*std::stod(line.substr(32, 7))*1e-4;
        nPointEpoch.elevation = DEG2RAD*std::stod(line.substr(39, 6))*1e-4;

        mEpoch.time = nPointEpoch.time;
        // surface temperature in [K]: in cstg file it is given in 0.1 degree Kelvin
        mEpoch.temperature = std::stod(line.substr(73, 4))*1e-1;
        // pressure in [Pa]: in cstg file it is given in 0.1 millibar = 10 Pa
        mEpoch.pressure = std::stod(line.substr(68, 5))*10;
        // relative humidity in percent
        mEpoch.humidity = std::stod(line.substr(77, 3));

        mEpoch.windSpeed = NAN;
        mEpoch.solarRadiation = NAN;
        mEpoch.precipitation = NAN;

        normalPointEpochStations.insert({stationMonumentNumber, nPointEpoch});
        metEpochStations.insert({stationMonumentNumber, mEpoch});

      } else
      {
        logStatus<<"File includes full-rate data. Please use merit2FullRate conversion program!<"<<fileNameNormalPoints<<">"<<Log::endl;
        break;
      }
    } // while


    // write results
    // -------------
    for(auto stationMonumentNr : stationMonumentNumbers)
    {
      SatelliteLaserRangingArc arcSlr;
      MeteorologicalArc arcMeteorological;

      for(const auto &npStat : normalPointEpochStations)
      {
        if(npStat.first == stationMonumentNr)
          arcSlr.push_back(npStat.second);
      }

      for(const auto &metStat : metEpochStations)
      {
        if(metStat.first == stationMonumentNr)
          arcMeteorological.push_back(metStat.second);
      }

      VariableList varList;
      varList.setVariable("station", stationMonumentNr);

      if(arcSlr.size())
      {
        logStatus<<"write normal point data to file <"<<fileNameNormalPoints(varList)<<">"<<Log::endl;
        InstrumentFile::write(fileNameNormalPoints(varList), arcSlr);
        Arc::printStatistics(arcSlr);
      }

      if(arcMeteorological.size())
      {
        logStatus<<"write meteorological data to file <"<<fileNameMeteorological(varList)<<">"<<Log::endl;
        InstrumentFile::write(fileNameMeteorological(varList), arcMeteorological);
        Arc::printStatistics(arcMeteorological);
      }
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
