/***********************************************/
/**
* @file cstg2NormalPoints.cpp
*
* @brief read SLR data from CSTG format.
*
* @author Barbara Suesser-Rechberger
* @date 2022-12-11
*
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
Converts \href{https://ilrs.gsfc.nasa.gov/data_and_products/data/npt/npt_format.html}{CSTG file} provided by the \href{https://ilrs.gsfc.nasa.gov/}{ILRS}
and writes an \file{instrument file (METEOROLOGICAL)}{instrument} including meteorological data like
temperature, air pressure and humidity as well as an \file{instrument file (SATELLITELASERRANGING)}{instrument}
including normal point data like range, accuracy, redundancy, wavelength and window size.
)";

/***********************************************/

#include "programs/program.h"
#include "base/string.h"
#include "inputOutput/file.h"
#include "files/fileInstrument.h"
#include <string_view>

/***** CLASS ***********************************/

/** @brief Read SLR data from CSTG format.
* @ingroup programsConversionGroup */
class Cstg2NormalPoints
{
  std::vector<std::string> splitString(const std::string& str);

public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(Cstg2NormalPoints, SINGLEPROCESS, "read SLR data from CSTG format for a given satellite", Conversion, Slr, Instrument)

/***********************************************/

void Cstg2NormalPoints::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName fileNameNormalPoints;
    FileName fileNameMeteorological;
    FileName fileNameIn;

    readConfig(config, "outputfileNormalPoints",    fileNameNormalPoints,   Config::MUSTSET,  "output/normalsPoints_{station}.dat", "variable {station} available");
    readConfig(config, "outputfileMeteorological",  fileNameMeteorological, Config::MUSTSET,  "output/meteorological_{station}.dat", "variable {station} available");
    readConfig(config, "inputfileSlrData",          fileNameIn,             Config::MUSTSET,  "", "SLR CSTG file");
    if(isCreateSchema(config)) return;

    // ==============================

    std::string stationMonumentNumber;
    std::string satelliteName;
    Time sessionDate;

    logStatus<<"read file <"<<fileNameIn<<">"<<Log::endl;
    InFile file(fileNameIn);

    std::vector<SatelliteLaserRangingEpoch> normalPointEpoch;
    std::vector<MeteorologicalEpoch> metEpoch;
    std::multimap<std::string, std::vector<MeteorologicalEpoch>> metEpochStations;
    std::multimap<std::string, std::vector<SatelliteLaserRangingEpoch>> normalPointEpochStations;
    std::vector<std::string> stationMonumentNumbers;
    std::string fileLine;

    while(std::getline(file, fileLine))
    {
      // Check for session marker
      if (fileLine == "99999")
      {
        SatelliteLaserRangingEpoch nPointEpoch;
        MeteorologicalEpoch mEpoch;
        std::string headerLine;
        std::getline(file, headerLine);

        // Extract header data
        Int year = std::stoi(headerLine.substr(7,2));
        year += (year <= 50) ? 2000 : 1900;
        // Session time in UTC
        Time sessionDate = date2time(year,1,1)+mjd2time(std::stoi(headerLine.substr(9,3))-1);

        // Store different station monument numbers
        std::string stationMonumentNumber = headerLine.substr(12,4);
        if(!stationMonumentNumbers.empty())
        {
          if(std::find(stationMonumentNumbers.begin(), stationMonumentNumbers.end(), stationMonumentNumber) == stationMonumentNumbers.end())
            stationMonumentNumbers.push_back(stationMonumentNumber);
        }
        else
        {
          stationMonumentNumbers.push_back(stationMonumentNumber);
        }
        Double transmitWavelength = std::stod(headerLine.substr(20,4));
        Int normalPointWindowIndicator = std::stoi(headerLine.substr(42,1));

        // Extract all data records, until next session marker ("99999")
        std::vector<std::string> dataRecords;
        while(file.good())
        {
          // Remember current position, this is needed if the extracted line have to be put back into stream (undo effect of getline)
          auto preReadPos = file.tellg();
          std::getline(file,fileLine);

          if(fileLine.empty())
            continue;
          if(fileLine != "99999")
            dataRecords.push_back(fileLine);
          else
          {
            // Set read pointer back to session marker (next getline will again extract the session marker)
            file.seekg(preReadPos, file.beg);
            break;
          }
        }

        for(auto dataRecord : dataRecords)
        {
          // In cstg file the seconds are given in the units of 0.1 us, GROOPS uses GPS time
          nPointEpoch.time = timeUTC2GPS(sessionDate + seconds2time(std::stod(dataRecord.substr(0, 12))*1e-7));
          nPointEpoch.range = 0.5*std::stod(dataRecord.substr(12, 12))*1e-12*LIGHT_VELOCITY;

          Double binRms = std::stod(dataRecord.substr(24, 7));
          Int count = std::stoi(dataRecord.substr(43, 4));
          if(count > 1)
            nPointEpoch.accuracy = 0.5*binRms*1e-12/std::sqrt(count-1)*LIGHT_VELOCITY;
          else
            nPointEpoch.accuracy = 0.5*binRms*1e-12*LIGHT_VELOCITY;

          nPointEpoch.redundancy = count;

          /* ILRS: The user of the data should interpret the value given as follows:
            3000 - 9999: units are 0.1 nanometer
            1000 - 2999: units are 1.0 nanometer
          */
          if(transmitWavelength >= 3000)
            nPointEpoch.wavelength = transmitWavelength*1e-10;
          else
            nPointEpoch.wavelength = transmitWavelength*1e-9;

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

          nPointEpoch.azmiuth = NAN;
          nPointEpoch.elevation = NAN;

          mEpoch.time = nPointEpoch.time;
          // surface temperature in [K]: in cstg file it is given in 0.1 degree Kelvin
          mEpoch.temperature = std::stod(dataRecord.substr(36, 4))*1e-1;
          // pressure in [Pa]: in cstg file it is given in 0.1 millibar = 10 Pa
          mEpoch.pressure = std::stod(dataRecord.substr(31, 5))*10;
          // relative humidity in percent
          mEpoch.humidity = std::stod(dataRecord.substr(40, 3));

          mEpoch.windSpeed = NAN;
          mEpoch.solarRadiation = NAN;
          mEpoch.precipitation = NAN;

          normalPointEpoch.push_back(nPointEpoch);
          metEpoch.push_back(mEpoch);

        }

        // End of session
        normalPointEpochStations.insert({stationMonumentNumber, normalPointEpoch});
        metEpochStations.insert({stationMonumentNumber, metEpoch});
        normalPointEpoch.clear();
        metEpoch.clear();
      }
    } //while

    // write results
    // -------------
    for(auto stationMonumentNr : stationMonumentNumbers)
    {
      SatelliteLaserRangingArc arcSlr;
      MeteorologicalArc arcMeteorological;

      for(const auto &npStat : normalPointEpochStations)
      {
        if(npStat.first == stationMonumentNr)
          for(const auto &np : npStat.second)
            arcSlr.push_back(np);
      }

      for(const auto &metStat : metEpochStations)
      {
        if(metStat.first == stationMonumentNr)
          for(const auto &met: metStat.second)
            arcMeteorological.push_back(met);
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

std::vector<std::string> Cstg2NormalPoints::splitString(const std::string& str)
{
  std::vector<std::string> recordLines;

  std::stringstream ss(str);
  std::string line;
  while(std::getline(ss, line, '\n'))
  {
    recordLines.push_back(line);
  }
  return recordLines;
}

/***********************************************/
