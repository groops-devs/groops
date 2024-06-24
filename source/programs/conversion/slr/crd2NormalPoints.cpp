/***********************************************/
/**
* @file crd2NormalPoints.cpp
*
* @brief read SLR data from CRD format.
*
* @author Torsten Mayer-Guerr
* @author Barbara Suesser-Rechberger
* @date 2022-05-07
*
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
Converts \href{https://ilrs.gsfc.nasa.gov/data_and_products/formats/crd.html}{CRD file}
and writes an \file{instrument file (METEOROLOGICAL)}{instrument} including meteorological data like
temperature, air pressure and humidity as well as an \file{instrument file (SATELLITELASERRANGING)}{instrument}
including normal point data like range, accuracy, redundancy, wavelength and window size.
)";

/***********************************************/

#include "programs/program.h"
#include "base/string.h"
#include "inputOutput/file.h"
#include "files/fileInstrument.h"

/***** CLASS ***********************************/

/** @brief Read SLR data from CRD format.
* @ingroup programsConversionGroup */
class Crd2NormalPoints
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(Crd2NormalPoints, SINGLEPROCESS, "read SLR data from CRD format for a given satellite", Conversion, Slr, Instrument)

/***********************************************/

void Crd2NormalPoints::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName fileNameNormalPoints;
    FileName fileNameMeteorological;
    FileName fileNameIn;

    readConfig(config, "outputfileNormalPoints",    fileNameNormalPoints,   Config::MUSTSET,  "output/normalsPoints_{station}.dat", "variable {station} available");
    readConfig(config, "outputfileMeteorological",  fileNameMeteorological, Config::MUSTSET,  "output/meteorological_{station}.dat", "variable {station} available");
    readConfig(config, "inputfileSlrData",          fileNameIn,            Config::MUSTSET,  "", "SLR CRD files");
    if(isCreateSchema(config)) return;

    // ==============================

    std::string stationMonumentNumber;
    Double transmitWavelength;
    std::string satelliteName;
    Time sessionDate;

    logStatus<<"read file <"<<fileNameIn<<">"<<Log::endl;
    InFile file(fileNameIn);
    std::string line;

    std::vector<SatelliteLaserRangingEpoch> normalPointEpoch;
    std::map<Time, MeteorologicalEpoch> metEpoch;
    std::multimap<std::string, std::map<Time, MeteorologicalEpoch>> metEpochStations;
    std::multimap<std::string, std::vector<SatelliteLaserRangingEpoch>> normalPointEpochStations;
    std::vector<std::string> stationMonumentNumbers;

    // Get a reference to MeteorologicalEpoch object at specified time. MeteorologicalEpoch object
    // will be default created (all NAN) if it does not already exist.
    auto metEpochGetOrInsert = [&metEpoch](const Time &time)-> MeteorologicalEpoch&
    {
      if(metEpoch.count(time) == 0)
      {
        metEpoch[time] = {};
        metEpoch[time].time = time;
        metEpoch[time].temperature = NAN;
        metEpoch[time].pressure = NAN;
        metEpoch[time].humidity = NAN;
        metEpoch[time].windSpeed = NAN;
        metEpoch[time].solarRadiation = NAN;
        metEpoch[time].precipitation = NAN;
      }
      return metEpoch[time];
    };

    while(std::getline(file, line))
    {
      if(line.empty())
        continue;
      std::string type = String::upperCase(line.substr(0,2));
      std::stringstream ss(line.substr(2));

      if(type == "H2") // station header
      {
        std::string stationName;
        ss>>stationName>>stationMonumentNumber;

        if(!stationMonumentNumbers.empty())
        {
          if(std::find(stationMonumentNumbers.begin(), stationMonumentNumbers.end(), stationMonumentNumber) == stationMonumentNumbers.end())
            stationMonumentNumbers.push_back(stationMonumentNumber);

        }
        else
        {
          stationMonumentNumbers.push_back(stationMonumentNumber);
        }
      }
      else if(type == "H3") // target header
      {
        ss>>satelliteName;
      }
      else if(type == "H4") // session header
      {
        Int dataType,year, month, day;

        ss>>dataType>>year>>month>>day;
        // Session time in UTC
        sessionDate = date2time(year, month, day);
      }
      else if(type == "H5") // Prediction Header
      {
      }
      else if(type == "H8") // end of session
      {
        normalPointEpochStations.insert({stationMonumentNumber, normalPointEpoch});
        metEpochStations.insert({stationMonumentNumber, metEpoch});

        // Clear the vector and map for the next session
        normalPointEpoch.clear();
        metEpoch.clear();
        continue;

      }
      else if(type == "C0")
      {
        Int detailType;
        ss>>detailType>>transmitWavelength;
      }
      else if(type == "11") // normal points record
      {
        Double      seconds, timeOfFlight;
        std::string systemId;
        Int         epoch, count;
        Double      window, binRms;
        SatelliteLaserRangingEpoch nPointEpoch;

        ss>>seconds>>timeOfFlight>>systemId>>epoch>>window>>count>>binRms;

        // Time in UTC
        Time time = sessionDate + seconds2time(seconds);

        switch(epoch)
        {
          case 0: time -= seconds2time(timeOfFlight);     break; // ground receive time (at SRP) (two-way)
          case 1: time -= seconds2time(0.5*timeOfFlight); break; // spacecraft bounce time (two-way)
          case 2: /* standard */                          break; // ground transmit time (at SRP) (two-way)
          default:
            logWarning<<epoch<<" Epoch Event (indicates the time event reference) not implemented"<<Log::endl;
        }

        // GROOPS uses GPS time
        nPointEpoch.time = timeUTC2GPS(time);
        nPointEpoch.range = 0.5*timeOfFlight*LIGHT_VELOCITY;

        if(count > 1)
          nPointEpoch.accuracy = 0.5*binRms*1e-12/std::sqrt(count-1)*LIGHT_VELOCITY;
        else
          nPointEpoch.accuracy = 0.5*binRms*1e-12*LIGHT_VELOCITY;

        nPointEpoch.redundancy = count;
        nPointEpoch.wavelength = transmitWavelength*1e-9;
        nPointEpoch.window = window;
        nPointEpoch.azmiuth = NAN;
        nPointEpoch.elevation = NAN;

        normalPointEpoch.push_back(nPointEpoch);
      }
      else if(type == "20") // meteorological record
      {
        Double seconds, surfacePressure, surfaceTemperature, relativeHumidity;

        ss>>seconds>>surfacePressure>>surfaceTemperature>>relativeHumidity;

        // GROOPS uses GPS time
        Time time = timeUTC2GPS(sessionDate + seconds2time(seconds));

        // Is neccessary due to record 20 and 21 may be in unspecified order.
        MeteorologicalEpoch& meteorologicalEpoch = metEpochGetOrInsert(time);
        // surface temperature in [K]
        meteorologicalEpoch.temperature = surfaceTemperature;
        // pressure in [Pa]: 1 millibar = 1 hPa = 100 Pa
        meteorologicalEpoch.pressure = surfacePressure*100;
        // relative humidity in percent
        meteorologicalEpoch.humidity = relativeHumidity;

      }
      else if(type == "21")
      {
        Double seconds, windSpeed;

        ss>>seconds>>windSpeed;

        // GROOPS uses GPS time
        Time time = timeUTC2GPS(sessionDate + seconds2time(seconds));

        // Is neccessary due to record 20 and 21 may be in unspecified order.
        MeteorologicalEpoch& meteorologicalEpoch = metEpochGetOrInsert(time);
        // Set wind speed at correct time
        meteorologicalEpoch.windSpeed = windSpeed;
      }
      else if(type == "30")
      {
        // Angles are not considered in the moment.
        logStatus<<"Record 30 available" <<Log::endl;
      }
      else
      {
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
          for(const auto &np : npStat.second)
            arcSlr.push_back(np);
      }

      for(const auto &metStat : metEpochStations)
      {
        if(metStat.first == stationMonumentNr)
          for(const auto &met: metStat.second)
            arcMeteorological.push_back(met.second);
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
