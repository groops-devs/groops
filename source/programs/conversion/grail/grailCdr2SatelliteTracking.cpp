/***********************************************/
/**
* @file grailCdr2SatelliteTracking.cpp
*
* @brief read CDR GRAIL data.
*
* @author Beate Klinger
* @date 2013-01-22
*
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program converts low-low satellite data measured by the K-band ranging system
from the GRAIL format into \file{instrument file (SATELLITETRACKING)}{instrument}.
The \config{inputfile}s contain also corrections for antenna offsets
and the so called light time correction.
The corrections can be stored in additional files in the same format as the observations.
If a phase break is found an artificial gap is created.
)";

/***********************************************/

#include "programs/program.h"
#include "base/string.h"
#include "files/fileInstrument.h"

/***** CLASS ***********************************/

/** @brief read CDR GRAIL data.
* @ingroup programsConversionGroup */
class GrailCdr2SatelliteTracking
{
  Double timeBias;

  void readFile(const FileName fileName,
                SatelliteTrackingArc &arc,
                SatelliteTrackingArc &arcAntenna,
                SatelliteTrackingArc &arcLight,
                SatelliteTrackingArc &arcTemperature);

public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(GrailCdr2SatelliteTracking, SINGLEPROCESS, "read CDR GRAIL data", Conversion, Instrument)

/***********************************************/

void GrailCdr2SatelliteTracking::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName outName, antCentrName, lighttimeName, temperatureName;
    std::vector<FileName> fileNames;

    readConfig(config, "outputfileSatelliteTracking", outName,         Config::OPTIONAL, "", "");
    readConfig(config, "outputfileAntCentr",          antCentrName,    Config::OPTIONAL, "", "");
    readConfig(config, "outputfileLighttime",         lighttimeName,   Config::OPTIONAL, "", "");
    readConfig(config, "outputfileTemperature",       temperatureName, Config::OPTIONAL, "", "");
    readConfig(config, "approximateTimeBias",         timeBias,        Config::DEFAULT,  "0.0", "[seconds]");
    readConfig(config, "inputfile",                   fileNames,       Config::MUSTSET,  "", "");
    if(isCreateSchema(config)) return;

    logStatus<<"read input files"<<Log::endl;
    SatelliteTrackingArc arcSST, arcAntenna, arcLight, arcTemperature;
    for(auto &fileName : fileNames)
    {
      logStatus<<"read file <"<<fileName<<">"<<Log::endl;
      readFile(fileName, arcSST, arcAntenna, arcLight, arcTemperature);
    }

    Arc::printStatistics(arcSST);

    // Daten speichern
    // ---------------
    if(!outName.empty())
    {
      logStatus<<"write data to <"<<outName<<">"<<Log::endl;
      InstrumentFile::write(outName, arcSST);
    }
    if(!antCentrName.empty())
    {
      logStatus<<"write antenna center correction to <"<<antCentrName<<">"<<Log::endl;
      InstrumentFile::write(antCentrName, arcAntenna);
    }
    if(!lighttimeName.empty())
    {
      logStatus<<"write light time correction to <"<<lighttimeName<<">"<<Log::endl;
      InstrumentFile::write(lighttimeName, arcLight);
    }
    if(!temperatureName.empty())
    {
      logStatus<<"write temperature correction to <"<<temperatureName<<">"<<Log::endl;
      InstrumentFile::write(temperatureName, arcTemperature);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GrailCdr2SatelliteTracking::readFile(const FileName fileName, SatelliteTrackingArc &arc, SatelliteTrackingArc &arcAntenna, SatelliteTrackingArc &arcLight, SatelliteTrackingArc &arcTemperature)
{
  try
  {
    SatelliteTrackingEpoch epoch, epochAntenna, epochLight, epochTemperature;

    InFile file(fileName);
    if(!file.good())
    {
      logWarning<<"cannot open file: "<<fileName.str()<<", continue..."<<Log::endl;
      return;
    }
    file.exceptions(std::ios::badbit|std::ios::failbit);

    // Header einlesen
    UInt numberOfRecords = 0;
    std::string line;
    for(;;)
    {
      getline(file, line);
      if(line.find("NUMBER OF DATA RECORDS")==0)
        numberOfRecords = String::toInt(line.substr(31, 10));
      // Header fertig ?
      if(line.find("END OF HEADER")==0)
        break;
    }

    // Eigentliche Daten einlesen
    for(UInt i=0; i<numberOfRecords; i++)
    {
      std::string line;
      try
      {
        getline(file, line);
      }
      catch(std::exception &/*e*/)
      {
        break;
      }

      std::stringstream ss(line);

      // time (TDB) >> seconds past 12:00:00 noon 01-Jan-2000
      Int32 seconds;
      ss>>seconds;
      epoch.time            = timeTT2GPS(mjd2time(51544.5)+seconds2time(seconds+timeBias));
      epochLight.time       = epoch.time;
      epochAntenna.time     = epoch.time;
      epochTemperature.time = epoch.time;

      // range, range rate, range acceleration
      ss>>epoch.range>>epoch.rangeRate>>epoch.rangeAcceleration;

      // ionospheric range correction
      Double iono;
      ss>>iono;

      // time of flight correction
      ss>>epochLight.range>>epochLight.rangeRate>>epochLight.rangeAcceleration;

      // Ka-band antenna offset correction
      ss>>epochAntenna.range>>epochAntenna.rangeRate>>epochAntenna.rangeAcceleration;

      Double dummy1, dummy2, dummy3, dummy4;
      ss>>dummy1>>dummy2>>dummy3>>dummy4;

      // data quality flags
      Byte flag0, flag1, flag2, flag3, flag4, flag5, flag6, flag7;
      ss>>flag0>>flag1>>flag2>>flag3>>flag4>>flag5>>flag6>>flag7;

      // Ka-band temperature correction
      ss>>dummy1>>epochTemperature.range>>epochTemperature.rangeRate>>epochTemperature.rangeAcceleration;

      Bool phasebreak = ((flag7 & 1) == 1);
      if(phasebreak)
      logWarning<<epoch.time.dateTimeStr()<<": flag="<<Int(flag7)<<Log::endl;

      // Sprünge in den Daten suchen
      if(arc.size() && (epoch.time <= arc.at(arc.size()-1).time))
          throw(Exception("epoch.time <= arc.at(arc.size()-1).time"));
      if(phasebreak || (arc.size() && ((epoch.time-arc.at(arc.size()-1).time).seconds()<=5.1) && (fabs(epoch.range - arc.at(arc.size()-1).range)>80)) ||
        (arc.size() && ((epoch.time-arc.at(arc.size()-1).time).seconds()<=5.1) && (fabs(epoch.rangeRate - arc.at(arc.size()-1).rangeRate)>0.5)) ||
        (arc.size() && ((epoch.time-arc.at(arc.size()-1).time).seconds()<=5.1) && (fabs(epoch.rangeAcceleration - arc.at(arc.size()-1).rangeAcceleration)>0.0009)))
      {
        logWarning<<epoch.time.dateTimeStr()<<": phase break? (epoch.range - arc.at(arc.size()-1).range) = "<<epoch.range - arc.at(arc.size()-1).range<<Log::endl;
        arc.remove(arc.size()-1);
        arcAntenna.remove(arcAntenna.size()-1);
        arcLight.remove(arcLight.size()-1);
        arcTemperature.remove(arcTemperature.size()-1);
        continue;
      }
      else
      {
        arc.push_back(epoch);
        arcAntenna.push_back(epochAntenna);
        arcLight.push_back(epochLight);
        arcTemperature.push_back(epochTemperature);
      }
    }
  }

  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
