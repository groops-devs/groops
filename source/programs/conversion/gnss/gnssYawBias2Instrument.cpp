/***********************************************/
/**
* @file gnssYawBias2Instrument.cpp
*
* @brief Convert yaw bias table from JPL to an InstrumentMisc file for a GPS PRN.
*
* Yaw bias file source: ftp://sideshow.jpl.nasa.gov/pub/gipsy_files/gipsy_params/yaw_bias_table.gz
*
* Extracts yaw bias from JPL yaw bias file based on PRN-SVN association in transmitter info file
* and writes an InstrumentMiscValue file for the PRN containing the biases.
*
* @author Sebastian Strasser
* @date 2016-11-14
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
Convert yaw bias table from JPL to instrument files per GPS PRN.

Yaw bias file source: \url{https://sideshow.jpl.nasa.gov/pub/gipsy_files/gipsy_params//yaw_bias_table.gz}

Extracts yaw bias from \config{inputfileYawBias} based on PRN-SVN association in \configFile{inputfiletransmitterInfo}{gnssStationInfo}
and writes an \configFile{outputfileInstrument}{instrument} (MISCVALUE) for the PRN containing the biases.

SVN 23 -3.5 deg yaw bias is changed to -0.5 deg.
Block IIF yaw bias is changed from -0.5 deg to -0.7 deg according to Kouba 2017 eclipse update manuscript.
)";

/***********************************************/

#include "programs/program.h"
#include "base/string.h"
#include "inputOutput/file.h"
#include "files/fileGnssStationInfo.h"
#include "files/fileInstrument.h"

/***** CLASS ***********************************/

/** @brief Convert yaw bias table from JPL to an InstrumentMisc file for a GPS PRN.
* @ingroup programsConversionGroup */
class GnssYawBias2Instrument
{
  class YawBias
  {
  public:
    UInt                svn;
    std::vector<Time>   timeStart;
    std::vector<Double> bias;
  };

  std::vector<YawBias> readYawBiasFile(FileName yawBiasFile) const;

public:
  void run(Config &config);
};

GROOPS_REGISTER_PROGRAM(GnssYawBias2Instrument, SINGLEPROCESS, "Convert yaw bias table from JPL to an InstrumentMisc file for a GPS PRN.", Conversion, Gnss)

/***********************************************/

void GnssYawBias2Instrument::run(Config &config)
{
  try
  {
    FileName fileNameInstrument, fileNameYawBias, fileNameTransmitterInfo;

    readConfig(config, "outputfileInstrument",     fileNameInstrument,      Config::MUSTSET,  "", "yaw bias instrument file");
    readConfig(config, "inputfileYawBias",         fileNameYawBias,         Config::MUSTSET,  "{groopsDataDir}/gnss/transmitterGPS/yawBias/config/yaw_bias_table", "JPL yaw bias table");
    readConfig(config, "inputfileTransmitterInfo", fileNameTransmitterInfo, Config::MUSTSET,  "{groopsDataDir}/gnss/transmitterGPS/transmitterInfo/igs/igs14/transmitterInfo_igs14.{prn}.xml", "GNSS transmitter info file");
    if(isCreateSchema(config)) return;

    logStatus << "read yaw bias file <" << fileNameYawBias << ">" << Log::endl;
    std::vector<YawBias> biasList = readYawBiasFile(fileNameYawBias);

    logStatus << "read transmitter info file <" << fileNameTransmitterInfo << ">" << Log::endl;
    GnssStationInfo transmitterInfo;
    readFileGnssStationInfo(fileNameTransmitterInfo, transmitterInfo);
    if(transmitterInfo.markerName != "GPS")
      throw(Exception("transmitter type not implemented: " + transmitterInfo.markerName));

    logStatus << "write instrument file <" << fileNameInstrument << ">" << Log::endl;

    MiscValueArc arc;
    MiscValueEpoch epoch; // initially no yaw bias
    epoch.time  = Time();
    epoch.value = 0;
    arc.push_back(epoch);
    for(UInt idAnt = 0; idAnt < transmitterInfo.antenna.size(); idAnt++)
    {
      // initial bias if PRN assignment has changed to a new satellite
      MiscValueEpoch epoch;
      epoch.time  = transmitterInfo.antenna.at(idAnt).timeStart;
      epoch.value = 0;
      arc.push_back(epoch);

      // check if bias data exist for this satellite
      UInt svn = std::stoi(transmitterInfo.antenna.at(idAnt).serial.substr(1));
      UInt idBias = NULLINDEX;
      for(UInt i = 0; i < biasList.size(); i++)
        if(biasList.at(i).svn == svn)
          idBias = i;

      if(idBias != NULLINDEX)
      {
        // find first bias relevant for time period
        UInt idx = 0;
        for(UInt i = 1; i < biasList.at(idBias).timeStart.size(); i++)
          if(biasList.at(idBias).timeStart.at(i) < transmitterInfo.antenna.at(idAnt).timeStart)
            idx++;
        arc.at(arc.size()-1).value = biasList.at(idBias).bias.at(idx);

        // further bias changes for the same satellite
        for(UInt i = idx+1; i < biasList.at(idBias).timeStart.size(); i++)
          if(biasList.at(idBias).timeStart.at(i) < transmitterInfo.antenna.at(idAnt).timeEnd)
          {
            MiscValueEpoch epoch;
            epoch.time  = biasList.at(idBias).timeStart.at(i);
            epoch.value = biasList.at(idBias).bias.at(i);
            arc.push_back(epoch);
          }
      }
    }

    InstrumentFile::write(fileNameInstrument, arc);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

std::vector<GnssYawBias2Instrument::YawBias> GnssYawBias2Instrument::readYawBiasFile(FileName yawBiasFile) const
{
  try
  {
    InFile file(yawBiasFile);

    std::vector<YawBias> biasList;
    Time timeJ2000 = mjd2time(J2000);

    std::string line;
    while(std::getline(file, line))
    {
      std::string lineID = line.substr(0,3);
      if(lineID != "GPS")
        continue;

      std::vector<std::string> strings = String::split(line, "\t ");

      YawBias bias;
      bias.svn       = std::stoi(strings.at(0).substr(3));
      UInt count     = std::stoi(strings.at(1));

      for(UInt i = 0; i < count; i++)
      {
        Time timeStart = timeJ2000 + seconds2time(std::stod(strings.at(2+2*i)));
        bias.timeStart.push_back(timeStart < Time() ? Time() : timeStart);

        Double biasValue = std::stod(strings.at(2+2*i+1));
        if(biasValue == 0)                          // no bias
          bias.bias.push_back(0);
        else if(bias.svn == 23 && biasValue == -3.5) // SVN 23 -3.5 deg ==> -0.5 deg
          bias.bias.push_back(-0.5);
        else if(bias.svn >= 62 && bias.svn <= 73 && biasValue == -2) // Block IIF -0.5 deg ==> -0.7 deg (see Kouba 2017 eclipse update)
          bias.bias.push_back(-0.7);
        else if(biasValue ==  1 || biasValue ==  2) // normal bias
          bias.bias.push_back(0.5);
        else if(biasValue == -1 || biasValue == -2) // anti-normal bias
          bias.bias.push_back(-0.5);
        else                                        // actual bias value
          bias.bias.push_back(biasValue);
      }

      biasList.push_back(bias);
    }

    return biasList;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
