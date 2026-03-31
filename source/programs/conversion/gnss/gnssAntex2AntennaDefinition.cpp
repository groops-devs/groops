/***********************************************/
/**
* @file gnssAntex2AntennaDefinition.cpp
*
* @brief Converts IGS ANTEX file to GNSS metadata and antenna definition files.
*
* @author Torsten Mayer-Guerr
* @author Sebastian Strasser
* @date 2012-11-23
*
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
Converts metadata and antenna definitions from the \href{https://files.igs.org/pub/data/format/antex14.txt}{IGS ANTEX format}.
to \configFile{antennaDefinition}{gnssAntennaDefinition}, \configFile{transmitterInfo}{platform}, and
\configFile{transmitterList}{stringList} files for the respective GNSS and for the list of ground station antennas.

The \file{transmitterInfo}{platform} files for GLONASS satellites should then be updated using \program{GnssGlonassFrequencyNumberUpdate}.
)";

/***********************************************/

#include "programs/program.h"
#include "base/string.h"
#include "inputOutput/file.h"
#include "files/filePlatform.h"
#include "files/fileStringTable.h"

/***** CLASS ***********************************/

/** @brief Converts IGS ANTEX file to GNSS metadata and antenna definition files.
* @ingroup programsConversionGroup */
class GnssAntex2AntennaDefinition
{
  Bool getLine(InFile &file, std::string &line, std::string &label, Bool throwException=FALSE) const;
  Bool testLabel(const std::string &labelInLine, const std::string &label, Bool optional=FALSE) const;
  static Bool equal(GnssAntennaDefinitionPtr antenna, GnssAntennaDefinitionPtr antenna2);
  void addAntenna(GnssAntennaDefinitionPtr antenna, std::vector<GnssAntennaDefinitionPtr> &antennaList);
  void addTransmitter(const PlatformGnssAntenna &platformAntenna, const std::string &markerName, const std::string &markerNumber, std::vector<Platform> &transmitterList);

public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(GnssAntex2AntennaDefinition, SINGLEPROCESS, "Converts IGS ANTEX file to GNSS metadata and antenna definition files.", Conversion, Gnss)

/***********************************************/

void GnssAntex2AntennaDefinition::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName fileNameAntennaStation, fileNameAntennaTransmitter, fileNameTransmitterInfo;
    FileName fileNameTransmitterListGps, fileNameTransmitterListGlonass, fileNameTransmitterListGalileo;
    FileName fileNameTransmitterListBeiDou, fileNameTransmitterListQzss, fileNameTransmitterListIrnss;
    FileName fileNameAntex;
    Time     timeStart;
    Bool     setZero;

    readConfig(config, "outputfileAntennaDefinitionStation",     fileNameAntennaStation,         Config::OPTIONAL, "",  "antenna center variations");
    readConfig(config, "outputfileAntennaDefinitionTransmitter", fileNameAntennaTransmitter,     Config::OPTIONAL, "",  "antenna center variations");
    readConfig(config, "outputfileTransmitterInfo",              fileNameTransmitterInfo,        Config::OPTIONAL, "",  "PRN is appended to file name");
    readConfig(config, "outputfileTransmitterListGps",           fileNameTransmitterListGps,     Config::OPTIONAL, "",  "list of PRNs");
    readConfig(config, "outputfileTransmitterListGlonass",       fileNameTransmitterListGlonass, Config::OPTIONAL, "",  "list of PRNs");
    readConfig(config, "outputfileTransmitterListGalileo",       fileNameTransmitterListGalileo, Config::OPTIONAL, "",  "list of PRNs");
    readConfig(config, "outputfileTransmitterListBeiDou",        fileNameTransmitterListBeiDou,  Config::OPTIONAL, "",  "list of PRNs");
    readConfig(config, "outputfileTransmitterListQzss",          fileNameTransmitterListQzss,    Config::OPTIONAL, "",  "list of PRNs");
    readConfig(config, "outputfileTransmitterListIrnss",         fileNameTransmitterListIrnss,   Config::OPTIONAL, "",  "list of PRNs");
    readConfig(config, "inputfileAntex",                         fileNameAntex,                  Config::MUSTSET,  "", "");
    readConfig(config, "timeStart",                              timeStart,                      Config::OPTIONAL, "",  "ignore older antenna definitions");
    readConfig(config, "createZeroModel",                        setZero,                        Config::DEFAULT,  "0", "create empty antenna patterns");
    if(isCreateSchema(config)) return;

    // ==============================

    std::vector<PlatformGnssAntenna> antennaInfos;

    // Open file
    // ---------
    InFile file(fileNameAntex);
    file.exceptions(std::ios::badbit|std::ios::failbit);

    // Read Header
    // -----------
    std::string line, label;
    do
    {
      getLine(file, line, label);
    }
    while(!testLabel(label, "END OF HEADER", TRUE));

    // Read data
    // ---------
    for(;;)
    {
      PlatformGnssAntenna      antennaInfo;
      GnssAntennaDefinitionPtr antenna = GnssAntennaDefinitionPtr(new GnssAntennaDefinition);
      antennaInfo.antennaDef = antenna;

      if(!getLine(file, line, label, FALSE))
        break;
      testLabel(label, "START OF ANTENNA");
      getLine(file, line, label);

      // BLOCK IIA           G01                 G032      1992-079A TYPE / SERIAL NO
      // JAVRINGANT_DMT  SCIS                                        TYPE / SERIAL NO
      testLabel(label, "TYPE / SERIAL NO");
      antenna->name   = String::trim(line.substr( 0,16));
      antenna->radome = String::trim(line.substr(16,4));
      if(antenna->radome == "NONE")
        antenna->radome.clear();
      antenna->serial    = String::trim(line.substr(20,20));
      antennaInfo.serial = String::trim(line.substr(40,10)); // SVN
      antennaInfo.radome = String::trim(line.substr(50,10)); // COSPAR

      getLine(file, line, label);
      testLabel(label, "METH / BY / # / DATE");

      getLine(file, line, label);
      testLabel(label, "DAZI");
      UInt azimutCount = 0;
      Double dazi = String::toDouble(line.substr(2, 6));
      if(dazi!=0)
        azimutCount = static_cast<UInt>(360./dazi);

      getLine(file, line, label);
      testLabel(label, "ZEN1 / ZEN2 / DZEN");
      Double zen1 = String::toDouble(line.substr(2, 6));
      Double zen2 = String::toDouble(line.substr(8, 6));
      Double dzen = String::toDouble(line.substr(14, 6));
      UInt zenCount = static_cast<UInt>((zen2-zen1)/dzen+1);
      getLine(file, line, label);
      testLabel(label, "# OF FREQUENCIES");
      UInt freqCount = String::toInt(line.substr(0, 6));
      getLine(file, line, label);

      // Optionals
      // ---------
      Int  year,month,day,hour,min;
      Double sec;
      if(testLabel(label, "VALID FROM", TRUE))
      {
        year  = String::toInt(line.substr(0, 6));
        month = String::toInt(line.substr(6, 6));
        day   = String::toInt(line.substr(12, 6));
        hour  = String::toInt(line.substr(18, 6));
        min   = String::toInt(line.substr(24, 6));
        sec   = String::toDouble(line.substr(33, 13));
        antennaInfo.timeStart = date2time(year, month, day, hour, min, sec);
        getLine(file, line, label);
      }
      else
        antennaInfo.timeStart = Time();

      if(testLabel(label, "VALID UNTIL", TRUE))
      {
        year  = String::toInt(line.substr(0, 6));
        month = String::toInt(line.substr(6, 6));
        day   = String::toInt(line.substr(12, 6));
        hour  = String::toInt(line.substr(18, 6));
        min   = String::toInt(line.substr(24, 6));
        sec   = String::toDouble(line.substr(33, 13));
        if(sec >= 59)
          sec = 60; // to prevent gaps, since VALID UNTIL is inclusive in ANTEX but timeEnd is exclusive in GROOPS
        antennaInfo.timeEnd = date2time(year, month, day, hour, min, sec);
        getLine(file, line, label);
      }
      else
        antennaInfo.timeEnd = date2time(2500,1,1,0,0,0);

      if(testLabel(label, "SINEX CODE", TRUE))
      {
        getLine(file, line, label);
        while (testLabel(label, "SINEX CODE", TRUE))
        {
          logWarning<<"Multiple lines with <SINEX CODE> detected for antenna "<<antenna->name<<" "<<antenna->serial<<" "<<antennaInfo.serial<<Log::endl;
          getLine(file, line, label);
        }
      }

      while(testLabel(label, "COMMENT", TRUE))
      {
        if(!antennaInfo.comment.empty())
          antennaInfo.comment += "\n";
        antennaInfo.comment += String::trim(line.substr(0,60));
        getLine(file, line, label);
      }

      antenna->patterns.resize(freqCount);

      for(UInt i=0; i<freqCount; i++)
      {
        testLabel(label, "START OF FREQUENCY");

        antenna->patterns.at(i).type = GnssType();
        switch(line.at(3))
        {
          case ' ': break;
          case 'G': antenna->patterns.at(i).type += GnssType::GPS;      break;
          case 'R': antenna->patterns.at(i).type += GnssType::GLONASS;  break;
          case 'E': antenna->patterns.at(i).type += GnssType::GALILEO;  break;
          case 'C': antenna->patterns.at(i).type += GnssType::BDS;      break;
          case 'S': antenna->patterns.at(i).type += GnssType::SBAS;     break;
          case 'I': antenna->patterns.at(i).type += GnssType::IRNSS;    break;
          case 'J': antenna->patterns.at(i).type += GnssType::QZSS;     break;
          default:
            logWarning<<"Unknown satellite system: '"<<line<<"'"<<Log::endl;
        }
        switch(String::toInt(line.substr(4, 2)))
        {
          case 0: break;
          case 1: antenna->patterns.at(i).type += GnssType::L1;  break;
          case 2: antenna->patterns.at(i).type += GnssType::L2;  break;
          case 3: antenna->patterns.at(i).type += GnssType::G3;  break;
          case 4: antenna->patterns.at(i).type += GnssType::G1a; break;
          case 5: antenna->patterns.at(i).type += GnssType::L5;  break;
          case 6: antenna->patterns.at(i).type += GnssType::E6;  break;
          case 7: antenna->patterns.at(i).type += GnssType::E5b; break;
          case 8: antenna->patterns.at(i).type += GnssType::E5;  break;
          case 9: antenna->patterns.at(i).type += GnssType::S9;  break;
          default:
            logWarning<<"Unknown frequency: '"<<line<<"'"<<Log::endl;
        }
        getLine(file, line, label);

        testLabel(label, "NORTH / EAST / UP");
        Double x = String::toDouble(line.substr(0, 10));
        Double y = String::toDouble(line.substr(10, 10));
        Double z = String::toDouble(line.substr(20, 10));
        // from millimeters to meters
        if(!setZero) antenna->patterns.at(i).offset = 1e-3 * Vector3d(x,y,z);

        antenna->patterns.at(i).dZenit  = Angle(dzen*DEG2RAD);
        if(azimutCount)
          antenna->patterns.at(i).pattern = Matrix(azimutCount, zenCount);
        else
          antenna->patterns.at(i).pattern = Matrix(1, zenCount);

        // read NOAZI
        getLine(file, line, label);
        for(UInt s=0; s<zenCount; s++)
          if(!setZero) antenna->patterns.at(i).pattern(0,s) = 1e-3 * String::toDouble(line.substr(8+8*s, 8));

        // Azimuth dependent values
        getLine(file, line, label);
        if(azimutCount)
          for(UInt z=0; z<azimutCount+1; z++)
          {
            for(UInt s=0; s<zenCount; s++)
              if(!setZero) antenna->patterns.at(i).pattern(z%azimutCount,s) = 1e-3 * String::toDouble(line.substr(8+8*s, 8));
            getLine(file, line, label);
          }

        testLabel(label, "END OF FREQUENCY");
        getLine(file, line, label);

        // skip optional frequency RMS block
        if(testLabel(label, "START OF FREQ RMS", TRUE))
        {
          for(;;)
          {
            getLine(file, line, label);
            if(testLabel(label, "END OF FREQ RMS", TRUE))
              break;
          }
          getLine(file, line, label);
        }
      }
      testLabel(label, "END OF ANTENNA");

      if((antennaInfo.timeEnd == date2time(2500,1,1,0,0,0)) || (antennaInfo.timeEnd > timeStart))
        antennaInfos.push_back(antennaInfo);
    } // for(antenna section)

    // ==============================

    // Satellites: Shift mean antenna offsets to PlatformEquipment::position
    // ---------------------------------------------------------------------
    for(auto &antennaInfo : antennaInfos)
      if(!antennaInfo.serial.empty()) // with SVN
      {
        // transformation to left-handed antenna system
        antennaInfo.local2antennaFrame = flipY() * rotaryZ(Angle(PI/2));

        // set the PlatformEquipment::position as the mean of offsets of all patterns
        auto antenna = antennaInfo.antennaDef;
        for(UInt i=0; i<antenna->patterns.size(); i++)
          antennaInfo.position += (1./antenna->patterns.size()) * antenna->patterns.at(i).offset;
        // shift the origin of offsets from the reference point (CoM) to PlatformEquipment::position
        for(UInt i=0; i<antenna->patterns.size(); i++)
          antenna->patterns.at(i).offset = antennaInfo.local2antennaFrame.transform(antenna->patterns.at(i).offset - antennaInfo.position);
      }

    // ==============================

    // Distinguish between different antennas for the same SVN by using different radome names
    // ---------------------------------------------------------------------------------------
    // sort SVNs, timeStart
    std::stable_sort(antennaInfos.begin(), antennaInfos.end(), [](auto &a1, auto &a2) {return (a1.serial != a2.serial) ? (a1.serial < a2.serial) : (a1.timeEnd < a2.timeEnd);});

    for(UInt i=antennaInfos.size()-1; i-->0;)
      if(!antennaInfos.at(i).serial.empty() && (antennaInfos.at(i).serial == antennaInfos.at(i+1).serial)) // with same SVN as before
      {
        if(equal(antennaInfos.at(i).antennaDef, antennaInfos.at(i+1).antennaDef))
          antennaInfos.at(i).radome = antennaInfos.at(i+1).radome; // copy COSPAR + valid until
        else
          antennaInfos.at(i).radome = antennaInfos.at(i).radome.substr(0, 9)+" (valid unil "s+antennaInfos.at(i).timeEnd.dateTimeStr()+")"s;
      }

    // ==============================

    std::vector<GnssAntennaDefinitionPtr> antennaListStation, antennaListTransmitter;
    std::vector<Platform>                 transmitterList;

    for(auto &antennaInfo : antennaInfos)
    {
      auto antenna = antennaInfo.antennaDef;

      if(!antennaInfo.serial.empty()) // satellites with SVN
      {
        // Satellites
        std::string prn  = antenna->serial;
        antenna->serial  = antennaInfo.serial; // SVN;
        antenna->radome  = antennaInfo.radome; // COSPAR;
        antennaInfo.name = antenna->name;      // BLOCK name
        addAntenna(antenna, antennaListTransmitter);

        if     (antenna->name.find("BLOCK")   != std::string::npos) addTransmitter(antennaInfo, "GPS",     prn, transmitterList);
        else if(antenna->name.find("GLONASS") != std::string::npos) addTransmitter(antennaInfo, "GLONASS", prn, transmitterList);
        else if(antenna->name.find("GALILEO") != std::string::npos) addTransmitter(antennaInfo, "GALILEO", prn, transmitterList);
        else if(antenna->name.find("BEIDOU")  != std::string::npos) addTransmitter(antennaInfo, "BEIDOU",  prn, transmitterList);
        else if(antenna->name.find("QZSS")    != std::string::npos) addTransmitter(antennaInfo, "QZSS",    prn, transmitterList);
        else if(antenna->name.find("IRNSS")   != std::string::npos) addTransmitter(antennaInfo, "IRNSS",   prn, transmitterList);
        else
          logWarning<<"Unknown satellite: "<<prn<<" "<<antenna->serial<<" ("<<antenna->name<<", "<<antenna->radome<<")"<<Log::endl;
      }
      else if(antennaInfo.timeStart == Time() && antennaInfo.timeEnd == date2time(2500,1,1,0,0,0)) // stations
      {
        addAntenna(antenna, antennaListStation);
      }
      else
        logWarning<<"antenna '"<<antenna->name<<"."<<antenna->serial<<"."<<antenna->radome<<"' with validity period -> skipped"<<Log::endl;
    }

    // ==============================

    // save
    // ----
    if(!fileNameAntennaStation.empty() && antennaListStation.size())
    {
      logStatus<<"save antennas to <"<<fileNameAntennaStation<<">"<< Log::endl;
      writeFileGnssAntennaDefinition(fileNameAntennaStation, antennaListStation);
    }

    if(!fileNameAntennaTransmitter.empty() && antennaListTransmitter.size())
    {
      logStatus<<"save antennas to <"<<fileNameAntennaTransmitter<<">"<< Log::endl;
      writeFileGnssAntennaDefinition(fileNameAntennaTransmitter, antennaListTransmitter);
    }

    if(!fileNameTransmitterInfo.empty() && transmitterList.size())
    {
      logStatus<<"save transmitter info to <"<<fileNameTransmitterInfo<<">"<< Log::endl;
      for(const auto &transmitter : transmitterList)
        writeFilePlatform(fileNameTransmitterInfo.appendBaseName("."+transmitter.markerNumber), transmitter);
    }

    auto writeTransmitterList = [&](const FileName &fileNameTransmitterList, const std::string &markerName)
    {
      if(!fileNameTransmitterList.empty() && transmitterList.size())
      {
        logStatus<<"save transmitter list to <"<<fileNameTransmitterList<<">"<< Log::endl;
        std::set<std::string> list;
        for(const auto &transmitter : transmitterList)
          if(transmitter.markerName == markerName)
            list.insert(transmitter.markerNumber);
        writeFileStringList(fileNameTransmitterList, std::vector<std::string>(list.begin(), list.end()));
      }
    };

    writeTransmitterList(fileNameTransmitterListGps,     "GPS");
    writeTransmitterList(fileNameTransmitterListGlonass, "GLONASS");
    writeTransmitterList(fileNameTransmitterListGalileo, "GALILEO");
    writeTransmitterList(fileNameTransmitterListBeiDou,  "BEIDOU");
    writeTransmitterList(fileNameTransmitterListQzss,    "QZSS");
    writeTransmitterList(fileNameTransmitterListIrnss,   "IRNSS");
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/

Bool GnssAntex2AntennaDefinition::getLine(InFile &file, std::string &line, std::string &label, Bool throwException) const
{
  try
  {
    std::getline(file, line);
    if(line.size()<80)
      line.resize(80,' ');
    label = line.substr(60,20);
    return TRUE;
  }
  catch(std::exception &/*e*/)
  {
    if(throwException)
      throw(Exception(std::string("IO Error: '")+line+"'"));
    line.clear();
    line.resize(80,' ');
    label = line.substr(60,20);
    return FALSE;
  }
}

/***********************************************/

Bool GnssAntex2AntennaDefinition::testLabel(const std::string &labelInLine, const std::string &label, Bool optional) const
{
  if(labelInLine.find(label)!=std::string::npos)
    return TRUE;
  if(optional)
    return FALSE;
  throw(Exception(std::string("In Line '")+labelInLine+"' label '"+label+"' expected\n"));
}

/***********************************************/

Bool GnssAntex2AntennaDefinition::equal(GnssAntennaDefinitionPtr antenna, GnssAntennaDefinitionPtr antenna2)
{
  try
  {
    // check if antennas are identical
    if(antenna->patterns.size() != antenna2->patterns.size())
      return FALSE;
    for(UInt idPattern=0; idPattern<antenna->patterns.size(); idPattern++)
    {
      if((antenna->patterns.at(idPattern).type              != antenna2->patterns.at(idPattern).type) ||
         (antenna->patterns.at(idPattern).pattern.rows()    != antenna2->patterns.at(idPattern).pattern.rows()) ||
         (antenna->patterns.at(idPattern).pattern.columns() != antenna2->patterns.at(idPattern).pattern.columns()))
        return FALSE;
      if(((antenna->patterns.at(idPattern).offset - antenna2->patterns.at(idPattern).offset).r() > 0.0001) ||
         (maxabs(antenna->patterns.at(idPattern).pattern - antenna2->patterns.at(idPattern).pattern) > 0.0001))
        return FALSE;
    }
    return TRUE;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssAntex2AntennaDefinition::addAntenna(GnssAntennaDefinitionPtr antenna, std::vector<GnssAntennaDefinitionPtr> &antennaList)
{
  try
  {
    auto iter = std::find_if(antennaList.begin(), antennaList.end(), [&](auto &a){return ((a->name == antenna->name) && (a->serial == antenna->serial) && (a->radome == antenna->radome));});
    if(iter == antennaList.end())
      antennaList.push_back(antenna);
    else if(!equal(antenna, *iter))
      logWarning<<"antenna '"<<antenna->name<<"."<<antenna->serial<<"."<<antenna->radome<<"' already in the list with different values -> skipped"<<Log::endl;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssAntex2AntennaDefinition::addTransmitter(const PlatformGnssAntenna &antennaInfo, const std::string &markerName, const std::string &markerNumber, std::vector<Platform> &transmitterList)
{
  try
  {
    const UInt idSat = std::distance(transmitterList.begin(), std::find_if(transmitterList.begin(), transmitterList.end(),
                                                                           [&](auto &t){return t.markerNumber == markerNumber;}));
    if(idSat >= transmitterList.size()) // new PRN?
    {
      Platform satelliteGps;
      satelliteGps.markerName   = markerName;
      satelliteGps.markerNumber = markerNumber;
      transmitterList.push_back(satelliteGps);
    }

    PlatformGnssReceiver receiverInfo;
    receiverInfo.name      = antennaInfo.name;
    receiverInfo.serial    = antennaInfo.serial;
    receiverInfo.timeStart = antennaInfo.timeStart;
    receiverInfo.timeEnd   = antennaInfo.timeEnd;
    transmitterList.at(idSat).equipments.push_back(std::make_shared<PlatformGnssAntenna>(antennaInfo));
    transmitterList.at(idSat).equipments.push_back(std::make_shared<PlatformGnssReceiver>(receiverInfo));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
