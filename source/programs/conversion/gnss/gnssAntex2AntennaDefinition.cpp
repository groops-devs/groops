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
to \configFile{transmitterInfo}{platform}, \configFile{antennaDefinition}{gnssAntennaDefinition}, \configFile{svnBlockTable}{stringTable},
and \configFile{transmitterList}{stringList} files for the respective GNSS and for the list of ground station antennas.
)";

/***********************************************/

#include "programs/program.h"
#include "base/string.h"
#include "inputOutput/file.h"
#include "files/fileGnssStationInfo.h"
#include "files/fileStringTable.h"

/***** CLASS ***********************************/

/** @brief Converts IGS ANTEX file to GNSS metadata and antenna definition files.
* @ingroup programsConversionGroup */
class GnssAntex2AntennaDefinition
{
  Bool getLine(InFile &file, std::string &line, std::string &label, Bool throwException=FALSE) const;
  Bool testLabel(const std::string &labelInLine, const std::string &label, Bool optional=TRUE) const;
  void addTransmitter(const GnssAntennaInfo &antennaInfo, const std::string &markerName, const std::string &markerNumber, std::vector<GnssStationInfo> &transmitterList);
  void addAntenna(const GnssAntennaDefinitionPtr &antenna, std::vector<GnssAntennaDefinitionPtr> &antennaList);

  enum Type
  {
    OTHER, STATION, GPS, GLONASS, GALILEO, BEIDOU, QZSS
  };

public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(GnssAntex2AntennaDefinition, SINGLEPROCESS, "Converts IGS ANTEX file to GNSS metadata and antenna definition files.", Conversion, Gnss)

/***********************************************/

void GnssAntex2AntennaDefinition::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName outFileNameAntennaStation, outFileNameAntennaTransmitter;
    FileName outFileNameTransmitterInfo;
    FileName outFileNameSvnBlockTableGps, outFileNameSvnBlockTableGlonass, outFileNameSvnBlockTableGalileo, outFileNameSvnBlockTableBeiDou, outFileNameSvnBlockTableQzss;
    FileName outFileNameTransmitterListGps, outFileNameTransmitterListGlonass, outFileNameTransmitterListGalileo, outFileNameTransmitterListBeiDou, outFileNameTransmitterListQzss;
    FileName inFileName;
    Time     timeStart;
    Bool     setZero;

    readConfig(config, "outputfileAntennaDefinitionStation",      outFileNameAntennaStation,         Config::OPTIONAL, "",  "antenna center variations");
    readConfig(config, "outputfileAntennaDefinitionTransmitter",  outFileNameAntennaTransmitter,     Config::OPTIONAL, "",  "antenna center variations");
    readConfig(config, "outputfileTransmitterInfo",               outFileNameTransmitterInfo,        Config::OPTIONAL, "",  "PRN is appended to file name");
    readConfig(config, "outputfileSvnBlockTableGps",              outFileNameSvnBlockTableGps,       Config::OPTIONAL, "",  "SVN to satellite block mapping");
    readConfig(config, "outputfileSvnBlockTableGlonass",          outFileNameSvnBlockTableGlonass,   Config::OPTIONAL, "",  "SVN to satellite block mapping");
    readConfig(config, "outputfileSvnBlockTableGalileo",          outFileNameSvnBlockTableGalileo,   Config::OPTIONAL, "",  "SVN to satellite block mapping");
    readConfig(config, "outputfileSvnBlockTableBeiDou",           outFileNameSvnBlockTableBeiDou,    Config::OPTIONAL, "",  "SVN to satellite block mapping");
    readConfig(config, "outputfileSvnBlockTableQzss",             outFileNameSvnBlockTableQzss,      Config::OPTIONAL, "",  "SVN to satellite block mapping");
    readConfig(config, "outputfileTransmitterListGps",            outFileNameTransmitterListGps,     Config::OPTIONAL, "",  "list of PRNs");
    readConfig(config, "outputfileTransmitterListGlonass",        outFileNameTransmitterListGlonass, Config::OPTIONAL, "",  "list of PRNs");
    readConfig(config, "outputfileTransmitterListGalileo",        outFileNameTransmitterListGalileo, Config::OPTIONAL, "",  "list of PRNs");
    readConfig(config, "outputfileTransmitterListBeiDou",         outFileNameTransmitterListBeiDou,  Config::OPTIONAL, "",  "list of PRNs");
    readConfig(config, "outputfileTransmitterListQzss",           outFileNameTransmitterListQzss,    Config::OPTIONAL, "",  "list of PRNs");
    readConfig(config, "inputfileAntex",                          inFileName,                        Config::MUSTSET,  "", "");
    readConfig(config, "timeStart",                               timeStart,                         Config::OPTIONAL, "",  "ignore older antenna definitions");
    readConfig(config, "createZeroModel",                         setZero,                           Config::DEFAULT,  "0", "create empty antenna patterns");
    if(isCreateSchema(config)) return;

    // ==============================

    std::vector<GnssAntennaDefinitionPtr> antennaListStation, antennaListGps, antennaListGlonass, antennaListGalileo, antennaListBeiDou, antennaListQzss;
    std::vector<GnssStationInfo>          transmitterListGps, transmitterListGlonass, transmitterListGalileo, transmitterListBeiDou, transmitterListQzss;

    // Open file
    // ---------
    InFile file(inFileName);
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
      std::string atxSVN, atxCOSPAR;
      GnssAntennaInfo antennaInfo;

      if(!getLine(file, line, label, FALSE))
        break;
      testLabel(label, "START OF ANTENNA");
      getLine(file, line, label);

      testLabel(label, "TYPE / SERIAL NO");
      antennaInfo.name   = String::trim(line.substr( 0,16));
      antennaInfo.radome = String::trim(line.substr(16,4));
      if(antennaInfo.radome == "NONE")
        antennaInfo.radome.clear();
      antennaInfo.serial = String::trim(line.substr(20,20));
      atxSVN    = String::trim(line.substr(40,10));
      atxCOSPAR = String::trim(line.substr(50,10));
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
        sec   = String::toDouble(line.substr(33, 43));
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
        sec   = String::toDouble(line.substr(33, 43));
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
      }

      while(testLabel(label, "COMMENT", TRUE))
      {
        if(!antennaInfo.comment.empty())
          antennaInfo.comment += "\n";
        antennaInfo.comment += String::trim(line.substr(0,60));
        getLine(file, line, label);
      }

      GnssAntennaDefinitionPtr antenna = GnssAntennaDefinitionPtr(new GnssAntennaDefinition);
      antennaInfo.antennaDef = antenna;
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

        // Azimut dependent values
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

      if(antennaInfo.timeEnd != date2time(2500,1,1,0,0,0) && antennaInfo.timeEnd <= timeStart)
        continue;

      Type type = OTHER;
      if(antennaInfo.name.find("BLOCK")   != std::string::npos) type = GPS;
      if(antennaInfo.name.find("GLONASS") != std::string::npos) type = GLONASS;
      if(antennaInfo.name.find("GALILEO") != std::string::npos) type = GALILEO;
      if(antennaInfo.name.find("BEIDOU")  != std::string::npos) type = BEIDOU;
      if(antennaInfo.name.find("QZSS")    != std::string::npos) type = QZSS;
      if(antennaInfo.timeStart == Time() && antennaInfo.timeEnd == date2time(2500,1,1,0,0,0)) type = STATION;

      // GNSS satellite? ==> create transmitter info
      if(type == GPS || type == GLONASS || type == GALILEO || type == BEIDOU || type == QZSS)
      {
        std::string prn = antennaInfo.serial;
        antennaInfo.serial = atxSVN;
        antennaInfo.radome = atxCOSPAR;

        // transformation to left-handed antenna system
        antennaInfo.local2antennaFrame = flipY() * rotaryZ(Angle(PI/2));

        // shift from GnssAntennaPattern::offset to GnssAntennaInfo::position
        for(UInt i=0; i<antenna->patterns.size(); i++)
          antennaInfo.position += (1./antenna->patterns.size()) * antenna->patterns.at(i).offset;
        for(UInt i=0; i<antenna->patterns.size(); i++)
          antenna->patterns.at(i).offset = antennaInfo.local2antennaFrame.transform(antenna->patterns.at(i).offset - antennaInfo.position);

        if(type == GPS)     addTransmitter(antennaInfo, "GPS",     prn, transmitterListGps);
        if(type == GLONASS) addTransmitter(antennaInfo, "GLONASS", prn, transmitterListGlonass);
        if(type == GALILEO) addTransmitter(antennaInfo, "GALILEO", prn, transmitterListGalileo);
        if(type == BEIDOU)  addTransmitter(antennaInfo, "BEIDOU",  prn, transmitterListBeiDou);
        if(type == QZSS)    addTransmitter(antennaInfo, "QZSS",    prn, transmitterListQzss);
      }

      antenna->name    = antennaInfo.name;
      antenna->serial  = antennaInfo.serial;
      antenna->radome  = antennaInfo.radome;
      antenna->comment = antennaInfo.comment;

      if(type == STATION) addAntenna(antenna, antennaListStation);
      if(type == GPS)     addAntenna(antenna, antennaListGps);
      if(type == GLONASS) addAntenna(antenna, antennaListGlonass);
      if(type == GALILEO) addAntenna(antenna, antennaListGalileo);
      if(type == BEIDOU)  addAntenna(antenna, antennaListBeiDou);
      if(type == QZSS)    addAntenna(antenna, antennaListQzss);
    } // for(antenna section)

    // ==============================

    // save
    // ----
    auto writeAntenna = [] (const FileName &outFileNameAntenna, const std::vector<GnssAntennaDefinitionPtr> &antennaList)
    {
      if(!outFileNameAntenna.empty() && antennaList.size())
      {
        logStatus<<"save antennas to <"<<outFileNameAntenna<<">"<< Log::endl;
        writeFileGnssAntennaDefinition(outFileNameAntenna, antennaList);
      }
    };

    // sort GPS antenna list according to SVN number (ascending)
    std::sort(antennaListGps.begin(), antennaListGps.end(), [](GnssAntennaDefinitionPtr a, GnssAntennaDefinitionPtr b) {return std::stoi(a->serial.substr(1,3)) < std::stoi(b->serial.substr(1,3));});

    writeAntenna(outFileNameAntennaStation, antennaListStation);
    std::vector<GnssAntennaDefinitionPtr> antennaListTransmitter;
    antennaListTransmitter.insert(antennaListTransmitter.end(), antennaListGps.begin(),     antennaListGps.end());
    antennaListTransmitter.insert(antennaListTransmitter.end(), antennaListGlonass.begin(), antennaListGlonass.end());
    antennaListTransmitter.insert(antennaListTransmitter.end(), antennaListGalileo.begin(), antennaListGalileo.end());
    antennaListTransmitter.insert(antennaListTransmitter.end(), antennaListBeiDou.begin(),  antennaListBeiDou.end());
    antennaListTransmitter.insert(antennaListTransmitter.end(), antennaListQzss.begin(),    antennaListQzss.end());
    writeAntenna(outFileNameAntennaTransmitter, antennaListTransmitter);

    auto writeTransmitterInfo = [] (const FileName &outFileNameTransmitterInfo, const std::vector<GnssStationInfo> &transmitterList)
    {
      if(!outFileNameTransmitterInfo.empty() && transmitterList.size())
      {
        logStatus<<"save transmitter info to <"<<outFileNameTransmitterInfo<<">"<< Log::endl;
        for(const auto &transmitter : transmitterList)
          writeFileGnssStationInfo(outFileNameTransmitterInfo.appendBaseName("."+transmitter.markerNumber), transmitter);
      }
    };

    std::vector<GnssStationInfo> transmitterList;
    transmitterList.insert(transmitterList.end(), transmitterListGps.begin(),     transmitterListGps.end());
    transmitterList.insert(transmitterList.end(), transmitterListGlonass.begin(), transmitterListGlonass.end());
    transmitterList.insert(transmitterList.end(), transmitterListGalileo.begin(), transmitterListGalileo.end());
    transmitterList.insert(transmitterList.end(), transmitterListBeiDou.begin(),  transmitterListBeiDou.end());
    transmitterList.insert(transmitterList.end(), transmitterListQzss.begin(),    transmitterListQzss.end());
    writeTransmitterInfo(outFileNameTransmitterInfo, transmitterList);

    auto writeSvnBlockTable = [] (const FileName &outFileNameSvnBlockTable, const std::vector<GnssAntennaDefinitionPtr> &antennaList)
    {
      if(!outFileNameSvnBlockTable.empty() && antennaList.size())
      {
        logStatus<<"save SVN-to-block table to <"<<outFileNameSvnBlockTable<<">"<< Log::endl;
        std::vector<std::vector<std::string>> table;
        for(const auto &antenna : antennaList)
           table.push_back({antenna->serial, antenna->name});
        std::sort(table.begin(), table.end(), [](const std::vector<std::string> &a, const std::vector<std::string> &b){ return a.at(0) < b.at(0); });
        writeFileStringTable(outFileNameSvnBlockTable, table);
      }
    };

    writeSvnBlockTable(outFileNameSvnBlockTableGps,     antennaListGps);
    writeSvnBlockTable(outFileNameSvnBlockTableGlonass, antennaListGlonass);
    writeSvnBlockTable(outFileNameSvnBlockTableGalileo, antennaListGalileo);
    writeSvnBlockTable(outFileNameSvnBlockTableBeiDou,  antennaListBeiDou);
    writeSvnBlockTable(outFileNameSvnBlockTableQzss,    antennaListQzss);

    auto writeTransmitterList = [] (const FileName &outFileNameTransmitterList, const std::vector<GnssStationInfo> &transmitterList)
    {
      if(!outFileNameTransmitterList.empty() && transmitterList.size())
      {
        logStatus<<"save transmitter list to <"<<outFileNameTransmitterList<<">"<< Log::endl;
        std::set<std::string> list;
        for(const auto &transmitter : transmitterList)
          list.insert(transmitter.markerNumber);
        writeFileStringList(outFileNameTransmitterList, std::vector<std::string>(list.begin(), list.end()));
      }
    };

    writeTransmitterList(outFileNameTransmitterListGps,     transmitterListGps);
    writeTransmitterList(outFileNameTransmitterListGlonass, transmitterListGlonass);
    writeTransmitterList(outFileNameTransmitterListGalileo, transmitterListGalileo);
    writeTransmitterList(outFileNameTransmitterListBeiDou,  transmitterListBeiDou);
    writeTransmitterList(outFileNameTransmitterListQzss,    transmitterListQzss);
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

void GnssAntex2AntennaDefinition::addTransmitter(const GnssAntennaInfo &antennaInfo, const std::string &markerName, const std::string &markerNumber, std::vector<GnssStationInfo> &transmitterList)
{
  try
  {
    UInt idSat;
    for(idSat=0; idSat<transmitterList.size(); idSat++)
      if(transmitterList.at(idSat).markerNumber == markerNumber)
        break;

    if(idSat>=transmitterList.size()) // new PRN?
    {
      GnssStationInfo satelliteGps;
      satelliteGps.markerName   = markerName;
      satelliteGps.markerNumber = markerNumber;
      transmitterList.push_back(satelliteGps);
    }

    transmitterList.at(idSat).antenna.push_back(antennaInfo);

    GnssReceiverInfo receiverInfo;
    receiverInfo.name      = antennaInfo.name;
    receiverInfo.serial    = antennaInfo.serial;
    receiverInfo.timeStart = antennaInfo.timeStart;
    receiverInfo.timeEnd   = antennaInfo.timeEnd;
    transmitterList.at(idSat).receiver.push_back(receiverInfo);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssAntex2AntennaDefinition::addAntenna(const GnssAntennaDefinitionPtr &antenna, std::vector<GnssAntennaDefinitionPtr> &antennaList)
{
  try
  {
    auto iter = std::find_if(antennaList.begin(), antennaList.end(), [&](auto &a){return ((a->name == antenna->name) && (a->serial == antenna->serial) && (a->radome == antenna->radome));});
    if(iter != antennaList.end())
    {
      // check if antennas are identical
      GnssAntennaDefinitionPtr antenna2 = *iter;
      Bool equal = (antenna2->patterns.size() == antenna->patterns.size());
      if(equal)
        for(UInt idPattern=0; idPattern<antenna->patterns.size(); idPattern++)
        {
          equal =  (antenna->patterns.at(idPattern).type              == antenna2->patterns.at(idPattern).type)
                && (antenna->patterns.at(idPattern).pattern.rows()    == antenna2->patterns.at(idPattern).pattern.rows())
                && (antenna->patterns.at(idPattern).pattern.columns() == antenna2->patterns.at(idPattern).pattern.columns());
          if(!equal)
            break;
          equal = ((antenna->patterns.at(idPattern).offset - antenna2->patterns.at(idPattern).offset).r() < 0.0001)
                && (maxabs(antenna->patterns.at(idPattern).pattern - antenna2->patterns.at(idPattern).pattern) < 0.0001);
          if(!equal)
            break;
        }
      if(!equal)
        logWarning<<"antenna '"<<antenna->name<<"."<<antenna->serial<<"."<<antenna->radome<<"' already in the list with different values -> skipped"<<Log::endl;
    }
    else
      antennaList.push_back(antenna);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
