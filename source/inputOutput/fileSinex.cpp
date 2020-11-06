/***********************************************/
/**
* @file fileSinex.cpp
*
* @brief SINEX file representation.
*
* @author Sebastian Strasser
* @date 2017-05-15
*
*/
/***********************************************/

#include "base/import.h"
#include "base/string.h"
#include "inputOutput/logging.h"
#include "inputOutput/file.h"
#include "inputOutput/system.h"
#include "files/fileGnssStationInfo.h"
#include "config/config.h"
#include "fileSinex.h"

/***********************************************/

Sinex::BlockType Sinex::blockType(const std::string &label)
{
  try
  {
    if(     label               == "FILE/REFERENCE")                  return BlockType::FILE_REFERENCE;
    else if(label               == "FILE/COMMENT")                    return BlockType::FILE_COMMENT;
    else if(label               == "INPUT/HISTORY")                   return BlockType::INPUT_HISTORY;
    else if(label               == "INPUT/FILES")                     return BlockType::INPUT_FILES;
    else if(label               == "INPUT/ACKNOWLEDGEMENTS" ||
            label               == "INPUT/ACKNOWLEDGMENTS")           return BlockType::INPUT_ACKNOWLEDGEMENTS;
    else if(label               == "NUTATION/DATA")                   return BlockType::NUTATION_DATA;
    else if(label               == "PRECESSION/DATA")                 return BlockType::PRECESSION_DATA;
    else if(label               == "SOURCE/ID")                       return BlockType::SOURCE_ID;
    else if(label               == "SITE/ID")                         return BlockType::SITE_ID;
    else if(label               == "SITE/DATA")                       return BlockType::SITE_DATA;
    else if(label               == "SITE/RECEIVER")                   return BlockType::SITE_RECEIVER;
    else if(label               == "SITE/ANTENNA")                    return BlockType::SITE_ANTENNA;
    else if(label               == "SITE/GPS_PHASE_CENTER")           return BlockType::SITE_GPS_PHASE_CENTER;
    else if(label               == "SITE/GAL_PHASE_CENTER")           return BlockType::SITE_GAL_PHASE_CENTER;
    else if(label               == "SITE/ECCENTRICITY")               return BlockType::SITE_ECCENTRICITY;
    else if(label               == "SATELLITE/ID")                    return BlockType::SATELLITE_ID;
    else if(label               == "SATELLITE/PHASE_CENTER")          return BlockType::SATELLITE_PHASE_CENTER;
    else if(label               == "BIAS/EPOCHS")                     return BlockType::BIAS_EPOCHS;
    else if(label               == "BIAS/DESCRIPTION")                return BlockType::BIAS_DESCRIPTION;
    else if(label               == "BIAS/SOLUTION")                   return BlockType::BIAS_SOLUTION;
    else if(label               == "SOLUTION/EPOCHS")                 return BlockType::SOLUTION_EPOCHS;
    else if(label               == "SOLUTION/DISCONTINUITY")          return BlockType::SOLUTION_DISCONTINUITY;
    else if(label               == "SOLUTION/STATISTICS")             return BlockType::SOLUTION_STATISTICS;
    else if(label               == "SOLUTION/ESTIMATE")               return BlockType::SOLUTION_ESTIMATE;
    else if(label               == "SOLUTION/APRIORI")                return BlockType::SOLUTION_APRIORI;
    else if(label.substr(0, 24) == "SOLUTION/MATRIX_ESTIMATE")        return BlockType::SOLUTION_MATRIX_ESTIMATE;
    else if(label.substr(0, 23) == "SOLUTION/MATRIX_APRIORI")         return BlockType::SOLUTION_MATRIX_APRIORI;
    else if(label               == "SOLUTION/NORMAL_EQUATION_VECTOR") return BlockType::SOLUTION_NORMAL_EQUATION_VECTOR;
    else if(label.substr(0, 31) == "SOLUTION/NORMAL_EQUATION_MATRIX") return BlockType::SOLUTION_NORMAL_EQUATION_MATRIX;
    else if(label               == "SATELLITE/IDENTIFIER")            return BlockType::SATELLITE_IDENTIFIER;
    else if(label               == "SATELLITE/PRN")                   return BlockType::SATELLITE_PRN;
    else if(label               == "SATELLITE/FREQUENCY_CHANNEL")     return BlockType::SATELLITE_FREQUENCY_CHANNEL;
    else if(label               == "SATELLITE/MASS")                  return BlockType::SATELLITE_MASS;
    else if(label               == "SATELLITE/COM")                   return BlockType::SATELLITE_COM;
    else if(label               == "SATELLITE/ECCENTRICITY")          return BlockType::SATELLITE_ECCENTRICITY;
    else if(label               == "SATELLITE/TX_POWER")              return BlockType::SATELLITE_TX_POWER;
    else
      throw(Exception("Unknown block type: " + label));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

std::string Sinex::format(Double value, UInt length, UInt precision)
{
  try
  {
    std::string s = value%("%"+length%"%i"s+"."+precision%"%i"s+"f");
    if(s.size() > length && s.substr(0,3) == "-0.")
      return "-." + s.substr(3, s.size());
    return s;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void Sinex::readConfigHeader(Config &config)
{
  try
  {
    Time        timeStart, timeEnd;
    std::string agencyCode, observationCode, constraintCode, solutionContent;
    std::string description, output, contact, software, hardware, input;
    std::vector<std::string> comment;
    FileName fileNameComment;

    readConfigSequence(config, "sinexHeader", Config::MUSTSET, "", "");
    readConfig(config, "agencyCode",              agencyCode,       Config::OPTIONAL, "TUG",    "identify the agency providing the data");
    readConfig(config, "timeStart",               timeStart,        Config::OPTIONAL, "",       "start time of the data");
    readConfig(config, "timeEnd",                 timeEnd,          Config::OPTIONAL, "",       "end time of the data ");
    readConfig(config, "observationCode",         observationCode,  Config::OPTIONAL, "C",      "technique used to generate the SINEX solution");
    readConfig(config, "constraintCode",          constraintCode,   Config::OPTIONAL, "2",      "0: tight constraint, 1: siginficant constraint, 2: unconstrained");
    readConfig(config, "solutionContent",         solutionContent,  Config::OPTIONAL, "",       "solution types contained in the SINEX solution (S O E T C A)");
    readConfig(config, "description",             description,      Config::OPTIONAL, "",       "organizitions gathering/alerting the file contents");
    readConfig(config, "contact",                 contact,          Config::OPTIONAL, "",       "Address of the relevant contact. e-mail");
    readConfig(config, "output",                  output,           Config::OPTIONAL, "",       "Description of the file contents");
    readConfig(config, "input",                   input,            Config::OPTIONAL, "",       "Brief description of the input used to generate this solution");
    readConfig(config, "software",                software,         Config::OPTIONAL, "GROOPS", "Software used to generate the file");
    readConfig(config, "hardware",                hardware,         Config::OPTIONAL, "",       "Computer hardware on which above software was run");
    readConfig(config, "inputfileComment",        fileNameComment,  Config::OPTIONAL, "",       "comments in the comment block from a file (truncated at 80 characters)");
    readConfig(config, "comment",                 comment,          Config::OPTIONAL, "",       "comments in the comment block");
    endSequence(config);
    if(isCreateSchema(config)) return;

    // header line
    Time   timeCurrent = System::now();
    std::stringstream ss;
    ss << "%=SNX 2.02 " << std::setw(3) << agencyCode.substr(0,3) << " " << time2str(timeCurrent) << " " << std::setw(3) << agencyCode.substr(0,3);
    ss << " " << time2str(timeStart) << " " << time2str(timeEnd) << " " << std::setw(1) << observationCode.substr(0,1) << " 00000";
    ss << " " << std::setw(1) << constraintCode.substr(0,1) << " " << solutionContent.substr(0,12);
    _header = ss.str();

    // reference block
    SinexTextPtr sinexTextReference = addBlock<SinexText>("FILE/REFERENCE");

    auto refLine = [] (const std::string &name, const std::string &info) -> std::string
    {
      std::stringstream ss;
      ss << std::left << std::setw(18) << name << " " << info.substr(0,60);
      return ss.str();
    };

    if(!description.empty()) sinexTextReference->readLine(refLine("DESCRIPTION", description));
    if(!output.empty())      sinexTextReference->readLine(refLine("OUTPUT",      output));
    if(!contact.empty())     sinexTextReference->readLine(refLine("CONTACT",     contact));
    if(!software.empty())    sinexTextReference->readLine(refLine("SOFTWARE",    software));
    if(!hardware.empty())    sinexTextReference->readLine(refLine("HARDWARE",    hardware));
    if(!input.empty())       sinexTextReference->readLine(refLine("INPUT",       input));

    // comment block
    SinexTextPtr sinexTextComment = addBlock<SinexText>("FILE/COMMENT");
    if(!fileNameComment.empty())
    {
      InFile commentFile(fileNameComment);
      std::string line;
      while(std::getline(commentFile, line))
        sinexTextComment->readLine(line);
    }
    for(const auto& line : comment)
      sinexTextComment->readLine(line);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
std::string Sinex::header() const
{
  try
  {
    UInt countParameter = 0;
    if(hasBlock("SOLUTION/APRIORI"))
      countParameter = getBlock<SinexSolutionVector>("SOLUTION/APRIORI")->size();
    std::string header = _header;
    if(header.size() > 65)
      header.replace(60, 5, countParameter%"%05i"s);
    return header;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void Sinex::readFile(const FileName &fileName)
{
  try
  {
    if(fileName.empty())
      throw(Exception("File name is empty."));

    InFile file(fileName);

    std::string line, blockLabel;
    BlockType blockType = UNKNOWN;
    if(file.peek() == '%')
      std::getline(file, _header);
    else
      logWarning << "mandatory header in first line (starting with %) is missing" << Log::endl;

    while(std::getline(file, line))
    {
      // trim from end
      line.erase(std::find_if(line.rbegin(), line.rend(), [](int ch) { return !std::isspace(ch); }).base(), line.end());

      // skip comments
      if(line.empty() || (line.at(0) == '*'))
        continue;

      // %ENDSNX
      if(line.at(0) == '%')
        break;

      // start data block
      if(line.at(0) == '+')
      {
        if(!blockLabel.empty())
          throw(Exception("New SINEX block starts unexpectedly: '" + line + "' in block '" + blockLabel + "'"));
        blockLabel = String::trim(line.substr(1));
        blockType  = Sinex::blockType(blockLabel);
        addBlock<SinexBlock>(blockLabel);
        continue;
      }

      // end data block
      if(line.at(0) == '-')
      {
        if(blockLabel != String::trim(line.substr(1)))
          throw(Exception("SINEX block ends unexpectedly: '" + line + "' in block '" + blockLabel + "'"));
        blockLabel.clear();
        blockType = UNKNOWN;
        continue;
      }

      // unknown line
      if(line.at(0) != ' ')
      {
        if(blockType != FILE_COMMENT)
          logWarning << "Unknown line identifier: '" << line << "'" << Log::endl;
        continue;
      }

      if(blockType == FILE_COMMENT)
        continue;

      // data lines
      if(_blocks.find(blockType) == _blocks.end())
        continue;
      _blocks.at(blockType)->readLine(line);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void Sinex::writeFile(const FileName &fileName) const
{
  try
  {
    if(fileName.empty())
      throw(Exception("File name is empty."));

    OutFile file(fileName);
    file << header() << std::endl;
    file << "*" << std::string(79, '-') << std::endl;
    for(const auto& block : _blocks)
    {
      if(block.second->writeBlock(file))
        file << "*" << std::string(79, '-') << std::endl;
    }
    file << "%ENDSNX";
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

std::string Sinex::time2str(Time time, Bool fourDigitYear)
{
  try
  {
    if(time == Time() || time == date2time(2500,1,1,0,0,0.))
      return (fourDigitYear ? "0000" : "00") + ":000:00000"s;

    // round to full second including rollover
    UInt   year, month, day, hour, minute;
    Double second;
    time.date(year, month, day, hour, minute, second);
    time = date2time(year, month, day, hour, minute, std::round(second)+0.1);
    time.date(year, month, day, hour, minute, second);

    std::stringstream ss;
    if(fourDigitYear)
      ss << year%"%04i"s << ":";
    else
      ss << (year%100)%"%02i"s << ":";
    ss << time.dayOfYear()%"%03i"s << ":";
    ss << std::round(time.mjdMod()*86400)%"%05i"s;
    return ss.str();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Time Sinex::str2time(const std::string &line, std::size_t pos, Bool fourDigitYear)
{
  try
  {
    UInt posOffset = fourDigitYear ? 2 : 0;
    UInt year = static_cast<UInt>(String::toInt(line.substr(pos+0, 2+posOffset)));
    UInt day  = static_cast<UInt>(String::toInt(line.substr(pos+3+posOffset, 3)));
    UInt sec  = static_cast<UInt>(String::toInt(line.substr(pos+7+posOffset, 5)));
    Time time;
    if(year != 0 || day != 0 || sec != 0)
    {
      if(!fourDigitYear)
        year += (year <= 50) ? 2000 : 1900;
      time = date2time(year,1,1) + mjd2time(day-1.) + seconds2time(static_cast<Double>(sec));
    }
    if(year == 0 && day == 0 && sec == 0)
      time = date2time(2500,1,1);
    return time;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Sinex::SinexTextPtr Sinex::addSiteIdBlock(const std::vector<GnssStationInfo> &stationInfos)
{
  try
  {
    auto rad2DegMinSec = [] (Double rad, Double &deg, Double &min, Double& sec)
    {
      deg = std::floor(rad*RAD2DEG);
      min = std::floor((rad*RAD2DEG-deg)*60);
      sec = ((rad*RAD2DEG-deg)*60-min)*60;
    };

    Ellipsoid ellipsoid;
    SinexTextPtr block = addBlock<Sinex::SinexText>("SITE/ID");
    for(const auto &stationInfo : stationInfos)
    {
      Angle lon, lat;
      Double height, lonDeg, lonMin, lonSec, latDeg, latMin, latSec;
      ellipsoid(stationInfo.approxPosition, lon, lat, height);
      rad2DegMinSec(lon < 0 ? lon+2*PI : lon, lonDeg, lonMin, lonSec);
      rad2DegMinSec(lat, latDeg, latMin, latSec);

      std::stringstream ss;
      ss << String::upperCase(resize(stationInfo.markerName, 4)) << "  A " << resize(stationInfo.markerNumber, 9) << " P " << resize(stationInfo.comment, 22) << " "
         << lonDeg%"%3i "s << lonMin%"%2i "s << lonSec%"%4.1f "s << latDeg%"%3i "s << latMin%"%2i "s << latSec%"%4.1f "s << height%"%7.1f"s;
      block->readLine(ss.str());
    }

    return block;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Sinex::SinexTextPtr Sinex::addSiteReceiverBlock(const std::vector<GnssStationInfo> &stationInfos, const Time &timeRef)
{
  try
  {
    SinexTextPtr block = addBlock<Sinex::SinexText>("SITE/RECEIVER");
    for(const auto &stationInfo : stationInfos)
    {
      const UInt idRecv = stationInfo.findReceiver(timeRef);
      if(idRecv == NULLINDEX)
      {
        logWarning << stationInfo.markerName << ": no receiver found at " << timeRef.dateTimeStr() << Log::endl;
        continue;
      }
      const GnssReceiverInfo recv = stationInfo.receiver.at(idRecv);
      std::stringstream ss;
      ss << String::upperCase(resize(stationInfo.markerName, 4)) << "  A    1 P " << Sinex::time2str(recv.timeStart) << " " << Sinex::time2str(recv.timeEnd) << " "
         << resize(recv.name, 20) << " " << (recv.serial.empty() ? "-----" : resize(recv.serial, 5)) << " " << (recv.version.empty() ? std::string(11, '-') : resize(recv.version, 11));
      block->readLine(ss.str());
    }

    return block;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Sinex::SinexTextPtr Sinex::addSiteAntennaBlock(const std::vector<GnssStationInfo> &stationInfos, const Time &timeRef)
{
  try
  {
    SinexTextPtr block = addBlock<Sinex::SinexText>("SITE/ANTENNA");
    for(const auto &stationInfo : stationInfos)
    {
      const UInt idAnt = stationInfo.findAntenna(timeRef);
      if(idAnt == NULLINDEX)
      {
        logWarning << stationInfo.markerName << ": no antenna found at " << timeRef.dateTimeStr() << Log::endl;
        continue;
      }
      const GnssAntennaInfo ant = stationInfo.antenna.at(idAnt);
      Angle roll, pitch, yaw;
      Rotary3d(ant.local2antennaFrame.matrix()).cardan(roll, pitch, yaw);
      std::stringstream ss;
      ss << String::upperCase(resize(stationInfo.markerName, 4)) << "  A    1 P " << Sinex::time2str(ant.timeStart) << " " << Sinex::time2str(ant.timeEnd) << " "
         << resize(ant.name, 15) << " " << (ant.radome.empty() ? "NONE" : resize(ant.radome, 4)) << " " << (ant.serial.empty() ? "-----" : resize(ant.serial, 5))
         << " " << (yaw*RAD2DEG)%"% 4i"s;
      block->readLine(ss.str());
    }

    return block;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Sinex::SinexTextPtr Sinex::addSiteGpsPhaseCenterBlock(const std::vector<GnssStationInfo> &stationInfos, const Time &timeRef, const std::string &antennaModel)
{
  try
  {
    SinexTextPtr block = addBlock<Sinex::SinexText>("SITE/GPS_PHASE_CENTER");
    for(const auto &stationInfo : stationInfos)
    {
      const UInt idAnt = stationInfo.findAntenna(timeRef);
      if(idAnt == NULLINDEX)
      {
        logWarning << stationInfo.markerName << ": no antenna found at " << timeRef.dateTimeStr() << Log::endl;
        continue;
      }
      const GnssAntennaInfo ant = stationInfo.antenna.at(idAnt);
      auto iterPattern1 = std::find_if(ant.antennaDef->pattern.begin(), ant.antennaDef->pattern.end(), [](const GnssAntennaPattern &pattern){ return pattern.type == GnssType::L1_G; });
      auto iterPattern2 = std::find_if(ant.antennaDef->pattern.begin(), ant.antennaDef->pattern.end(), [](const GnssAntennaPattern &pattern){ return pattern.type == GnssType::L2_G; });
      if(iterPattern1 == ant.antennaDef->pattern.end() || iterPattern2 == ant.antennaDef->pattern.end())
      {
        logWarning << stationInfo.markerName <<": GPS phase center not found for antenna " << ant.name << " " << ant.radome << " " << ant.serial << Log::endl;
        continue;
      }
      std::stringstream ss;
      ss << resize(ant.name, 15) << " " << (ant.radome.empty() ? "NONE" : resize(ant.radome, 4)) << " " << (ant.serial.empty() ? "-----" : resize(ant.serial, 5)) << " "
         << format(iterPattern1->offset.z()) << " " << format(iterPattern1->offset.x()) << " " << format(iterPattern1->offset.y()) << " "
         << format(iterPattern2->offset.z()) << " " << format(iterPattern2->offset.x()) << " " << format(iterPattern2->offset.y()) << " " << resize(antennaModel, 10);
      block->readLine(ss.str());
    }

    return block;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Sinex::SinexTextPtr Sinex::addSiteGalileoPhaseCenterBlock(const std::vector<GnssStationInfo> &stationInfos, const Time &timeRef, const std::string &antennaModel)
{
  try
  {
    SinexTextPtr block = addBlock<Sinex::SinexText>("SITE/GAL_PHASE_CENTER");
    for(const auto &stationInfo : stationInfos)
    {
      const UInt idAnt = stationInfo.findAntenna(timeRef);
      if(idAnt == NULLINDEX)
      {
        logWarning << stationInfo.markerName << ": no antenna found at " << timeRef.dateTimeStr() << Log::endl;
        continue;
      }
      const GnssAntennaInfo ant = stationInfo.antenna.at(idAnt);
      auto iterPattern1 = std::find_if(ant.antennaDef->pattern.begin(), ant.antennaDef->pattern.end(), [](const GnssAntennaPattern &pattern){ return pattern.type == GnssType::L1_E; });
      auto iterPattern5 = std::find_if(ant.antennaDef->pattern.begin(), ant.antennaDef->pattern.end(), [](const GnssAntennaPattern &pattern){ return pattern.type == GnssType::L5_E; });
      auto iterPattern6 = std::find_if(ant.antennaDef->pattern.begin(), ant.antennaDef->pattern.end(), [](const GnssAntennaPattern &pattern){ return pattern.type == GnssType::L6_E; });
      auto iterPattern7 = std::find_if(ant.antennaDef->pattern.begin(), ant.antennaDef->pattern.end(), [](const GnssAntennaPattern &pattern){ return pattern.type == GnssType::L7_E; });
      auto iterPattern8 = std::find_if(ant.antennaDef->pattern.begin(), ant.antennaDef->pattern.end(), [](const GnssAntennaPattern &pattern){ return pattern.type == GnssType::L8_E; });
      if(iterPattern1 == ant.antennaDef->pattern.end() || iterPattern5 == ant.antennaDef->pattern.end() || iterPattern6 == ant.antennaDef->pattern.end() ||
         iterPattern7 == ant.antennaDef->pattern.end() || iterPattern8 == ant.antennaDef->pattern.end())
      {
//        logWarning << station <<": Galileo phase center not found for antenna " << ant.name << " " << ant.radome << " " << ant.serial << Log::endl;
        continue;
      }
      std::stringstream ss;
      ss << resize(ant.name, 15) << " " << (ant.radome.empty() ? "NONE" : resize(ant.radome, 4)) << " " << (ant.serial.empty() ? "-----" : resize(ant.serial, 5)) << " "
         << format(iterPattern1->offset.z()) << " " << format(iterPattern1->offset.x()) << " " << format(iterPattern1->offset.y()) << " "
         << format(iterPattern5->offset.z()) << " " << format(iterPattern5->offset.x()) << " " << format(iterPattern5->offset.y()) << " " << resize(antennaModel, 10) << " ";
      block->readLine(ss.str());
      ss.str("");
      ss << resize(ant.name, 15) << " " << (ant.radome.empty() ? "NONE" : resize(ant.radome, 4)) << " " << (ant.serial.empty() ? "-----" : resize(ant.serial, 5)) << " "
         << format(iterPattern6->offset.z()) << " " << format(iterPattern6->offset.x()) << " " << format(iterPattern6->offset.y()) << " "
         << format(iterPattern7->offset.z()) << " " << format(iterPattern7->offset.x()) << " " << format(iterPattern7->offset.y()) << " " << resize(antennaModel, 10) << " ";
      block->readLine(ss.str());
      ss.str("");
      ss << resize(ant.name, 15) << " " << (ant.radome.empty() ? "NONE" : resize(ant.radome, 4)) << " " << (ant.serial.empty() ? "-----" : resize(ant.serial, 5)) << " "
         << format(iterPattern8->offset.z()) << " " << format(iterPattern8->offset.x()) << " " << format(iterPattern8->offset.y()) << " "
         << std::string(20, ' ') << " " << resize(antennaModel, 10) << " ";
      block->readLine(ss.str());
    }

    return block;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Sinex::SinexTextPtr Sinex::addSiteEccentricityBlock(const std::vector<GnssStationInfo> &stationInfos, const Time &timeRef)
{
  try
  {
    SinexTextPtr block = addBlock<Sinex::SinexText>("SITE/ECCENTRICITY");
    for(const auto &stationInfo : stationInfos)
    {
      const UInt idAnt = stationInfo.findAntenna(timeRef);
      if(idAnt == NULLINDEX)
      {
        logWarning << stationInfo.markerName << ": no antenna found at " << timeRef.dateTimeStr() << Log::endl;
        continue;
      }
      const GnssAntennaInfo ant = stationInfo.antenna.at(idAnt);
      std::stringstream ss;
      ss << String::upperCase(resize(stationInfo.markerName, 4)) << "  A    1 P " << Sinex::time2str(ant.timeStart) << " "
         << Sinex::time2str(ant.timeEnd) << " UNE " << ant.position.z()%"%8.4f "s << ant.position.x()%"%8.4f "s << ant.position.y()%"%8.4f "s;
      block->readLine(ss.str());
    }

    return block;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Sinex::SinexTextPtr Sinex::addSatelliteIdBlock(const std::vector<GnssStationInfo> &transmitterInfos, const Time &timeRef)
{
  try
  {
    SinexTextPtr block = addBlock<Sinex::SinexText>("SATELLITE/ID");
    std::map<std::string, Time> svn2TimeStart, svn2TimeEnd;
    for(const auto &transmitterInfo : transmitterInfos)
      for(const auto &ant : transmitterInfo.antenna)
      {
        svn2TimeStart[ant.serial] = (svn2TimeStart[ant.serial] == Time() ? ant.timeStart : std::min(svn2TimeStart[ant.serial], ant.timeStart));
        svn2TimeEnd[ant.serial]   = std::max(svn2TimeStart[ant.serial], ant.timeEnd);
      }

    for(const auto &transmitterInfo : transmitterInfos)
    {
      const UInt idAnt = transmitterInfo.findAntenna(timeRef);
      if(idAnt == NULLINDEX)
      {
        logWarning << transmitterInfo.markerNumber << ": no antenna found at " << timeRef.dateTimeStr() << Log::endl;
        continue;
      }
      GnssAntennaInfo ant = transmitterInfo.antenna.at(idAnt);
      std::stringstream ss;
      ss << ant.serial << " " << transmitterInfo.markerNumber.substr(1,2) << " " << ant.radome << " P " << Sinex::time2str(svn2TimeStart[ant.serial])
         << " " << Sinex::time2str(svn2TimeEnd[ant.serial]) << " " << ant.name;
      block->readLine(ss.str());
    }

    return block;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Sinex::SinexTextPtr Sinex::addSatellitePhaseCenter(const std::vector<GnssStationInfo> &transmitterInfos, const Time &timeRef, const std::string &antennaModel)
{
  try
  {
    SinexTextPtr block = addBlock<Sinex::SinexText>("SATELLITE/PHASE_CENTER");
    auto antennaOffsetStr = [&](const GnssAntennaInfo &info, const GnssAntennaPattern &pattern)
    {
      Vector3d offset = info.position + info.local2antennaFrame.inverseTransform(pattern.offset);
      return pattern.type.str().substr(1,1) + " " + format(offset.z()) + " " + format(offset.x()) + " " + format(offset.y());
    };
    for(const auto &transmitterInfo : transmitterInfos)
    {
      const UInt idAnt = transmitterInfo.findAntenna(timeRef);
      if(idAnt == NULLINDEX)
      {
        logWarning << transmitterInfo.markerNumber << ": no antenna found at " << timeRef.dateTimeStr() << Log::endl;
        continue;
      }
      const GnssAntennaInfo ant = transmitterInfo.antenna.at(idAnt);
      std::vector<const GnssAntennaPattern*> patterns;
      for(const auto &pattern : ant.antennaDef->pattern)
        if(pattern.type == GnssType::PHASE)
          patterns.push_back(&pattern);
      std::sort(patterns.begin(), patterns.end(), [](const auto &p1, const auto &p2){ return p1->type < p2->type; });
      for(UInt i = 0; i < patterns.size(); i += 2)
      {
        std::stringstream ss;
        ss << ant.serial << " " << antennaOffsetStr(ant, *patterns.at(i)) << " " << (i+1 < patterns.size() ? antennaOffsetStr(ant, *patterns.at(i+1)) : std::string(22, ' '))
           << " " << resize(antennaModel, 10) << " A " << (patterns.at(i)->pattern.rows() > 1 ? "F" : "E");
        block->readLine(ss.str());
      }
    }

    return block;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}


/***********************************************/

std::string Sinex::SinexText::header() const
{
  try
  {
    if(blockType() == FILE_REFERENCE)
      return "*INFO_TYPE_________ INFO________________________________________________________";
    if(blockType() == INPUT_ACKNOWLEDGEMENTS)
      return "*AGY ______________________________FULL_DESCRIPTION_____________________________";
    if(blockType() == INPUT_HISTORY)
      return "*_VERSION_ CRE __CREATION__ OWN _DATA_START_ __DATA_END__ T PARAM S ____TYPE____";
    if(blockType() == INPUT_FILES)
      return "*OWN __CREATION__ ___________FILENAME__________ ___________DESCRIPTION__________";
    if(blockType() == SITE_ID)
      return "*SITE PT __DOMES__ T _STATION DESCRIPTION__ _LONGITUDE_ _LATITUDE__ HEIGHT_";
    if(blockType() == SITE_RECEIVER)
      return "*SITE PT SOLN T _DATA START_ __DATA_END__ ___RECEIVER_TYPE____ _S/N_ _FIRMWARE__";
    if(blockType() == SITE_ANTENNA)
      return "*SITE PT SOLN T _DATA START_ __DATA_END__ ____ANTENNA_TYPE____ _S/N_";
    if(blockType() == SITE_GPS_PHASE_CENTER)
      return "*____ANTENNA_TYPE____ _S/N_ _L1_U_ _L1_N_ _L1_E_ _L2_U_ _L2_N_ _L2_E_ _ANTMODEL_";
    if(blockType() == SITE_GAL_PHASE_CENTER)
      return "*____ANTENNA_TYPE____ _S/N_ _L1_U_ _L1_N_ _L1_E_ _L5_U_ _L5_N_ _L5_E_ _ANTMODEL_\n"s +
             "*____ANTENNA_TYPE____ _S/N_ _L6_U_ _L6_N_ _L6_E_ _L7_U_ _L7_N_ _L7_E_ _ANTMODEL_\n"s +
             "*____ANTENNA_TYPE____ _S/N_ _L8_U_ _L8_N_ _L8_E_ ____________________ _ANTMODEL_"s;
    if(blockType() == SITE_ECCENTRICITY)
      return "*SITE PT SOLN T _DATA START_ __DATA_END__ AXE _ECC_U__ _ECC_N__ _ECC_E__";
    if(blockType() == SOLUTION_EPOCHS)
      return "*SITE PT SOLN T _DATA_START_ __DATA_END__ _MEAN_EPOCH_";
    if(blockType() == SATELLITE_ID)
      return "*SVN_ PR COSPAR_ID T _DATA_START_ __DATA_END__ ______ANTENNA_______";
    if(blockType() == SATELLITE_PHASE_CENTER)
      return "*SVN_ L SATA_Z SATA_X SATA_Y L SATA_Z SATA_X SATA_Y _ANTMODEL_ T M";
    if(blockType() == BIAS_DESCRIPTION)
      return "*KEYWORD________________________________ VALUE(S)_______________________________";
    if(blockType() == BIAS_SOLUTION)
      return "*BIAS SVN_ PRN STATION__ OBS1 OBS2 BIAS_START____ BIAS_END______ UNIT __ESTIMATED_VALUE____ _STD_DEV___ __ESTIMATED_SLOPE____ _STD_DEV___";
    return "";
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}


/***********************************************/

Bool Sinex::SinexText::writeBlock(std::ostream &file) const
{
  try
  {
    if(!file.good())
      throw(Exception("Cannot write SINEX block: '" + label() + "'"));

    if(!_lines.size())
      return FALSE;

    file << "+" << label() << std::endl;
    if(!header().empty())
      file << header() << std::endl;
    for(const auto& line : _lines)
      if(line.length()>0)
        file << (line.at(0) != ' ' ? " " : "") << line.substr(0,(line.at(0) != ' ' ? 79 : 80)) << std::endl;
      else
        file<<std::endl;
    file << "-" << label() << std::endl;

    return TRUE;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void Sinex::SinexSolutionEpochs::readLine(const std::string &line)
{
  try
  {
    Epoch epoch;
    epoch.siteCode        = String::trim(line.substr(1, 4));
    epoch.pointCode       = String::trim(line.substr(6, 2));
    epoch.solutionId      = String::trim(line.substr(9, 4));
    epoch.observationCode = String::trim(line.substr(14, 1));
    epoch.timeStart       = str2time(line, 16);
    epoch.timeEnd         = str2time(line, 29);
    epoch.timeMean        = str2time(line, 42);
    _epochs.push_back(epoch);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Bool Sinex::SinexSolutionEpochs::writeBlock(std::ostream &file) const
{
  try
  {
    if(!file.good())
      throw(Exception("Cannot write SINEX block: '" + label() + "'"));

    if(!_epochs.size())
      return FALSE;

    file << "+" << label() << std::endl;
    file << header() << std::endl;
    for(const auto& epoch : _epochs)
    {
      file << " " << std::left  << std::setw(4)  << (epoch.siteCode.empty()        ? "----" : epoch.siteCode);
      file << " " << std::right << std::setw(2)  << (epoch.pointCode.empty()       ? "--"   : epoch.pointCode);
      file << " " << std::right << std::setw(4)  << (epoch.solutionId.empty()      ? "----" : epoch.solutionId);
      file << " " << std::right << std::setw(1)  << (epoch.observationCode.empty() ? "-"    : epoch.observationCode);
      file << " " << std::right << std::setw(12) << time2str(epoch.timeStart);
      file << " " << std::right << std::setw(12) << time2str(epoch.timeEnd);
      file << " " << std::right << std::setw(12) << time2str(epoch.timeMean);
      file << std::endl;
    }
    file << "-" << label() << std::endl;

    return TRUE;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

std::string Sinex::SinexSolutionVector::header() const
{
  try
  {
    if(blockType() == SOLUTION_APRIORI)
      return "*INDEX _TYPE_ CODE PT SOLN _REF_EPOCH__ UNIT S ____APRIORI_VALUE____ __STD_DEV__";
    if(blockType() == SOLUTION_ESTIMATE)
      return "*INDEX _TYPE_ CODE PT SOLN _REF_EPOCH__ UNIT S ___ESTIMATED_VALUE___ __STD_DEV__";
    if(blockType() == SOLUTION_NORMAL_EQUATION_VECTOR)
      return "*INDEX _TYPE_ CODE PT SOLN _REF_EPOCH__ UNIT S ___RIGHT_HAND_SIDE___";
    return "*";
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void Sinex::SinexSolutionVector::readLine(const std::string &line)
{
  try
  {
    Parameter parameter;
    parameter.parameterIndex = static_cast<UInt>(String::toInt(line.substr(1, 5)));
    parameter.parameterType  = String::trim(line.substr(7, 6));
    parameter.siteCode       = String::trim(line.substr(14, 4));
    parameter.pointCode      = String::trim(line.substr(19, 2));
    parameter.solutionId     = String::trim(line.substr(22, 4));
    parameter.time           = str2time(line, 27);
    parameter.unit           = String::trim(line.substr(40, 4));
    parameter.constraintCode = String::trim(line.substr(45, 1));
    parameter.value          = String::toDouble(line.substr(47, 21));
    if(blockType() != SOLUTION_NORMAL_EQUATION_VECTOR && line.length() > 68)
      parameter.sigma        = String::toDouble(line.substr(69, 11));
    _parameters.resize(parameter.parameterIndex);
    _parameters.back() = parameter;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Bool Sinex::SinexSolutionVector::writeBlock(std::ostream &file) const
{
  try
  {
    if(!file.good())
      throw(Exception("Cannot write SINEX block: '" + label() + "'"));

    if(!_parameters.size())
      return FALSE;

    file << "+" << label() << std::endl;
    file << header() << std::endl;
    for(const auto& param : _parameters)
    {
      file << " " << std::right << std::setw(5)  << param.parameterIndex;
      file << " " << std::left  << std::setw(6)  << param.parameterType;
      file << " " << std::left  << std::setw(4)  << (param.siteCode.empty()   ? "----" : param.siteCode);
      file << " " << std::right << std::setw(2)  << (param.pointCode.empty()  ? "--"   : param.pointCode);
      file << " " << std::right << std::setw(4)  << (param.solutionId.empty() ? "----" : param.solutionId);
      file << " " << std::right << std::setw(12) << time2str(param.time);
      file << " " << std::left  << std::setw(4)  << param.unit;
      file << " " << std::right << std::setw(1)  << param.constraintCode;
      file << " " << std::right << std::setw(21) << std::setprecision(14) << std::scientific << param.value;
      if(blockType() != SOLUTION_NORMAL_EQUATION_VECTOR)
        file << " " << std::right << std::setw(11) << std::setprecision(4)  << std::scientific << param.sigma;
      file << std::endl;
    }
    file << "-" << label() << std::endl;

    return TRUE;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Vector Sinex::SinexSolutionVector::vector() const
{
  try
  {
    Vector vector(size());
    for(UInt i = 0; i < size(); i++)
      vector(i) = _parameters.at(i).value;
    return vector;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

std::vector<ParameterName> Sinex::SinexSolutionVector::parameterNames() const
{
  try
  {
    std::vector<ParameterName> parameterNames(size());
    for(UInt i = 0; i < size(); i++)
    {
      // spherical harmonics coefficients
      if(_parameters.at(i).parameterType == "CN" || _parameters.at(i).parameterType == "SN")
      {
        std::string type = "sphericalHarmonics." + std::string(_parameters.at(i).parameterType == "CN" ? "c_" : "s_")
                         + String::trim(_parameters.at(i).siteCode) + "_" // degree
                         + String::trim(_parameters.at(i).solutionId);     // order
        parameterNames.at(i) = ParameterName("", type);
      }
      // other parameters
      else
      {
        std::string object = (!_parameters.at(i).siteCode.empty() && _parameters.at(i).siteCode != "----") ? String::trim(_parameters.at(i).siteCode) : "";
        if(_parameters.at(i).parameterType.substr(0,3) != "SAT")
          std::transform(object.begin(), object.end(), object.begin(), ::tolower);

        std::string type;
        if(     _parameters.at(i).parameterType == "STAX")   type = "position.x";
        else if(_parameters.at(i).parameterType == "STAY")   type = "position.y";
        else if(_parameters.at(i).parameterType == "STAZ")   type = "position.z";
        else if(_parameters.at(i).parameterType == "VELX")   type = "velocity.x";
        else if(_parameters.at(i).parameterType == "VELY")   type = "velocity.y";
        else if(_parameters.at(i).parameterType == "VELZ")   type = "velocity.z";
        else if(_parameters.at(i).parameterType == "XGC")    type = "geocenter.x";
        else if(_parameters.at(i).parameterType == "YGC")    type = "geocenter.y";
        else if(_parameters.at(i).parameterType == "ZGC")    type = "geocenter.z";
        else if(_parameters.at(i).parameterType == "LOD")    type = "LOD";
        else if(_parameters.at(i).parameterType == "UT")     type = "UT1";
        else if(_parameters.at(i).parameterType == "XPO")    type = "polarMotion.xp";
        else if(_parameters.at(i).parameterType == "YPO")    type = "polarMotion.yp";
        else if(_parameters.at(i).parameterType == "XPOR")   type = "polarMotionRate.xp";
        else if(_parameters.at(i).parameterType == "YPOR")   type = "polarMotionRate.yp";
        else if(_parameters.at(i).parameterType == "NUT_X")  type = "nutation.X";
        else if(_parameters.at(i).parameterType == "NUT_Y")  type = "nutation.Y";
        else if(_parameters.at(i).parameterType == "NUTR_X") type = "nutationRate.X";
        else if(_parameters.at(i).parameterType == "NUTR_Y") type = "nutationRate.Y";
        else if(_parameters.at(i).parameterType == "SAT__X") type = "position.x";
        else if(_parameters.at(i).parameterType == "SAT__Y") type = "position.y";
        else if(_parameters.at(i).parameterType == "SAT__Z") type = "position.z";
        else if(_parameters.at(i).parameterType == "SAT_VX") type = "velocity.x";
        else if(_parameters.at(i).parameterType == "SAT_VY") type = "velocity.y";
        else if(_parameters.at(i).parameterType == "SAT_VZ") type = "velocity.z";
        else if(_parameters.at(i).parameterType == "SATA_X") type = "antennaCenterVariations.xOffset";
        else if(_parameters.at(i).parameterType == "SATA_Y") type = "antennaCenterVariations.yOffset";
        else if(_parameters.at(i).parameterType == "SATA_Z") type = "antennaCenterVariations.zOffset";
        else                                                 type = _parameters.at(i).parameterType; // not all types implemented yet, see SINEX documentation
        parameterNames.at(i) = ParameterName(object, type);
      }
    }
    return parameterNames;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Sinex::SinexSolutionMatrix::SinexSolutionMatrix(const std::string &label, const UInt size) : SinexBlock(label.substr(0, label.find_first_of(' ')))
{
  try
  {
    if(size == 0)
      throw(Exception("Matrix size must be greater than zero: " + label));

    const std::string baseLabel  = label.substr(0, label.find_first_of(' '));
    const std::string triangular = label.length() > baseLabel.length() ? label.substr(baseLabel.length()+1, 1) : "";
    if(triangular == "L")
      _matrix = Matrix(size, Matrix::SYMMETRIC, Matrix::LOWER);
    else if(triangular == "U")
      _matrix = Matrix(size, Matrix::SYMMETRIC, Matrix::UPPER);
    else
      throw(Exception("Undefined SINEX matrix type for: " + label));

    if(baseLabel != "SOLUTION/NORMAL_EQUATION_MATRIX")
    {
      const std::string type = label.substr(baseLabel.length()+3);
      if(type == "CORR")
        _type = CORRELATION;
      else if(type == "COVA")
        _type = COVARIANCE;
      else if(type == "INFO")
        _type = INFORMATION;
      else
        throw(Exception("Undefined SINEX matrix type for: " + label));
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

std::string Sinex::SinexSolutionMatrix::label() const
{
  try
  {
    std::string label = SinexBlock::label();

    if(_matrix.getType() == Matrix::SYMMETRIC && _matrix.isUpper())
      label += " U";
    else if(_matrix.getType() == Matrix::SYMMETRIC && !_matrix.isUpper())
      label += " L";
    else
      throw(Exception("undefined SINEX matrix type for: " + label));

    if(label.substr(0, label.length()-2) != "SOLUTION/NORMAL_EQUATION_MATRIX")
    {
      if(type() == CORRELATION)
        label += " CORR";
      else if(type() == COVARIANCE)
        label += " COVA";
      else if(type() == INFORMATION)
        label += " INFO";
      else
        throw(Exception("undefined SINEX matrix type for: " + label));
    }

    return label;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void Sinex::SinexSolutionMatrix::readLine(const std::string &line)
{
  try
  {
    UInt i = static_cast<UInt>(String::toInt(line.substr(1, 5))) - 1;
    UInt j = static_cast<UInt>(String::toInt(line.substr(7, 5))) - 1;
    for(UInt k = 0; k < 3; k++)
      if(line.length() >= 13+k*22+21)
        _matrix(i,j+k) += String::toDouble(line.substr(13+k*22, 21));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Bool Sinex::SinexSolutionMatrix::writeBlock(std::ostream &file) const
{
  try
  {
    if(!file.good())
      throw(Exception("Cannot write SINEX block: '" + label() + "'"));

    if(!_matrix.size())
      return FALSE;

    file << "+" << label() << std::endl;
    file << header() << std::endl;
    const UInt size    = _matrix.rows();
    const Bool isUpper = _matrix.isUpper();
    for(UInt i = 0; i < size; i++)
      for(UInt j = (isUpper ? i : 0); j < (isUpper ? size : i+1); j++)
        if(_matrix(i,j) != 0.)
        {
          file << " " << std::setw(5) << i+1 << " " << std::setw(5) << j+1;
          for(UInt k = 0; k < 3; k++)
          {
            if(j < (isUpper ? size : i+1) && _matrix(i,j) != 0.)
              file << " " << std::right << std::setw(21) << std::setprecision(14) << _matrix(i, (k<2) ? j++ : j);
            else
              break;
          }
          file << std::endl;
        }
    file << "-" << label() << std::endl;

    return TRUE;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void Sinex::SinexSolutionStatistics::readLine(const std::string &line)
{
  try
  {
    _values[String::trim(line.substr(1, 30))] = String::toDouble(line.substr(32, 22));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Bool Sinex::SinexSolutionStatistics::writeBlock(std::ostream &file) const
{
  try
  {
    if(!file.good())
      throw(Exception("Cannot write SINEX block: '" + label() + "'"));

    if(!_values.size())
      return FALSE;

    file << "+" << label() << std::endl;
    file << header() << std::endl;
    for(const auto& value : _values)
      file << " " << std::left << std::setw(30) << value.first << " " << (std::ceil(value.second) == value.second ? value.second%"%22i"s : value.second%"%22.15e"s) << std::endl;
    file << "-" << label() << std::endl;

    return TRUE;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Double Sinex::SinexSolutionStatistics::value(const std::string &name) const
{
  try
  {
    auto iter = _values.find(name);
    if(iter == _values.end())
      throw(Exception("SINEX solution statistics not found: " + name));

    return _values.at(name);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
