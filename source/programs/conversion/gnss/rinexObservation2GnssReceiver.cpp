/***********************************************/
/**
* @file rinexObservation2GnssReceiver.cpp
*
* @brief Converts RINEX or Compact RINEX files to GROOPS GnssReceiver Instrument file.
*
* @author Sebastian Strasser
* @date 2018-07-30
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
Converts \href{https://files.igs.org/pub/data/format/rinex_4.00.pdf}{RINEX} (version 2, 3, and 4) and
\href{https://terras.gsi.go.jp/ja/crx2rnx.html}{Compact RINEX} observation files to
\file{GnssReceiver Instrument file}{instrument}.

In case of \href{https://files.igs.org/pub/data/format/rinex211.txt}{RINEX v2.x} observation files
containing GLONASS satellites, a mapping from PRN
to frequency number must be provided via \config{inputfileMatrixPrn2FrequencyNumber}
in the form of a \file{matrix file}{matrix} with columns: GLONASS PRN, mjdStart, mjdEnd, frequencyNumber.
Source for mapping: \url{http://semisys.gfz-potsdam.de/semisys/api/?symname=2002&format=json&satellite=GLO}.
RINEX v3+ observation files already contain this information.

\configClass{useType}{gnssType} and \configClass{ignoreType}{gnssType} can be used to filter
the observation types that will be exported.

If \configFile{inputfileStationInfo}{platform} is set, RINEX antenna and receiver info
will be cross-checked with the provided file and warnings are raised in case of differences.

A list of semi-codeless GPS receivers (observing C2D instead of C2W) can be provided via
\configFile{inputfileSemiCodelessReceivers}{stringList} in ASCII format with one receiver name per line.
Observation types will be automatically corrected for these receivers.

Some LEO satellites use special RINEX observation types, either from the unofficial RINEX v2.20
or custom ones. These can be provided via \configFile{inputfileSpecialObservationTypes}{stringTable}
in ASCII format. The file must  must contain a table with two columns, the first being the special type,
and the second being the equivalent RINEX v3 type.

%Example for RINEX v2.20:

%\begin{tabular}{ll}
%LA & L1C \\
%L1 & L1W \\
%L2 & L2W \\
%SA & S1C \\
%S1 & S1W \\
%S2 & S2W \\
%\end{tabular}
)";

/***********************************************/

#include "programs/program.h"
#include "base/string.h"
#include "inputOutput/file.h"
#include "files/filePlatform.h"
#include "files/fileMatrix.h"
#include "files/fileInstrument.h"
#include "files/fileStringTable.h"

/***** CLASS ***********************************/

/** @brief Converts RINEX or Compact RINEX files to GROOPS GnssReceiver Instrument file.
* @ingroup programsConversionGroup */
class RinexObservation2GnssReceiver
{
  struct FrequencyNumberInterval
  {
    FrequencyNumberInterval(const Time &timeStart_, const Time &timeEnd_, Int frequencyNumber_) : timeStart(timeStart_), timeEnd(timeEnd_), frequencyNumber(frequencyNumber_) {}

    Time timeStart, timeEnd;
    Int  frequencyNumber;
  };

  Double               rinexVersion, compactRinexVersion;
  Time                 timeOfFirstObs;
  Platform             stationInfo;
  Platform             stationInfoRinex;
  PlatformGnssAntenna  antennaRinex;
  PlatformGnssReceiver receiverRinex;
  GnssReceiverArc      receiverArc;

  Bool isSemiCodelessReceiver = FALSE;
  std::vector<std::string> semiCodelessReceivers;

  std::vector<GnssType> useType, ignoreType;
  std::map<std::string, std::vector<FrequencyNumberInterval>> prn2FrequencyNumbers;
  std::map<std::string, std::string> specialTypes;
  std::map<Char, std::vector<GnssType>> system2ObsTypes; // system, obsTypes

  void readHeader(InFile &file, UInt lineCount=MAX_UINT);
  void readObservationData(InFile &file);
  void readCompactObservationData(InFile &file);
  void checkStationInfo();
  void addEpoch(const Time &time, const std::vector<GnssType> &satNumber, const std::vector<Vector> &obs, Double clockError=0.);
  Bool getLine(InFile &file, std::string &line, std::string &label) const;
  Bool testLabel(const std::string &labelInLine, const std::string &label, Bool optional=TRUE) const;
  Time readEpochTime(const std::string &line) const;
  std::vector<GnssType> getSystemObsTypes(const GnssType &prn, const Time &time) const;
  Double readOptionalDouble(const std::string &line, size_t pos, size_t len);

  static Bool isGpsAntiSpoofingEnabled(const Time &time);

public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(RinexObservation2GnssReceiver, SINGLEPROCESS, "Converts RINEX or Compact RINEX files to GROOPS GnssReceiver Instrument file.", Conversion, Gnss, Instrument)

/***********************************************/

void RinexObservation2GnssReceiver::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName fileNameOut, fileNameInPrn2FrequencyNumber, fileNameStationInfo, fileNameSemiCodelessReceivers, fileNameSpecialObservationTypes;
    std::vector<FileName> fileNameInObs;

    readConfig(config, "outputfileGnssReceiver",             fileNameOut,                       Config::MUSTSET,  "", "");
    readConfig(config, "inputfileRinexObservation",          fileNameInObs,                     Config::MUSTSET,  "", "RINEX or Compact RINEX observation files");
    readConfig(config, "inputfileMatrixPrn2FrequencyNumber", fileNameInPrn2FrequencyNumber,     Config::OPTIONAL, "{groopsDataDir}/gnss/transmitter/glonassPrnSvn2FrequencyNumber.txt", "(required for RINEX v2 files containing GLONASS observations) GROOPS matrix with columns: GLONASS PRN, SVN, mjdStart, mjdEnd, frequencyNumber");
    readConfig(config, "inputfileStationInfo",               fileNameStationInfo,               Config::OPTIONAL, "", "used to determine semi-codeless receivers and to cross-check antenna and receiver info");
    readConfig(config, "inputfileSemiCodelessReceivers",     fileNameSemiCodelessReceivers,     Config::OPTIONAL, "{groopsDataDir}/gnss/receiverStation/semiCodelessReceivers.txt", "ASCII list with one receiver name per line");
    readConfig(config, "inputfileSpecialObservationTypes",   fileNameSpecialObservationTypes,   Config::OPTIONAL, "", "ASCII table mapping special observation types to RINEX 3 types, e.g.: LA L1C");
    readConfig(config, "useType",                            useType,                           Config::OPTIONAL, "", "only use observations that match any of these patterns");
    readConfig(config, "ignoreType",                         ignoreType,                        Config::OPTIONAL, "", "ignore observations that match any of these patterns");
    if(isCreateSchema(config)) return;

    if(!fileNameStationInfo.empty())
      readFilePlatform(fileNameStationInfo, stationInfo);

    if(!fileNameSemiCodelessReceivers.empty())
      readFileStringList(fileNameSemiCodelessReceivers, semiCodelessReceivers);

    if(!fileNameSpecialObservationTypes.empty())
    {
      std::vector<std::vector<std::string>> table;
      readFileStringTable(fileNameSpecialObservationTypes, table);
      for(const auto & row : table)
      {
        if(row.size() >= 2)
          specialTypes[row.at(0)] = row.at(1);
        else
          logWarning<<"invalid special type mapping for "<<row.at(0)<<Log::endl;
      }
    }

    if(!fileNameInPrn2FrequencyNumber.empty())
    {
      Matrix A;
      readFileMatrix(fileNameInPrn2FrequencyNumber, A);
      for(UInt i = 0; i < A.rows(); i++)
        prn2FrequencyNumbers[A(i,0)%"R%02i"s].push_back(FrequencyNumberInterval(mjd2time(A(i,2)), mjd2time(A(i,3)), static_cast<Int>(A(i,4))));
    }

    // read RINEX files
    for(const auto &fileName : fileNameInObs)
    {
      try
      {
      logStatus<<"read RINEX observation file <"<<fileName<<">"<<Log::endl;
      InFile file(fileName);
      file.exceptions(std::ios::badbit|std::ios::failbit);

      // read header
      std::string line, label;
      getLine(file, line, label);
      compactRinexVersion = 0;
      if(testLabel(label, "CRINEX VERS   / TYPE"))
      {
        compactRinexVersion = String::toDouble(line.substr(0, 20));
        getLine(file, line, label);
        testLabel(label, "CRINEX PROG / DATE", FALSE);
        getLine(file, line, label);
      }

      testLabel(label, "RINEX VERSION / TYPE", FALSE);
      rinexVersion = String::toDouble(line.substr(0, 9));
      if(rinexVersion<2)
        throw(Exception("Can only read RINEX files starting from RINEX version 2.0"));
      if(line.at(20)!='O')
        throw(Exception("File must contain observation data"));

      // clean up data from previous file
      system2ObsTypes.clear();
      stationInfoRinex = Platform();
      antennaRinex     = PlatformGnssAntenna();
      receiverRinex    = PlatformGnssReceiver();

      readHeader(file);
      checkStationInfo();
      if(compactRinexVersion > 0)
        readCompactObservationData(file);
      else
        readObservationData(file);
      }
      catch(std::exception &e)
      {
        logWarning<<e.what()<<"; continue..."<<Log::endl;
      }
    }

    if(receiverArc.size() == 0)
      throw(Exception("empty arc"));

    receiverArc.sort();
    receiverArc.removeDuplicateEpochs(/*keepFirst*/FALSE);

    logStatus<<"write gnss data to <"<<fileNameOut<<">"<<Log::endl;
    Arc::printStatistics(receiverArc);
    InstrumentFile::write(fileNameOut, receiverArc);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void RinexObservation2GnssReceiver::readHeader(InFile &file, UInt lineCount)
{
  try
  {
    std::string line, label;
    UInt unknownLineCount = 0;

    for(UInt idLine=0; idLine<lineCount; idLine++)
    {
      if(!getLine(file, line, label))
        throw(Exception("error while reading RINEX header"));
      if(std::all_of(line.begin(), line.end(), isspace))
        continue;
      // ====================================
      if(testLabel(label, "END OF HEADER"))
        break;
      // ====================================
      else if(testLabel(label, "PGM / RUN BY / DATE"))
      {
      }
      // ====================================
      else if(testLabel(label, "COMMENT"))
      {
      }
      // ====================================
      else if(testLabel(label, "MARKER NAME"))
      {
        std::string markerName = String::trim(line.substr(0,4)); //trim(line.substr(0,60));
        std::transform(markerName.begin(), markerName.end(), markerName.begin(), ::toupper);
        if((!stationInfoRinex.markerName.empty()) && (stationInfoRinex.markerName != markerName))
          logWarning<<"Marker name changed '"<<markerName<<"' != '"<<stationInfoRinex.markerName<<"'"<<Log::endl;
        stationInfoRinex.markerName = markerName;
      }
      // ====================================
      else if(testLabel(label, "MARKER TYPE"))
      {
        std::string comment = String::trim(line.substr(0,60));
        if((!stationInfoRinex.comment.empty()) && (stationInfoRinex.comment != comment))
          logWarning<<"Marker type changed '"<<comment<<"' != '"<<stationInfoRinex.comment<<"'"<<Log::endl;
        stationInfoRinex.comment = comment;
      }
      // ====================================
      else if(testLabel(label, "MARKER NUMBER"))
      {
        std::string markerNumber = String::trim(line.substr(0,20));
        if((!stationInfoRinex.markerNumber.empty()) && (stationInfoRinex.markerNumber != markerNumber))
          logWarning<<"Marker number changed '"<<markerNumber<<"' != '"<<stationInfoRinex.markerNumber<<"'"<<Log::endl;
        stationInfoRinex.markerNumber = markerNumber;
      }
      // ====================================
      else if(testLabel(label, "OBSERVER / AGENCY"))
      {
      }
      // ====================================
      else if(testLabel(label, "REC # / TYPE / VERS"))
      {
        receiverRinex.serial  = String::trim(line.substr(0,20));
        receiverRinex.name    = String::trim(line.substr(20,20));
        receiverRinex.version = String::trim(line.substr(40,20));
        std::transform(receiverRinex.name.begin(), receiverRinex.name.end(), receiverRinex.name.begin(), ::toupper); // convert to upper case
      }
      // ====================================
      else if(testLabel(label, "ANT # / TYPE"))
      {
        antennaRinex.serial = String::trim(line.substr(0,20));
        antennaRinex.name   = String::trim(line.substr(20,16));
        antennaRinex.radome = String::trim(line.substr(36,4));
        if(antennaRinex.radome == "NONE")
          antennaRinex.radome.clear();
      }
      // ====================================
      else if(testLabel(label, "APPROX POSITION XYZ"))
      {
        stationInfoRinex.approxPosition.x() = readOptionalDouble(line, 0,14);
        stationInfoRinex.approxPosition.y() = readOptionalDouble(line,14,14);
        stationInfoRinex.approxPosition.z() = readOptionalDouble(line,28,14);
      }
      // ====================================
      else if(testLabel(label, "ANTENNA: DELTA H/E/N"))
      {
        antennaRinex.position.z() =  readOptionalDouble(line, 0,14);
        antennaRinex.position.y() =  readOptionalDouble(line,14,14);
        antennaRinex.position.x() =  readOptionalDouble(line,28,14);
      }
      // ====================================
      else if(testLabel(label, "ANTENNA: DELTA X/Y/Z"))
      {
        antennaRinex.position.x() = readOptionalDouble(line, 0,14);
        antennaRinex.position.y() = readOptionalDouble(line,14,14);
        antennaRinex.position.z() = readOptionalDouble(line,28,14);
      }
      // ====================================
      else if(testLabel(label, "ANTENNA: B.SIGHT XYZ"))
      {
      }
      // ====================================
      else if(testLabel(label, "ANTENNA: ZERODIR AZI "))
      {
      }
      // ====================================
      else if(testLabel(label, "ANTENNA: ZERODIR XYZ"))
      {
      }
      // ====================================
      else if(testLabel(label, "CENTER OF MASS: XYZ"))
      {
      }
      // ====================================
      else if(testLabel(label, "WAVELENGTH FACT L1/2"))
      {
        Double factorL1 = String::toInt(line.substr(0, 6));
        Double factorL2 = String::toInt(line.substr(6, 6));
        if((factorL1!=1)||(factorL2!=1))
        {
          logInfo<<"'"<<line<<"'"<<Log::endl;
          throw(Exception("not full cycle ambiguities"));
        }
      }
      // ====================================
      else if(testLabel(label, "# / TYPES OF OBSERV") || // version 2
              testLabel(label, "SYS / # / OBS TYPES"))   // version 3
      {
        const Char system = line[0] != ' ' ? line[0] : '*';
        const Int typeCount = String::toInt(line.substr(1, 5));

        std::stringstream ss(line.substr(7,53));
        for(Int i = 0; i < typeCount; i++)
        {
          std::string type;
          if(!(ss >> type)) // with possible continuation lines
          {
            getLine(file, line, label);
            idLine++;
            ss = std::stringstream(line.substr(7,53));
            ss >> type;
          }

          if(specialTypes.find(type) != specialTypes.end())
            type = specialTypes.at(type);

          if(rinexVersion < 3 && type.size() == 2)
          {
            if(type[0] == 'P')
              type = "C"s+type[1]+type[0]; // e.g. P1/P2 ==> C1P/C2P
            else
              type += '?'; // e.g. C1/L1/S2 ==> C1?/L1?/S2?
          }

          if(rinexVersion <= 3.02 && system == 'C' && type[1] == '1')
            type[1] = '2'; // version 3.02: BeiDou C1C/L1I/... ==> C2C/L2I/...

          system2ObsTypes[system].push_back(GnssType(type+system+"**"));
        }
      }
      // ====================================
      else if(testLabel(label, "INTERVAL"))
      {
        //interval = String::toDouble(line.substr(0, 10));
      }
      // ====================================
      else if(testLabel(label, "TIME OF FIRST OBS"))
      {
        Int year   = String::toInt(line.substr(0, 6));
        Int month  = String::toInt(line.substr(6, 6));
        Int day    = String::toInt(line.substr(12, 6));
        Int hour   = String::toInt(line.substr(18, 6));
        Int min    = String::toInt(line.substr(24, 6));
        Double sec = String::toDouble(line.substr(30, 13));
        timeOfFirstObs = date2time(year, month, day, hour, min, sec);
        if((line.substr(48,3)!="   ")&&(line.substr(48,3)!="GPS"))
          logWarning<<"not GPS time"<<Log::endl;
      }
      // ====================================
      else if(testLabel(label, "TIME OF LAST OBS"))
      {
      }
      // ====================================
      else if(testLabel(label, "RCV CLOCK OFFS APPL"))
      {
        Int flag = String::toInt(line.substr(0, 6));
        if(flag!=0)
          logWarning<<"RCV CLOCK OFFS APPL"<<Log::endl;
      }
      // ====================================
      else if(testLabel(label, "SYS / PHASE SHIFT"))
      {
        if(rinexVersion < 3.01 && line.substr(1,59).find_first_not_of(' ') != std::string::npos)
          logWarning<<"SYS / PHASE SHIFT contains corrections that may have not been applied"<<Log::endl;
      }
      // ====================================
      else if(testLabel(label, "GLONASS COD/PHS/BIS"))
      {
      }
      // ====================================
      else if(testLabel(label, "GLONASS SLOT / FRQ #"))
      {
        const Int satCount = String::toInt(line.substr(0, 3));

        std::stringstream ss(line.substr(4,56));
        for(Int i = 0; i < satCount; i++)
        {
          std::string prn;
          if(!(ss >> prn)) // with possible continuation lines
          {
            getLine(file, line, label);
            idLine++;
            ss = std::stringstream(line.substr(4,56));
            ss >> prn;
          }
          Int freqNo;
          ss >> freqNo;
          prn2FrequencyNumbers[prn] = {FrequencyNumberInterval(timeOfFirstObs, date2time(2500,1,1), freqNo)};
        }
      }
      // ====================================
      else if(testLabel(label, "LEAP SECONDS"))
      {
      }
      // ====================================
      else if(testLabel(label, "# OF SATELLITES"))
      {
      }
      // ====================================
      else if(testLabel(label, "PRN / # OF OBS"))
      {
      }
      // ====================================
      else if(testLabel(label, "CORR TO SYSTEM TIME"))
      {
      }
      // ====================================
      else if(testLabel(label, "SIGNAL STRENGTH UNIT"))
      {
      }
      // ====================================
      else if(testLabel(label, "DOI"))
      {
      }
      // ====================================
      else if(testLabel(label, "LICENCE OF USE"))
      {
      }
      // ====================================
      else if(testLabel(label, "STATION INFORMATION"))
      {
      }
      // ====================================
      else
      {
        if(++unknownLineCount < 20)
        {
          logWarning<<"Unknown header label:"<<Log::endl;
          logWarning<<"'"<<line<<"'"<<Log::endl;
        }
        else if(unknownLineCount == 20)
          logWarning<<"skipping additional warnings about unknown header labels"<<Log::endl;
      }
      // ====================================
    }

    // check if receiver measures semi-codeless (GPS L2)
    auto recv = stationInfo.findEquipment<PlatformGnssReceiver>(timeOfFirstObs);
    if(recv && std::find(semiCodelessReceivers.begin(), semiCodelessReceivers.end(), recv->name) != semiCodelessReceivers.end())
      isSemiCodelessReceiver = TRUE;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void RinexObservation2GnssReceiver::readObservationData(InFile &file)
{
  std::string line, label;

  try
  {
    while(getLine(file, line, label))
    {
      if(std::all_of(line.begin(), line.end(), isspace))
        continue;

      if(label.substr(0,7) == "COMMENT")
      {
        logWarning<<"ignoring comment line '"<<line<<"'"<<Log::endl;
        continue;
      }

      const Int  epochFlag = String::toInt(line.substr(rinexVersion < 3 ? 26 : 29, 3));
      const UInt satCount  = String::toInt(line.substr(rinexVersion < 3 ? 29 : 32, 3));

      // events?
      if((epochFlag>=2)&&(epochFlag!=6))
      {
        readHeader(file, satCount);
        continue;
      }

      const Time time = readEpochTime(line);
      const Double clockOffset = String::toDouble(line.substr(rinexVersion < 3 ? 68 : 41, rinexVersion < 3 ? 12 :15));

      // read observed satellites
      std::vector<GnssType> satNumber(satCount);
      if(rinexVersion < 3)
      {
        const UInt maxSatCountPerLine = 12;
        for(UInt idSat = 0; idSat < satCount; idSat++)
        {
          if(idSat > 0 && idSat%maxSatCountPerLine == 0) // with possible continuation lines
            getLine(file, line, label);
          std::string prn = line.substr(32+3*(idSat%maxSatCountPerLine), 3);
          if(prn.at(0) == ' ') prn.at(0) = 'G';
          if(prn.at(1) == ' ') prn.at(1) = '0';
          satNumber.at(idSat) = GnssType("***"+prn);
        }
      }

      // read observations
      const UInt maxObsCountPerLine = rinexVersion < 3 ? 5 : MAX_UINT;
      std::vector<Vector> obs(satCount);
      for(UInt idSat = 0; idSat < satCount; idSat++)
      {
        getLine(file, line, label);

        if(rinexVersion >= 3)
          satNumber.at(idSat) = GnssType("***"+line.substr(0,3));

        const UInt obsCount = getSystemObsTypes(satNumber.at(idSat), time).size();
        obs.at(idSat) = Vector(obsCount);
        if(rinexVersion >= 3)
          line.resize(3+16*obsCount, ' ');
        for(UInt idType = 0; idType < obsCount; idType++)
        {
          if(idType > 0 && idType%maxObsCountPerLine == 0) // with possible continuation lines
            getLine(file, line, label);
          obs.at(idSat)(idType) = String::toDouble(line.substr((rinexVersion >= 3 ? 3 : 0)+16*(idType%maxObsCountPerLine), 14));

          // TODO: LLI and signal strength
        }
      }

      if(epochFlag==6)
      {
        logWarning<<"skipping cycle slip records at "+time.dateTimeStr()<<Log::endl;
        continue;
      }

      addEpoch(time, satNumber, obs, clockOffset);
    }
  }
  catch(std::exception &e)
  {
    logError<<"'"<<line<<"'"<<Log::endl;
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

// Compact RINEX reference: Hatanaka, Y. (2008) A Compression Format and Tools for GNSS Observation Data, Bulletin of the Geographical Survey Institute, 55, 21-30
void RinexObservation2GnssReceiver::readCompactObservationData(InFile &file)
{
  std::string line, label;

  try
  {
    std::string epochLine;
    std::vector<Double> clockOffSetSeries;
    std::map<GnssType, std::vector<std::vector<Double>>> prn2TypeObsSeries; // prn, obsType, observation series
    while(getLine(file, line, label))
    {
      if(label.substr(0,7) == "COMMENT")
      {
        logWarning<<"ignoring comment line '"<<line<<"'"<<Log::endl;
        continue;
      }

      if(line[0] == '&' || line[0] == '>') // epoch line initialization
        epochLine = line;
      else // differential epoch line
      {
        epochLine.resize(std::max(line.size(), epochLine.size()), ' ');
        for(UInt i = 0; i < line.size(); i++)
        {
          if(line[i] == '&')
            epochLine[i] = ' ';
          else if(line[i] != ' ')
            epochLine[i] = line[i];
        }
        line = epochLine;
      }

      const Int  epochFlag = String::toInt(line.substr(rinexVersion < 3 ? 26 : 29, 3));
      const UInt satCount  = String::toInt(line.substr(rinexVersion < 3 ? 29 : 32, 3));

      // events?
      if((epochFlag>=2)&&(epochFlag!=6))
      {
        readHeader(file, satCount);
        continue;
      }

      Time time = readEpochTime(line);

      // read observed satellites
      std::vector<GnssType> satNumber(satCount);
      for(UInt idSat = 0; idSat < satCount; idSat++)
      {
        std::string prn = line.substr((rinexVersion < 3 ? 32 : 41)+3*(idSat), 3);
        if(prn.at(0) == ' ') prn.at(0) = 'G';
        if(prn.at(1) == ' ') prn.at(1) = '0';
        satNumber.at(idSat) = GnssType("***"+prn);
      }

      // read clock offset
      getLine(file, line, label);
      const Bool isEmpty = std::all_of(line.begin(), line.end(), isspace);
      if(!isEmpty && line[1] == '&') // initializing value
        clockOffSetSeries = {static_cast<Double>(line[0]-'0'), std::stod(line.substr(2))}; // {differential order, state at epoch}
      else if(!isEmpty) // differential value
      {
        if(clockOffSetSeries.size() < clockOffSetSeries.at(0)+2)
          clockOffSetSeries.push_back(std::stod(line)); // expand series by one if smaller than order+2
        std::partial_sum(clockOffSetSeries.rbegin(), clockOffSetSeries.rend()-1, clockOffSetSeries.rbegin()); // update series using cumulative sum
      }

      // read observations
      std::vector<Vector> obs(satCount);
      for(UInt idSat = 0; idSat < satCount; idSat++)
      {
        getLine(file, line, label);

        const UInt obsCount = getSystemObsTypes(satNumber.at(idSat), time).size();
        obs.at(idSat) = Vector(obsCount);

        prn2TypeObsSeries[satNumber.at(idSat)].resize(obsCount);

        std::stringstream ss(line);
        for(UInt idType = 0; idType < obsCount; idType++)
        {
          if(ss.get() == ' ' || !ss.good()) // check first character
            continue; // skip missing data

          const Bool isInitializing = (ss.get()== '&' ? TRUE : FALSE); // check second character
          if(ss.good())
            ss.unget(); // if second character existed, put it back into stream
          else
            ss.clear(); // else clear fail states
          ss.unget();   // put first character back into stream

          std::vector<Double> *obsSeries = &prn2TypeObsSeries[satNumber.at(idSat)].at(idType);
          if(isInitializing)
          {
            obsSeries->assign(2, 0);           // {differential order, state at epoch}
            obsSeries->at(0) = ss.get() - '0'; // set order
            ss.get();                          // skip the &
            ss >> obsSeries->at(1);            // set initial state
          }
          else // differential value
          {
            if(!obsSeries->size())
              throw(Exception("error in compact RINEX file (differential value given but arc is not initialized)"));
            if(obsSeries->size() < obsSeries->at(0)+2)
              obsSeries->push_back(0); // expand series by one if smaller than order+2
            ss >> obsSeries->back();
            std::partial_sum(obsSeries->rbegin(), obsSeries->rend()-1, obsSeries->rbegin()); // update series using cumulative sum
          }

          if(ss.good())
            ss.get(); // skip the trailing space

          obs.at(idSat).at(idType) = obsSeries->at(1)/1000;
        }

        // TODO: LLI and signal strength
      }

      if(epochFlag==6)
      {
        logWarning<<"skipping cycle slip records at "+time.dateTimeStr()<<Log::endl;
        continue;
      }

      addEpoch(time, satNumber, obs, clockOffSetSeries.size() ? clockOffSetSeries.at(1) : 0);
    }
  }
  catch(std::exception &e)
  {
    logError<<"'"<<line<<"'"<<Log::endl;
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void RinexObservation2GnssReceiver::checkStationInfo()
{
  try
  {
    if((!stationInfo.markerName.empty()) && (stationInfo.markerName != stationInfoRinex.markerName))
      logWarning<<timeOfFirstObs.dateTimeStr()<<" Marker name differs '"<<stationInfo.markerName<<"' != '"<<stationInfoRinex.markerName<<"'"<<Log::endl;
    if((!stationInfo.markerNumber.empty()) && (stationInfo.markerNumber != stationInfoRinex.markerNumber))
      logWarning<<timeOfFirstObs.dateTimeStr()<<" Marker number differs '"<<stationInfo.markerNumber<<"' != '"<<stationInfoRinex.markerNumber<<"'"<<Log::endl;

    // test antenna
    if(stationInfo.equipments.size())
    {
      auto antenna = stationInfo.findEquipment<PlatformGnssAntenna>(timeOfFirstObs);
      if(antenna)
      {
        if((!antennaRinex.name.empty()) && (antenna->name != antennaRinex.name))
          logWarning<<timeOfFirstObs.dateTimeStr()<<" Antenna name differs '"<<antenna->name<<"' != '"<<antennaRinex.name<<"'"<<Log::endl;
        /*
        if((!antennaRinex.serial.empty()) && (antenna->serial != antennaRinex.serial))
         */
        if((!antennaRinex.serial.empty()) && (antennaRinex.serial.compare(0,antenna->serial.length(),antenna->serial)))
          logWarning<<timeOfFirstObs.dateTimeStr()<<" Antenna serial differs '"<<antenna->serial<<"' != '"<<antennaRinex.serial<<"'"<<Log::endl;
        if((!antennaRinex.radome.empty()) && (antenna->radome != antennaRinex.radome))
          logWarning<<timeOfFirstObs.dateTimeStr()<<" Antenna radome differs '"<<antenna->radome<<"' != '"<<antennaRinex.radome<<"'"<<Log::endl;
        if((antenna->position-antennaRinex.position).r()>=0.001 || std::isnan(antennaRinex.position.r()))
          logWarning<<timeOfFirstObs.dateTimeStr()<<" Antenna delta position differs"<<Log::endl;

        antennaRinex = *antenna;
      }
      else
        logWarning << timeOfFirstObs.dateTimeStr() << ": Antenna not found in stationInfo: '"<<antennaRinex.name<<"'"<<Log::endl;
    }

    // test receiver
    if(stationInfo.equipments.size())
    {
      auto recv = stationInfo.findEquipment<PlatformGnssReceiver>(timeOfFirstObs);
      if(recv)
      {
        if((!receiverRinex.name.empty()) && (recv->name != receiverRinex.name))
          logWarning<<timeOfFirstObs.dateTimeStr()<<" Receiver name differs '"<<recv->name<<"' != '"<<receiverRinex.name<<"'"<<Log::endl;
        if((!receiverRinex.serial.empty()) && (recv->serial != receiverRinex.serial))
          logWarning<<timeOfFirstObs.dateTimeStr()<<" Receiver serial differs '"<<recv->serial<<"' != '"<<receiverRinex.serial<<"'"<<Log::endl;
        if((!receiverRinex.version.empty()) && (recv->version != receiverRinex.version))
          logWarning<<timeOfFirstObs.dateTimeStr()<<" Receiver version differs '"<<recv->version<<"' != '"<<receiverRinex.version<<"'"<<Log::endl;

        receiverRinex = *recv;
      }
      else
        logWarning << timeOfFirstObs.dateTimeStr() << ": Receiver not found in stationInfo: '"<<receiverRinex.name<<"'"<<Log::endl;
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void RinexObservation2GnssReceiver::addEpoch(const Time &time, const std::vector<GnssType> &satNumber, const std::vector<Vector> &obs, Double clockError)
{
  try
  {
    GnssReceiverEpoch epoch;
    epoch.time = time;
    epoch.clockError = clockError;

    // collect types and satellites that occur at this epoch
    const UInt satCount = satNumber.size();
    std::set<GnssType> occuringTypes;
    std::vector<Bool> isOccuringSat(satCount, FALSE);
    for(UInt idSat=0; idSat<satCount; idSat++)
    {
      const std::vector<GnssType> obsTypes = getSystemObsTypes(satNumber.at(idSat), time);
      for(UInt idType=0; idType<obsTypes.size(); idType++)
      {
        const GnssType type = obsTypes.at(idType) + satNumber.at(idSat);
        Bool use = !useType.size() ? TRUE : FALSE;
        if(GnssType::index(useType, type) != NULLINDEX)
          use = TRUE;
        if(GnssType::index(ignoreType, type) != NULLINDEX)
          use = FALSE;
        if(!use || obs.at(idSat)(idType) == 0.0)
          continue;
        occuringTypes.insert(type & GnssType::NOPRN);
        isOccuringSat.at(idSat) = TRUE;
      }
    }

    // add collected types to epoch
    for(const auto &type : occuringTypes)
      epoch.obsType.push_back(type);

    // add satellites and observations to epoch
    for(UInt idSat=0; idSat<satCount; idSat++)
      if(isOccuringSat.at(idSat))
      {
        GnssType sat = satNumber.at(idSat);

        // check GLONASS frequency number
        if((sat & GnssType::SYSTEM) == GnssType::GLONASS)
        {
          const std::string prn = sat.str().substr(3,3);
          try
          {
            auto iter = std::find_if(prn2FrequencyNumbers.at(prn).begin(), prn2FrequencyNumbers.at(prn).end(), [&](const auto &interval){ return time.isInInterval(interval.timeStart, interval.timeEnd); });
            if(iter != prn2FrequencyNumbers.at(prn).end())
              sat.setFrequencyNumber(iter->frequencyNumber);
            else
              throw(Exception("no frequency number found for "+prn+" at "+time.dateTimeStr()+" (frequency number file must be provided for RINEX v2.x)"));
          }
          catch(std::exception &e)
          {
            if(sat.PRN.prn() > 24) // ignore test satellites with PRN 26, 27 etc. if frequency number is unknown
              continue;
            GROOPS_RETHROW(e);
          }
        }

        epoch.satellite.push_back(sat);

        const std::vector<GnssType> obsTypes = getSystemObsTypes(sat, time);
        for(const auto &type : occuringTypes)
        {
          if((type & GnssType::SYSTEM) != (sat & GnssType::SYSTEM))
            continue;

          auto iter = std::find(obsTypes.begin(), obsTypes.end(), type);
          if(iter == obsTypes.end())
            continue;

          const UInt idType = std::distance(obsTypes.begin(), iter);
          epoch.observation.push_back(obs.at(idSat)(idType));
          if(type == GnssType::PHASE && obs.at(idSat)(idType) != 0.0)
            epoch.observation.back() *= (type+sat).wavelength(); // cycles -> meters
        }
      }

    // skip empty epochs
    if(epoch.satellite.size())
      receiverArc.push_back(epoch);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Bool RinexObservation2GnssReceiver::getLine(InFile &file, std::string &line, std::string &label) const
{
  try
  {
    std::getline(file, line);
    if(line.size() && line.back() == '\r')
      line.pop_back();
    if(line.size()<80)
      line.resize(80,' ');
    label = line.substr(60,20);
    return TRUE;
  }
  catch(...)
  {
    line.clear();
    line.resize(80,' ');
    label = line.substr(60,20);
    return FALSE;
  }
}

/***********************************************/

Bool RinexObservation2GnssReceiver::testLabel(const std::string &labelInLine, const std::string &label, Bool optional) const
{
  if(labelInLine.find(label)!=std::string::npos)
    return TRUE;
  if(optional)
    return FALSE;
  throw(Exception(std::string("In line '")+labelInLine+"' label '"+label+"' expected\n"));
}

/***********************************************/

Time RinexObservation2GnssReceiver::readEpochTime(const std::string &line) const
{
  try
  {
    Int year   = String::toInt(line.substr(rinexVersion < 3 ?  1 :  2, rinexVersion < 3 ? 2 : 4));
    Int month  = String::toInt(line.substr(rinexVersion < 3 ?  3 :  7, 3));
    Int day    = String::toInt(line.substr(rinexVersion < 3 ?  6 : 10, 3));
    Int hour   = String::toInt(line.substr(rinexVersion < 3 ?  9 : 13, 3));
    Int minute = String::toInt(line.substr(rinexVersion < 3 ? 12 : 16, 3));
    Double sec = String::toDouble(line.substr(rinexVersion < 3 ? 15 : 19, 11));
    if(rinexVersion < 3)
      year += (year<80 ? 2000 : 1900);
    return date2time(year, month, day, hour, minute, sec);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

std::vector<GnssType> RinexObservation2GnssReceiver::getSystemObsTypes(const GnssType &prn, const Time &time) const
{
  try
  {
    std::vector<GnssType> obsTypes = system2ObsTypes.at(rinexVersion < 3 ? '*' : prn.prnStr()[0]);

    // GPS: replace observation types
    if(prn == GnssType::GPS)
      for(auto &&type : obsTypes)
      {
        if(type == (GnssType::RANGE + GnssType::L1 + GnssType::UNKNOWN_ATTRIBUTE))
          type = GnssType::RANGE + GnssType::L1 + GnssType::C; // C1? ==> C1C
        if(type == GnssType::P)
          type = (type & GnssType::TYPE) + (type & GnssType::FREQUENCY) + GnssType::W; // **P ==> **W
        if(isSemiCodelessReceiver && type == (GnssType::RANGE + GnssType::L2 + GnssType::W) && isGpsAntiSpoofingEnabled(time))
          type = (type & GnssType::TYPE) + (type & GnssType::FREQUENCY) + GnssType::D; // C2W ==> C2D for semicodeless receivers if anti-spoofing is enabled
      }

    // GLONASS: replace observation types
    if(prn == GnssType::GLONASS)
      for(auto &&type : obsTypes)
        if(type == (GnssType::RANGE + GnssType::UNKNOWN_ATTRIBUTE) && (type == GnssType::G1 || type == GnssType::G2))
          type = GnssType::RANGE + (type & GnssType::FREQUENCY) + GnssType::C; // C1?/C2? ==> C1C/C2C

    return obsTypes;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

Double RinexObservation2GnssReceiver::readOptionalDouble(const std::string &line, size_t pos, size_t len)
{
  try
  {
    return String::toDouble(line.substr(pos, len));
  }
  catch(...)
  {
    return NAN_EXPR;
  }
}

/***********************************************/

Bool RinexObservation2GnssReceiver::isGpsAntiSpoofingEnabled(const Time &time)
{
  try
  {
    const Double mjd = time.mjd();

    if((                   mjd < 49383.)      ||
    (mjd >= 49826.87499 && mjd < 49847.83334) ||
    (mjd >= 49887.      && mjd < 49909.)      ||
    (mjd >= 50000.      && mjd < 50022.)      ||
    (mjd >= 50481.      && mjd < 50503.))
      return FALSE;

    return TRUE;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/
