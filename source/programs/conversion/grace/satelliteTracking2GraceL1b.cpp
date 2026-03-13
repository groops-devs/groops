/***********************************************/
/**
* @file satelliteTracking2GraceL1b.cpp
*
* @brief Read GROOPS data and convert to GRACE L1B data.
*
* @author Felix Oehlinger
* @date 2026-03-03
*
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program converts low-low satellite tracking data (KBR or LRI) from
the GROOPS format \file{instrument file (SATELLITETRACKING)}{instrument}
to the GRACE SDS format (KBR1B or LRI1B).
It reads the satellite tracking data and optionally corrections
(antenna offsets and light time corrections) and flags into one \config{outputfile}.

The text file \config{inputfileHeader} is placed at the beginning of the \config{outputfile}.
The \reference{text parser}{general.parser:text} is applied so that all variables can be used.
In addition, the times of the data are available with the variables \verb|{epochmin}|, \verb|{epochmax}|,
and \verb|{epochcount}|.

See also \program{GraceL1b2SatelliteTracking}.
)";

/***********************************************/

#include "programs/program.h"
#include "parser/stringParser.h"
#include "parser/dataVariables.h"
#include "files/fileInstrument.h"

/***** CLASS ***********************************/

/** @brief Read GROOPS L1B data.
* @ingroup programsConversionGroup */
class SatelliteTracking2GraceL1b
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(SatelliteTracking2GraceL1b, SINGLEPROCESS, "read GROOPS L1B data (KBR1B or LRI1B)", Conversion, Grace, Instrument)

/*************************************************************/

void SatelliteTracking2GraceL1b::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName    fileNameOut;
    FileName    fileNameHeader, fileNameSst, fileNameAntCentr, fileNameLighttime, fileNameSnr, fileNameIonoCorr;
    std::string satelliteId;

    readConfig(config, "outputfile",                 fileNameOut,       Config::MUSTSET,  "", "KBR1B or LRI1B");
    readConfig(config, "inputfileHeader",            fileNameHeader,    Config::OPTIONAL, "", "YAML Header, {epochmin}, {epochmax}, {epochcount} available");
    readConfig(config, "inputfileSatelliteTracking", fileNameSst,       Config::MUSTSET,  "", "SATELLITETRACKING");
    readConfig(config, "inputfileIonoCorr",          fileNameIonoCorr,  Config::OPTIONAL, "", "MISCVALUE");
    readConfig(config, "inputfileLighttime",         fileNameLighttime, Config::OPTIONAL, "", "SATELLITETRACKING");
    readConfig(config, "inputfileAntCentr",          fileNameAntCentr,  Config::OPTIONAL, "", "SATELLITETRACKING");
    readConfig(config, "inputfileSNR",               fileNameSnr,       Config::OPTIONAL, "", "MISCVALUES(K_A_SNR, Ka_A_SNR, K_B_SNR, Ka_B_SNR, qualflg)");
    if(isCreateSchema(config)) return;

    // =============================================

    logStatus<<"write SDS output <"<<fileNameOut<<">"<<Log::endl;
    SatelliteTrackingArc arcSst      = InstrumentFile::read(fileNameSst);
    SatelliteTrackingArc arcLight    = InstrumentFile::read(fileNameLighttime);
    SatelliteTrackingArc arcAntCentr = InstrumentFile::read(fileNameAntCentr);
    MiscValueArc         arcIonoCorr = InstrumentFile::read(fileNameIonoCorr);
    MiscValuesArc        arcSnr      = InstrumentFile::read(fileNameSnr);
    Arc::checkSynchronized({arcSst, arcIonoCorr, arcLight, arcAntCentr, arcSnr});

    OutFile file(fileNameOut);
    if(!fileNameHeader.empty())
    {
      VariableList varList = config.getVarList();
      addDataVariables("epoch", arcSst.times(), varList);

      InFile fileHeader(fileNameHeader);
      std::string line;
      while(std::getline(fileHeader, line))
        file<<StringParser::parse(line, varList)<<std::endl;
    }
    else
      file<<"# End of YAML header"<<std::endl;

    Single::forEach(arcSst.size(), [&](UInt i)
    {
      SatelliteTrackingEpoch sst, light, antCentr;
      MiscValueEpoch         ionoCorr;
      MiscValuesEpoch        snr(5);

      sst = arcSst.at(i);
      if(arcLight   .size()) light    = arcLight.at(i);
      if(arcAntCentr.size()) antCentr = arcAntCentr.at(i);
      if(arcIonoCorr.size()) ionoCorr = arcIonoCorr.at(i);
      if(arcSnr     .size()) snr      = arcSnr.at(i);

      file<<static_cast<UInt>(std::round((arcSst.at(i).time-mjd2time(51544.5)).seconds()));
      file<<sst.range%" %.16g"s<<sst.rangeRate%" %.16g"s<<sst.rangeAcceleration%" %.16g"s;
      file<<ionoCorr.value%" %.16g"s;
      file<<light.range%" %.16g"s<<light.rangeRate%" %.16g"s<<light.rangeAcceleration%" %.16g"s;
      file<<antCentr.range%" %.16g"s<<antCentr.rangeRate%" %.16g"s<<antCentr.rangeAcceleration%" %.16g"s;
      file<<snr.values.at(0)%" %i"s<<snr.values.at(1)%" %i"s<<snr.values.at(2)%" %i"s<<snr.values.at(3)%" %i  "s;
      const Byte qualflg = static_cast<Byte>(snr.values.at(4));
      for(UInt k=8; k-->0;)
        file<<(((qualflg>>k) & 1) ? '1' : '0');
      file<<std::endl;
    });
  }
  catch(std::exception &e)
  {
      GROOPS_RETHROW(e)
  }
}

/***********************************************/
