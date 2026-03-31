/***********************************************/
/**
* @file starCamera2GraceL1b.cpp
*
* @brief Read GROOPS L1B data and convert to SDS file format.
*
* @author Felix Öhlinger
* @date 2026-03-06
*
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program converts orientation data measured by a star camera (SRF to CRF)
from the GROOPS format \file{instrument file (STARCAMERA)}{instrument} to the GRACE SDS format (SCA1B).

It reads one \config{inputfileStarCamera} and optionally one
\config{inputfileStarCameraFlags} containing \file{MISCVALUES}{instrument}(sca\_id, qual\_rss, qualflg),
and writes one SDS output file.

The text file \config{inputfileHeader} is placed at the beginning of the \config{outputfile}.
The \reference{text parser}{general.parser:text} is applied so that all variables can be used.
In addition, the times of the data are available with the variables \verb|{epochmin}|, |\verb|{epochmax}|,
and \verb|{epochcount}|.

See also \program{GraceL1b2StarCamera}.
)";

/***********************************************/

#include "programs/program.h"
#include "parser/stringParser.h"
#include "parser/dataVariables.h"
#include "files/fileInstrument.h"

/***** CLASS ***********************************/

/** @brief Read GROOPS L1B data and convert to SDS data format.
* @ingroup programsConversionGroup */
class StarCamera2GraceL1b
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(StarCamera2GraceL1b, SINGLEPROCESS, "convert GROOPS L1B data to SDS data (SCA1B)", Conversion, Grace, Instrument)

/*************************************************************/

void StarCamera2GraceL1b::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName    fileNameOut;
    FileName    fileNameHeader, fileNameSca, fileNameFlags;
    std::string satelliteId;

    readConfig(config, "outputfile",               fileNameOut,    Config::MUSTSET,  "", "SCA1B");
    readConfig(config, "inputfileHeader",          fileNameHeader, Config::OPTIONAL, "", "YAML Header, {epochmin}, {epochmax}, {epochcount} available");
    readConfig(config, "inputfileStarCamera",      fileNameSca,    Config::MUSTSET,  "", "STARCAMERA");
    readConfig(config, "inputfileStarCameraFlags", fileNameFlags,  Config::OPTIONAL, "", "MISCVALUES(sca_id, qual_rss, qualflg)");
    readConfig(config, "satelliteId",              satelliteId,    Config::MUSTSET,  "", "A, B, C or D");
    if(isCreateSchema(config)) return;

        // =============================================

    logStatus<<"write SDS output <"<<fileNameOut<<">"<<Log::endl;
    StarCameraArc arcSca   = InstrumentFile::read(fileNameSca);
    MiscValuesArc arcFlags = InstrumentFile::read(fileNameFlags);
    Arc::checkSynchronized({arcSca, arcFlags});

    OutFile file(fileNameOut);
    if(!fileNameHeader.empty())
    {
      VariableList varList = config.getVarList();
      addDataVariables("epoch", arcSca.times(), varList);

      InFile fileHeader(fileNameHeader);
      std::string line;
      while(std::getline(fileHeader, line))
        file<<StringParser::parse(line, varList)<<std::endl;
    }
    else
      file<<"# End of YAML header"<<std::endl;

    Single::forEach(arcSca.size(), [&](UInt i)
    {
      const Vector q = arcSca.at(i).rotary.quaternion();
      MiscValuesEpoch flags(3);
      if(arcFlags.size()) flags = arcFlags.at(i);

      file<<static_cast<UInt>(std::round((arcSca.at(i).time-mjd2time(51544.5)).seconds()))<<" "<<satelliteId;
      file<<flags.values.at(0)%" %i"s; // sca_id
      for(UInt k=0; k<4; k++)
        file<<q(k)%" %20.16g"s;
      file<<flags.values.at(1)%" %.16g "s; // qual_rss
      const Byte qualflg  = static_cast<Byte>(flags.values(2));
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
