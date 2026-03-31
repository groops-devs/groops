/***********************************************/
/**
* @file accelerometer2GraceL1b.cpp
*
* @brief Read GROOPS accelerometer data and write GRACE L1B data.
*
* @author Felix Oehlinger
* @date 2025-01-21
*
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program converts accelerometer data from the \file{instrument file (ACCELEROMETER)}{instrument}
format into GRACE SDS format.

The text file \config{inputfileHeader} is placed at the beginning of the \config{outputfile}.
The \reference{text parser}{general.parser:text} is applied so that all variables can be used.
In addition, the times of the data are available with the variables \verb|{epochmin}|, \verb|{epochmax}|,
and \verb|{epochcount}|.

See also \program{GraceL1b2Accelerometer}.
)";

/***********************************************/

#include "programs/program.h"
#include "parser/stringParser.h"
#include "parser/dataVariables.h"
#include "files/fileInstrument.h"

/***** CLASS ***********************************/

/** @brief Read GROOPS accelerometer data and write GRACE L1B data.
* @ingroup programsConversionGroup */
class Accelerometer2GraceL1b
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(Accelerometer2GraceL1b, SINGLEPROCESS, "convert GROOPS ACC data into SDS L1B data", Conversion, Grace, Instrument)

/*************************************************************/

void Accelerometer2GraceL1b::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName    fileNameOut;
    FileName    fileNameHeader, fileNameAcc, fileNameAng, fileNameFlags;
    std::string satelliteId;

    readConfig(config, "outputfile",                    fileNameOut,    Config::MUSTSET,  "", "ACT1B");
    readConfig(config, "inputfileHeader",               fileNameHeader, Config::OPTIONAL, "", "YAML Header, {epochmin}, {epochmax}, {epochcount} available");
    readConfig(config, "inputfileAccelerometer",        fileNameAcc,    Config::MUSTSET,  "", "ACCELEROMETER");
    readConfig(config, "inputfileAngularAccelerometer", fileNameAng,    Config::OPTIONAL, "", "ACCELEROMETER");
    readConfig(config, "inputfileFlags",                fileNameFlags,  Config::OPTIONAL, "", "MISCVALUES(qualflg, acl_res.x, acl_res.y, acl_res.z)");
    readConfig(config, "satelliteId",                   satelliteId,    Config::MUSTSET,  "", "A, B, C or D");
    if(isCreateSchema(config)) return;

    // =============================================

    logStatus<<"write SDS output <"<<fileNameOut<<">"<<Log::endl;
    AccelerometerArc accArc   = InstrumentFile::read(fileNameAcc);
    AccelerometerArc angArc   = InstrumentFile::read(fileNameAng);
    MiscValuesArc    flagsArc = InstrumentFile::read(fileNameFlags);
    Arc::checkSynchronized({accArc, angArc, flagsArc});

    OutFile file(fileNameOut);
    if(!fileNameHeader.empty())
    {
      VariableList varList = config.getVarList();
      addDataVariables("epoch", accArc.times(), varList);

      InFile fileHeader(fileNameHeader);
      std::string line;
      while(std::getline(fileHeader, line))
        file<<StringParser::parse(line, varList)<<std::endl;
    }
    else
      file<<"# End of YAML header"<<std::endl;

    Single::forEach(accArc.size(), [&](UInt i)
    {
      Vector3d acc, ang, res;
      Byte     qualflg  = 0;

      acc = accArc.at(i).acceleration;
      if(angArc.size())
        ang = angArc.at(i).acceleration;
      if(flagsArc.size())
      {
        qualflg  = static_cast<Byte>(flagsArc.at(i).values.at(0));
        res      = Vector3d(flagsArc.at(i).values.row(1, 3));
      }

      file<<static_cast<UInt>(std::round((accArc.at(i).time-mjd2time(51544.5)).seconds()))<<" "<<satelliteId;
      file<<acc.x()%" %.16g"s<<acc.y()%" %.16g"s<<acc.z()%" %.16g"s;
      file<<ang.x()%" %.16g"s<<ang.y()%" %.16g"s<<ang.z()%" %.16g"s;
      file<<res.x()%" %.16g"s<<res.y()%" %.16g"s<<res.z()%" %.16g  "s;
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
