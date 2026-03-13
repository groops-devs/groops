/***********************************************/
/**
* @file orbit2GraceL1b.cpp
*
* @brief Read GROOPS L1B data and convert to SDS file format.
*
* @author Torsten Mayer-Guerr
* @date 2026-03-12
*
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program converts an \file{instrument file (ORBIT)}{instrument} specified in the celestial reference frame (CRF)
to the GRACE/GRACE-FO SDS format (GNV1B, GNI1B). If \configClass{earthRotation}{earthRotationType} is provided,
the orbit is rotated into the terrestrial reference frame as required for the GNV1B product; otherwise,
a GNI1B product is written.

The text file \config{inputfileHeader} is placed at the beginning of the \config{outputfile}.
The \reference{text parser}{general.parser:text} is applied so that all variables can be used.
In addition, the times of the data are available with the variables \verb|{epochmin}|, \verb|{epochmax}|,
and \verb|{epochcount}|.

See also \program{GraceL1b2Orbit}.
)";

/***********************************************/

#include "programs/program.h"
#include "parser/stringParser.h"
#include "parser/dataVariables.h"
#include "files/fileInstrument.h"
#include "classes/earthRotation/earthRotation.h"

/***** CLASS ***********************************/

/** @brief Read GROOPS L1B data and convert to SDS data format.
* @ingroup programsConversionGroup */
class Orbit2GraceL1b
{
public:
    void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(Orbit2GraceL1b, SINGLEPROCESS, "convert GROOPS L1B data to SDS data", Conversion, Grace, Instrument)

/*************************************************************/

void Orbit2GraceL1b::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName         fileNameOut;
    FileName         fileNameHeader, fileNameOrbit;
    std::string      satelliteId;
    EarthRotationPtr earthRotation;

    readConfig(config, "outputfile",      fileNameOut,    Config::MUSTSET,  "",     "GNV1B/GNI1B");
    readConfig(config, "inputfileHeader", fileNameHeader, Config::OPTIONAL, "",     "YAML Header, {epochmin}, {epochmax}, {epochcount} available");
    readConfig(config, "inputfileOrbit",  fileNameOrbit,  Config::MUSTSET,  "",     "");
    readConfig(config, "satelliteId",     satelliteId,    Config::MUSTSET,  "",     "A, B, C or D");
    readConfig(config, "earthRotation",   earthRotation,  Config::OPTIONAL, "file", "rotation into Earth fixed frame");
    if(isCreateSchema(config)) return;

    // =============================================


    logStatus<<"write SDS output <"<<fileNameOut<<">"<<Log::endl;
    OrbitArc orbit = InstrumentFile::read(fileNameOrbit);
    OutFile file(fileNameOut);
    if(!fileNameHeader.empty())
    {
      VariableList varList = config.getVarList();
      addDataVariables("epoch", orbit.times(), varList);

      InFile fileHeader(fileNameHeader);
      std::string line;
      while(std::getline(fileHeader, line))
        file<<StringParser::parse(line, varList)<<std::endl;
    }
    else
      file<<"# End of YAML header"<<std::endl;

    Single::forEach(orbit.size(), [&](UInt i)
    {
      Vector3d p = orbit.at(i).position;
      Vector3d v = orbit.at(i).velocity;

      if(earthRotation)
      {
        const Rotary3d rotation = earthRotation->rotaryMatrix(orbit.at(i).time);
        if(v.r() > 0)
          v = rotation.rotate(v - crossProduct(earthRotation->rotaryAxis(orbit.at(i).time), p));
        p = rotation.rotate(p);
      }

      file<<static_cast<UInt>(std::round((orbit.at(i).time-mjd2time(51544.5)).seconds()))<<" "<<satelliteId<<" "<<((earthRotation) ? "E" : "I");
      file<<p.x()%" %18.16g"s<<p.y()%" %18.16g"s<<p.z()%" %18.16g"s;
      file<<" 1e+33 1e+33 1e+33";
      file<<v.x()%" %18.16g"s<<v.y()%" %18.16g"s<<v.z()%" %18.16g"s;
      file<<" 1e+33 1e+33 1e+33";
      file<<"  10000000"<<std::endl;
    });
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
