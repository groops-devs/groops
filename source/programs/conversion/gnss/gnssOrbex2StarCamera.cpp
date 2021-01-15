/***********************************************/
/**
* @file gnssOrbex2StarCamera.cpp
*
* @brief Converts satellite attitude from ORBEX to StarCamera format.
*
* @author Sebastian Strasser
* @date 2020-09-29
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
Converts GNSS satellite attitude from \href{http://acc.igs.org/misc/proposal_orbex_april2019.pdf}{ORBEX file format}
(quaternions) to \file{instrument file (STARCAMERA)}{instrument}.
The resulting star camera files contain the rotation from satellite body frame to TRF, or to CRF in case
\configClass{earthRotation}{earthRotationType} is provided.

See also \program{GnssAttitude2Orbex}.
)";

/***********************************************/

#include "programs/program.h"
#include "base/string.h"
#include "classes/earthRotation/earthRotation.h"
#include "files/fileInstrument.h"

/***** CLASS ***********************************/

/** @brief Converts satellite attitude from ORBEX to StarCamera format.
* @ingroup programsGroup */
class GnssOrbex2StarCamera
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(GnssOrbex2StarCamera, SINGLEPROCESS, "Converts satellite attitude from ORBEX to StarCamera format.", Conversion, Gnss, Instrument)

/***********************************************/

void GnssOrbex2StarCamera::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName fileNameOutStarCamera, fileNameInOrbex;
    std::vector<std::string> identifiers;
    EarthRotationPtr         earthRotation;

    readConfig(config, "outputfileStarCamera", fileNameOutStarCamera, Config::MUSTSET,  "",  "rotation from body frame to TRF/CRF, identifier is appended to each file");
    readConfig(config, "inputfileOrbex",       fileNameInOrbex,       Config::MUSTSET,  "",  "");
    readConfig(config, "identifier",           identifiers,           Config::OPTIONAL, "",  "(empty = all) satellite identifier, e.g. G23 or E05");
    readConfig(config, "earthRotation",        earthRotation,         Config::OPTIONAL, "",  "rotation from TRF to CRF");
    if(isCreateSchema(config)) return;

    logStatus<<"read ORBEX file <"<<fileNameInOrbex<<">"<<Log::endl;
    InFile file(fileNameInOrbex);
    std::map<std::string, StarCameraArc> starCameraArcs;
    std::string line;
    Bool isDataBlock = FALSE;
    Time time;
    Rotary3d crf2trf;
    while(std::getline(file, line))
    {
      if(String::startsWith(line, "+EPHEMERIS/DATA"))
      {
        isDataBlock = TRUE;
        continue;
      }
      if(String::startsWith(line, "-EPHEMERIS/DATA"))
        isDataBlock = FALSE;
      if(!isDataBlock || String::startsWith(line, "*"))
        continue;

      if(String::startsWith(line, "##"))
      {
        time = date2time(String::toInt(line.substr(3,4)), String::toInt(line.substr(8,2)), String::toInt(line.substr(11,2)),
                         String::toInt(line.substr(14,2)), String::toInt(line.substr(17,2)), String::toDouble(line.substr(20,15)));
        if(earthRotation)
          crf2trf = earthRotation->rotaryMatrix(time);
      }

      if(String::startsWith(line, " ATT"))
      {
        std::string id = line.substr(5,3);
        if(identifiers.size() && std::find(identifiers.begin(), identifiers.end(), id) == identifiers.end())
          continue;

        Rotary3d trf2sat(Vector({String::toDouble(line.substr(24,19)), String::toDouble(line.substr(44,19)),
                             String::toDouble(line.substr(64,19)), String::toDouble(line.substr(84,19))}));
        StarCameraEpoch epoch;
        epoch.time = time;
        epoch.rotary = inverse(trf2sat * crf2trf); // (trf/crf2sat --> sat2trf/crf)
        starCameraArcs[id].push_back(epoch);
      }
    }

    logStatus<<"write star camera files <"<<fileNameOutStarCamera<<">"<<Log::endl;
    for(const auto &sat : starCameraArcs)
      InstrumentFile::write(fileNameOutStarCamera.appendBaseName("."+sat.first), sat.second);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
