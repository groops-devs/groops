/***********************************************/
/**
* @file variational2Orbit.cpp
*
* @brief Extracts orbit from variational file.
*
* @author Sebastian Strasser
* @date 2016-07-13
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
Extracts the reference \configFile{outputfileOrbit}{instrument}
from \configFile{inputfileVariational}{variationalEquation}.
)";

/***********************************************/

#include "programs/program.h"
#include "files/fileInstrument.h"
#include "files/fileVariationalEquation.h"

/***** CLASS ***********************************/

/** @brief Extracts orbit from variational file.
* @ingroup programsGroup */
class Variational2Orbit
{
 public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(Variational2Orbit, SINGLEPROCESS, "Extracts orbit from variational file.", Misc, VariationalEquation, Instrument)

/***********************************************/

void Variational2Orbit::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName fileNameOutOrbit, fileNameInVariational;

    renameDeprecatedConfig(config, "outputFileOrbit",      "outputfileOrbit",      date2time(2020, 7, 15));
    renameDeprecatedConfig(config, "inputFileVariational", "inputfileVariational", date2time(2020, 7, 15));

    readConfig(config, "outputfileOrbit",      fileNameOutOrbit,      Config::MUSTSET, "", "output orbit (instrument) file");
    readConfig(config, "inputfileVariational", fileNameInVariational, Config::MUSTSET, "", "input variational file");
    if(isCreateSchema(config)) return;

    FileVariationalEquation fileVariational(fileNameInVariational);
    std::list<Arc> arcList;
    for(UInt arcNo = 0; arcNo < fileVariational.arcCount(); arcNo++)
      arcList.push_back(fileVariational.readArc(arcNo).orbitArc());

    logStatus << "write orbit to file <" << fileNameOutOrbit << ">" << Log::endl;
    InstrumentFile::write(fileNameOutOrbit, arcList);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
