/***********************************************/
/**
* @file instrumentStarCamera2RotaryMatrix.cpp
*
* @brief Compute rotary matrix.
*
* @author Torsten Mayer-Guerr
* @date 2020-02-02
*
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
Write \configFile{inputfileStarCamera}{instrument} rotations
as \configFile{outputfileInstrument}{instrument} rotary matrices
(for each epoch $xx, xy, xz, yx, yy, yz, zx, zy, zz$).
)";

/***********************************************/

#include "programs/program.h"
#include "files/fileInstrument.h"

/***** CLASS ***********************************/

/** @brief Compute rotary matrix.
* @ingroup programsGroup */
class InstrumentStarCamera2RotaryMatrix
{
public:
  void run(Config &config);
};

GROOPS_REGISTER_PROGRAM(InstrumentStarCamera2RotaryMatrix, PARALLEL, "Compute rotary matrix", Instrument)

/***********************************************/

void InstrumentStarCamera2RotaryMatrix::run(Config &config)
{
  try
  {
    FileName fileNameOut, fileNameStarCamera;

    readConfig(config, "outputfileInstrument", fileNameOut,        Config::MUSTSET, "", "xx, xy, xz, yx, yy, yz, zx, zy, zz (MISCVALUES)");
    readConfig(config, "inputfileStarCamera",  fileNameStarCamera, Config::MUSTSET, "", "");
    if(isCreateSchema(config)) return;

    logStatus<<"read star camera data <"<<fileNameStarCamera<<">"<<Log::endl;
    InstrumentFile file(fileNameStarCamera);

    std::vector<Arc> arcList(file.arcCount());
    Parallel::forEach(arcList, [&](UInt arcNo)
    {
      StarCameraArc starCamera = file.readArc(arcNo);
      MiscValuesArc arc;
      for(UInt i=0; i<starCamera.size(); i++)
      {
        MiscValuesEpoch epoch(6);
        epoch.time   = starCamera.at(i).time;
        epoch.values = flatten(starCamera.at(i).rotary.matrix().trans());
        arc.push_back(epoch);
      }
      return arc;
    }); // forEach

    if(Parallel::isMaster())
    {
      logStatus<<"write rotary matrix to file <"<<fileNameOut<<">"<<Log::endl;
      InstrumentFile::write(fileNameOut, arcList);
      Arc::printStatistics(arcList);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
