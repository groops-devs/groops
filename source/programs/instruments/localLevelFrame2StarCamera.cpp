/***********************************************/
/**
* @file localLevelFrame2StarCamera.cpp
*
* @brief Compute rotation matrix from local level frame (north, east, down) to TRF and save as StarCamera file.
*
* @author Sebastian Strasser
* @date 2017-06-28
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
Compute rotation (\file{StarCamera file}{instrument}) from local level frame (ellipsoidal north, east, down)
to TRF for positions given in \configFile{inputfileInstrument}{instrument} (first 3 data columns).
)";

/***********************************************/

#include "programs/program.h"
#include "files/fileInstrument.h"

/***********************************************/

/** @brief Rotation from local level frame (ellipsoidal north, east, down) to TRF.
* @ingroup programsGroup */
class LocalLevelFrame2StarCamera
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(LocalLevelFrame2StarCamera, PARALLEL, "Rotation from local level frame (ellipsoidal north, east, down) to TRF.", Instrument, Simulation)
GROOPS_RENAMED_PROGRAM(SimulateLocalLevelFrame, LocalLevelFrame2StarCamera, date2time(2018, 7, 4))

/***********************************************/

void LocalLevelFrame2StarCamera::run(Config &config, Parallel::CommunicatorPtr comm)
{
  try
  {
    Double a, f;
    FileName fileNameStarCamera;
    FileName fileNameInstrument;
    Bool constantOrigin;

    readConfig(config, "outputfileStarCamera", fileNameStarCamera, Config::MUSTSET,  "",  "rotation matrix from local level frame (ellipsoidal north, east, down) to TRF");
    readConfig(config, "inputfileInstrument",  fileNameInstrument, Config::MUSTSET,  "",  "origin of local level frame");
    readConfig(config, "constantOriginPerArc", constantOrigin,     Config::DEFAULT,  "0", "use constant origin for all epochs of an arc (median position)");
    readConfig(config, "R",                    a,                  Config::DEFAULT,  STRING_DEFAULT_GRS80_a, "reference radius for ellipsoidal coordinates");
    readConfig(config, "inverseFlattening",    f,                  Config::DEFAULT,  STRING_DEFAULT_GRS80_f, "reference flattening for ellipsoidal coordinates, 0: spherical coordinates");
    if(isCreateSchema(config)) return;

    logStatus<<"compute local level frame rotation"<<Log::endl;
    Ellipsoid ellipsoid(a, f);
    InstrumentFile epochFile(fileNameInstrument);
    std::vector<Arc> arcList(epochFile.arcCount());
    Parallel::forEach(arcList, [&](UInt arcNo)
    {
      Arc      arc = epochFile.readArc(arcNo);
      Matrix   A   = arc.matrix();
      Vector3d origin(median(A.column(1)), median(A.column(2)), median(A.column(3)));

      StarCameraArc scaArc;
      for(UInt k=0; k<A.rows(); k++)
      {
        // Local (north, east, down) -> global TRF
        StarCameraEpoch epoch;
        epoch.time   = arc.at(k).time;
        epoch.rotary = localNorthEastDown((constantOrigin ? origin : Vector3d(A.slice(k,1,1,3))), ellipsoid);
        scaArc.push_back(epoch);
      }
      return scaArc;
    }, comm);

    if(Parallel::isMaster(comm))
    {
      logStatus<<"write star camera data to file <"<<fileNameStarCamera<<">"<<Log::endl;
      InstrumentFile::write(fileNameStarCamera, arcList);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
