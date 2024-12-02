/***********************************************/
/**
* @file starCamera2Orbex.cpp
*
* @brief Converts satellite attitude from StarCamera to ORBEX format.
*
* @author Andre Hauschild
* @date 2024-10-10
*
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
Converts satellite attitude from \file{instrument file (STARCAMERA)}{instrument} to \href{http://acc.igs.org/misc/ORBEX009.pdf}{ORBEX file format}
(quaternions).

If \configClass{earthRotation}{earthRotationType} is provided, the output file contains quaternions for rotation from TRF to satellite
body frame (IGS/ORBEX convention), otherwise the rotation is from CRF to satellite body frame.

See also \program{GnssAttitude2Orbex}, \program{GnssOrbex2StarCamera}, \program{SimulateStarCameraGnss}.
)";

/***********************************************/

#include "programs/program.h"
#include "base/polynomial.h"
#include "classes/earthRotation/earthRotation.h"
#include "classes/timeSeries/timeSeries.h"
#include "files/fileInstrument.h"
#include "files/fileStringTable.h"
#include "inputOutput/system.h"

//#include "inputOutput/file.h"

/***** CLASS ***********************************/

/** @brief Convert satellite attitude from star camera instrument data to ORBEX file format (quaternions).
* @ingroup programsConversionGroup */
class starCamera2Orbex
{
  class Record
  {
    public:
      std::string identifier;
      Vector      quaternion;

      Record(const std::string &identifier, const_MatrixSliceRef quaternion) : identifier(identifier), quaternion(quaternion) {}
  };

public:
  class Satellite
  {
    public:
    FileName            fileNameStarCamera;
    std::string         identifier;
    StarCameraArc       starCamera;
    std::vector<Time>   times;
    UInt                idEpoch = 0;
  };

public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(starCamera2Orbex, SINGLEPROCESS, "Convert satellite attitude from star camera instrument file to ORBEX file format (quaternions).", Conversion, Instrument)

/***********************************************/

template<> Bool readConfig(Config &config, const std::string &name, starCamera2Orbex::Satellite &var, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation)
{
  if(!readConfigSequence(config, name, mustSet, defaultValue, annotation))
    return FALSE;
  readConfig(config, "inputfileStarCamera", var.fileNameStarCamera, Config::MUSTSET,  "starcamera.{satellite}.dat", "instrument file containing star camera data");
  readConfig(config, "identifier",          var.identifier,         Config::MUSTSET,  "", "string identifier (e.g. GNSS PRN: G01)");
  endSequence(config);
  return TRUE;
}

/***********************************************/

void starCamera2Orbex::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {

    FileName fileNameOrbex, fileNameAttitude;
    std::vector<Satellite>   satellites;
    std::vector<FileName> fileNameTransmitterList;
    std::string variablePrn, coordSystem, contact, createdBy, description, inputData;
    std::vector<std::string> comments;
    TimeSeriesPtr timeSeriesPtr;
    EarthRotationPtr earthRotation;
    UInt interpolationDegree;
    Polynomial polynomial;

    readConfig(config, "outputfileOrbex",          fileNameOrbex,           Config::MUSTSET,  "",      "ORBEX file");
    readConfig(config, "satellite",                satellites,              Config::MUSTSET,  "",      "");
    readConfig(config, "timeSeries",               timeSeriesPtr,           Config::DEFAULT,  "",      "resample to these epochs (otherwise input file epochs are used)");
    readConfig(config, "earthRotation",            earthRotation,           Config::OPTIONAL, "file",  "rotate data into Earth-fixed frame");
    readConfig(config, "interpolationDegree",      interpolationDegree,     Config::MUSTSET,  "7",     "for attitude and Earth rotation interpolation");
    readConfig(config, "description",              description,             Config::MUSTSET,  "",      "description of file contents");
    readConfig(config, "createdBy",                createdBy,               Config::MUSTSET,  "",      "name of agency");
    readConfig(config, "inputData",                inputData,               Config::MUSTSET,  "p",     "description of input data (see ORBEX description)");
    readConfig(config, "contact",                  contact,                 Config::MUSTSET,  "",      "email address");
    readConfig(config, "referenceFrame",           coordSystem,             Config::MUSTSET,  "IGb14", "reference frame used in file");
    readConfig(config, "comment",                  comments,                Config::OPTIONAL, "",      "");
    if(isCreateSchema(config)) return;

    std::vector<Time> times = timeSeriesPtr->times();

    std::set<Time> timesSet; // unique set of times covering all satellites
    std::map<Time, std::vector<Record>> times2Records;

    auto iter = satellites.begin();
    while(iter != satellites.end())
    {
      try
      {
        logStatus<<"read star camera file <"<<iter->fileNameStarCamera<<">"<<Log::endl;
        iter->starCamera = InstrumentFile::read(iter->fileNameStarCamera);

        for(UInt idEpoch=0; idEpoch<iter->starCamera.size(); idEpoch++)
        {
          iter->starCamera.at(idEpoch).rotary = inverse(iter->starCamera.at(idEpoch).rotary);
          if(earthRotation)
            iter->starCamera.at(idEpoch).rotary *= inverse(earthRotation->rotaryMatrix(iter->starCamera.at(idEpoch).time));
        }
        Matrix quaternion = iter->starCamera.matrix().column(1, 4);

        // (optionally) resample attitude
        std::vector<Time> prnTimes = iter->starCamera.times();
        if(times.size())
        {
          Polynomial polynomial(prnTimes, interpolationDegree);
          quaternion = polynomial.interpolate(times, quaternion);
          prnTimes = times;
        }

        for(UInt idEpoch=0; idEpoch<quaternion.rows(); idEpoch++)
          times2Records[prnTimes.at(idEpoch)].push_back(Record(iter->identifier, (quaternion.row(idEpoch)/norm(quaternion.row(idEpoch))).trans())); // normalize quaternions

      }
      catch(std::exception &e)
      {
        logWarning << e.what() << " continue..." << Log::endl;
        iter = satellites.erase(iter);
        continue;
      }
      iter->times = iter->starCamera.times();
      std::copy(iter->times.begin(), iter->times.end(), std::inserter(timesSet, timesSet.end()));
      iter++;
    }

    if (timesSet.empty())
      throw(Exception("No star camera input data found."));

    times.clear();
    std::copy(timesSet.begin(), timesSet.end(), std::inserter(times, times.end()));

    logStatus<<"writing ORBEX file to <"<<fileNameOrbex<<">"<<Log::endl;
    OutFile outfile(fileNameOrbex);
    outfile << "%=ORBEX  0.09" << std::endl;
    outfile << "%%" << std::endl;
    outfile << "*------------------------------------------------------------------------------------------------------" << std::endl;
    outfile << "+FILE/DESCRIPTION" << std::endl;
    outfile << " DESCRIPTION         " << description << std::endl;
    outfile << " CREATED_BY          " << createdBy << std::endl;
    outfile << " CREATION_DATE       " << System::now()%"%y %m %d %H %M %S"s << std::endl;
    outfile << " INPUT_DATA          " << inputData << std::endl;
    outfile << " CONTACT             " << contact << std::endl;
    outfile << " TIME_SYSTEM         GPS" << std::endl;
    outfile << " START_TIME          " << times.front()%"%y %m %d %H %M %012.9S000"s << std::endl;// ATTENTION: "fake" 15.12 precision due to max Time precision
    outfile << " END_TIME            " << times.back() %"%y %m %d %H %M %012.9S000"s << std::endl;// ATTENTION: "fake" 15.12 precision due to max Time precision
    outfile << " EPOCH_INTERVAL      " << medianSampling(times).seconds()%"%9.3f"s << std::endl;
    outfile << " COORD_SYSTEM        " << coordSystem << std::endl;
    outfile << " FRAME_TYPE          " << (earthRotation ? "ECEF" : "ECI") << std::endl;
    outfile << " LIST_OF_REC_TYPES   ATT" << std::endl;
    outfile << "-FILE/DESCRIPTION" << std::endl;
    if(comments.size())
    {
      outfile << "*------------------------------------------------------------------------------------------------------" << std::endl;
      for(auto &comment : comments)
        outfile << "* " << comment << std::endl;
    }
    outfile << "*------------------------------------------------------------------------------------------------------" << std::endl;
    outfile << "+SATELLITE/ID_AND_DESCRIPTION" << std::endl;
    for(const auto &iter : satellites)
      outfile << " " << iter.identifier << std::endl;
    outfile << "-SATELLITE/ID_AND_DESCRIPTION" << std::endl;
    outfile << "*------------------------------------------------------------------------------------------------------" << std::endl;
    outfile << "+EPHEMERIS/DATA" << std::endl;
    outfile << "*ATT RECORDS: TRANSFORMATION FROM TERRESTRIAL FRAME COORDINATES (T) TO SAT. BODY FRAME ONES (B) SUCH AS" << std::endl;
    outfile << "*                                 (0,B) = q.(0,T).trans(q)" << std::endl;
    outfile << "*REC ID_              N ___q0_(scalar)_____ ____q1__x__________ ____q2__y__________ ____q3__z__________" << std::endl;
    for(const auto &epoch : times2Records)
    {
      outfile << epoch.first%"## %y %m %d %H %M %012.9S000 "s << epoch.second.size()%"%3i"s << std::endl; // ATTENTION: "fake" 15.12 precision due to max Time precision
      for(const auto &record : epoch.second)
      {
        outfile << " ATT " << std::left << std::setw(17) << record.identifier << std::right << record.quaternion.size();
        for(UInt i = 0; i < record.quaternion.size(); i++)
          outfile << record.quaternion(i)%" %19.16f"s;
        outfile << std::endl;
      }
    }
    outfile << "-EPHEMERIS/DATA" << std::endl;
    outfile << "*------------------------------------------------------------------------------------------------------" << std::endl;
    outfile << "%END_ORBEX" << std::endl;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
