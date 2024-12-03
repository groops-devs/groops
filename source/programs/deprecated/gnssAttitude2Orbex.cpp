/***********************************************/
/**
* @file gnssAttitude2Orbex.cpp
*
* @brief DEPRECATED. Please use StarCamera2Orbex instead.
*
* @author Sebastian Strasser
* @date 2019-05-29
*
* @deprecated Please use StarCamera2Orbex instead.
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
DEPRECATED. Please use \program{StarCamera2Orbex} instead.
)";

/***********************************************/

#include "programs/program.h"
#include "base/polynomial.h"
#include "classes/earthRotation/earthRotation.h"
#include "classes/timeSeries/timeSeries.h"
#include "files/fileInstrument.h"
#include "files/fileStringTable.h"
#include "inputOutput/system.h"

/***** CLASS ***********************************/

/** @brief DEPRECATED. Please use StarCamera2Orbex instead.
* @ingroup programsConversionGroup */
class GnssAttitude2Orbex
{
  class Record
  {
  public:
    std::string prn;
    Vector      quaternion;

    Record(const std::string &prn, const_MatrixSliceRef quaternion) : prn(prn), quaternion(quaternion) {}
  };

public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(GnssAttitude2Orbex, SINGLEPROCESS, "DEPRECATED. Please use StarCamera2Orbex instead.", Deprecated)

/***********************************************/

void GnssAttitude2Orbex::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName fileNameOrbex, fileNameAttitude;
    std::vector<FileName> fileNameTransmitterList;
    std::string variablePrn, coordSystem, contact, createdBy, description, inputData;
    std::vector<std::string> comments;
    TimeSeriesPtr timeSeriesPtr;
    EarthRotationPtr earthRotation;
    UInt interpolationDegree;
    Polynomial polynomial;

    readConfig(config, "outputfileOrbex",          fileNameOrbex,           Config::MUSTSET,  "",      "ORBEX file");
    readConfig(config, "inputfileTransmitterList", fileNameTransmitterList, Config::MUSTSET,  "",      "ASCII list with transmitter PRNs");
    readConfig(config, "inputfileAttitude",        fileNameAttitude,        Config::MUSTSET,  "attitude.{prn}.dat", "instrument file containing attitude");
    readConfig(config, "variablePrn",              variablePrn,             Config::DEFAULT,  "prn",   "loop variable for PRNs from transmitter list");
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

    logWarning<<"DEPRECATED. Please use StarCamera2Orbex instead."<<Log::endl;

    std::vector<Time> times = timeSeriesPtr->times();

    std::vector<std::string> transmitterList;
    for(const auto &fileName : fileNameTransmitterList)
    {
      logStatus<<"reading transmitter list from <"<<fileName<<">"<<Log::endl;
      readFileStringList(fileName, transmitterList);
    }

    VariableList fileNameVariableList;
    std::map<Time, std::vector<Record>> times2Records;
    for(const auto &prn : transmitterList)
    {
      fileNameVariableList.setVariable(variablePrn, prn);
      StarCameraArc arc = InstrumentFile::read(fileNameAttitude(fileNameVariableList));

      for(UInt idEpoch=0; idEpoch<arc.size(); idEpoch++)
      {
        arc.at(idEpoch).rotary = inverse(arc.at(idEpoch).rotary);
        if(earthRotation)
          arc.at(idEpoch).rotary *= inverse(earthRotation->rotaryMatrix(arc.at(idEpoch).time));
      }
      Matrix quaternion = arc.matrix().column(1, 4);

      // (optionally) resample attitude
      std::vector<Time> prnTimes = arc.times();
      if(times.size())
      {
        Polynomial polynomial(prnTimes, interpolationDegree);
        quaternion = polynomial.interpolate(times, quaternion);
        prnTimes = times;
      }

      for(UInt idEpoch=0; idEpoch<quaternion.rows(); idEpoch++)
        times2Records[prnTimes.at(idEpoch)].push_back(Record(prn, (quaternion.row(idEpoch)/norm(quaternion.row(idEpoch))).trans())); // normalize quaternions
    }

    times.clear();
    for(const auto &epoch : times2Records)
      times.push_back(epoch.first);

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
    for(const auto &prn : transmitterList)
      outfile << " " << prn << std::endl;
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
        outfile << " ATT " << record.prn << std::string(14, ' ') << record.quaternion.size();
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
