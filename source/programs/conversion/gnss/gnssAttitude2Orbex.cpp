/***********************************************/
/**
* @file gnssAttitude2Orbex.cpp
*
* @brief Convert attitude of GNSS satellites to ORBEX file format (quaternions).
*
* @author Sebastian Strasser
* @date 2019-05-29
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
Convert attitude of GNSS satellites to \href{http://acc.igs.org/misc/proposal_orbex_april2019.pdf}{ORBEX file format} (quaternions).

If \configClass{earthRotation}{earthRotationType} is provided, the output file contains quaternions for rotation from TRF to satellite
body frame (IGS/ORBEX convention), otherwise the rotation is from CRF to satellite body frame.

See also \program{GnssOrbex2StarCamera}, \program{SimulateStarCameraGnss}.
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

/** @brief Convert attitude of GNSS satellites to ORBEX file format (quaternions).
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

GROOPS_REGISTER_PROGRAM(GnssAttitude2Orbex, SINGLEPROCESS, "Convert attitude of GNSS satellites to ORBEX file format (quaternions).", Conversion, Gnss)

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
    readConfig(config, "timeSeries",               timeSeriesPtr,           Config::OPTIONAL, "",      "resample to these epochs (otherwise input file epochs are used)");
    readConfig(config, "earthRotation",            earthRotation,           Config::OPTIONAL, "",      "rotate data into Earth-fixed frame");
    readConfig(config, "interpolationDegree",      interpolationDegree,     Config::MUSTSET,  "7",     "for attitude and Earth rotation interpolation");
    readConfig(config, "description",              description,             Config::MUSTSET,  "",      "description of file contents");
    readConfig(config, "createdBy",                createdBy,               Config::MUSTSET,  "",      "name of agency");
    readConfig(config, "inputData",                inputData,               Config::MUSTSET,  "p",     "description of input data (see ORBEX description)");
    readConfig(config, "contact",                  contact,                 Config::MUSTSET,  "",      "email address");
    readConfig(config, "referenceFrame",           coordSystem,             Config::MUSTSET,  "IGb14", "reference frame used in file");
    readConfig(config, "comment",                  comments,                Config::OPTIONAL, "",      "");
    if(isCreateSchema(config)) return;

    std::vector<Time> times;
    if(timeSeriesPtr)
    {
      times = timeSeriesPtr->times();
      polynomial = Polynomial(interpolationDegree);
    }

    std::vector<std::string> transmitterList;
    for(const auto &fileName : fileNameTransmitterList)
    {
      logStatus<<"reading transmitter list from <"<<fileName<<">"<<Log::endl;
      readFileStringList(fileName, transmitterList);
    }

    VariableList fileNameVariableList;
    addVariable(variablePrn, fileNameVariableList);
    std::map<Time, std::vector<Record>> times2Records;
    for(const auto &prn : transmitterList)
    {
      fileNameVariableList[variablePrn]->setValue(prn);

      StarCameraArc arc = InstrumentFile::read(fileNameAttitude(fileNameVariableList));
      Matrix quaternion(arc.size(), 4);

      for(UInt idEpoch = 0; idEpoch < arc.size(); idEpoch++)
      {
        Rotary3d rot = inverse(arc.at(idEpoch).rotary);
        if(earthRotation)
          rot *= inverse(earthRotation->rotaryMatrix(arc.at(idEpoch).time));
        copy(rot.quaternion().trans(), quaternion.row(idEpoch));
      }

      // ensure consistent sign of quaternions
      for(UInt idEpoch=1; idEpoch<quaternion.rows(); idEpoch++)
        if(inner(quaternion.row(idEpoch), quaternion.row(idEpoch-1))<0)
          quaternion.row(idEpoch) *= -1;

      // (optionally) resample attitude
      std::vector<Time> prnTimes = arc.times();
      if(times.size())
      {
        quaternion = polynomial.interpolate(times, prnTimes, quaternion);
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
