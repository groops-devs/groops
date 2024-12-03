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
Converts satellite attitude from \file{instrument file (STARCAMERA)}{instrument} to
\href{http://acc.igs.org/misc/ORBEX009.pdf}{ORBEX file format} (quaternions).

If \configClass{earthRotation}{earthRotationType} is provided, the output file contains quaternions
for rotation from TRF to satellite body frame (IGS/ORBEX convention),
otherwise the rotation is from CRF to satellite body frame.

See also \program{GnssOrbex2StarCamera}, \program{SimulateStarCameraGnss}.
)";

/***********************************************/

#include "programs/program.h"
#include "base/polynomial.h"
#include "classes/earthRotation/earthRotation.h"
#include "classes/timeSeries/timeSeries.h"
#include "files/fileInstrument.h"
#include "inputOutput/system.h"

/***** CLASS ***********************************/

/** @brief Converts satellite attitude from StarCamera to ORBEX format.
* @ingroup programsConversionGroup */
class StarCamera2Orbex
{
public:
  class Satellite
  {
  public:
    FileName    fileNameAttitude;
    std::string identifier, description;
  };

  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(StarCamera2Orbex, SINGLEPROCESS, "Converts satellite attitude from StarCamera to ORBEX format.", Conversion, Instrument)

/***********************************************/

template<> Bool readConfig(Config &config, const std::string &name, StarCamera2Orbex::Satellite &var, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation)
{
  if(!readConfigSequence(config, name, mustSet, defaultValue, annotation))
    return FALSE;
  readConfig(config, "inputfileStarCamera", var.fileNameAttitude, Config::MUSTSET,  "", "");
  readConfig(config, "identifier",          var.identifier,       Config::MUSTSET,  "", "string identifier (e.g. GNSS PRN: G01)");
  readConfig(config, "description",         var.description,      Config::OPTIONAL, "", "e.g. BLOCK IIR-B, GRACE");
  endSequence(config);
  return TRUE;
}

/***********************************************/

void StarCamera2Orbex::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName                 fileNameOrbex;
    std::vector<Satellite>   satellites;
    EarthRotationPtr         earthRotation;
    TimeSeriesPtr            timeSeriesPtr;
    UInt                     interpolationDegree;
    std::string              coordSystem, contact, createdBy, description, inputData;
    std::vector<std::string> comments;

    readConfig(config, "outputfileOrbex",     fileNameOrbex,       Config::MUSTSET,  "",      "ORBEX file");
    readConfig(config, "satellite",           satellites,          Config::MUSTSET,  "",      "");
    readConfig(config, "earthRotation",       earthRotation,       Config::OPTIONAL, "file",  "rotate data into Earth-fixed frame");
    readConfig(config, "timeSeries",          timeSeriesPtr,       Config::DEFAULT,  "",      "resample to these epochs (otherwise input file epochs are used)");
    readConfig(config, "interpolationDegree", interpolationDegree, Config::MUSTSET,  "7",     "for attitude and Earth rotation interpolation");
    readConfig(config, "description",         description,         Config::MUSTSET,  "",      "description of file contents");
    readConfig(config, "createdBy",           createdBy,           Config::MUSTSET,  "",      "name of agency");
    readConfig(config, "inputData",           inputData,           Config::MUSTSET,  "p",     "description of input data (see ORBEX description)");
    readConfig(config, "contact",             contact,             Config::MUSTSET,  "",      "email address");
    readConfig(config, "referenceFrame",      coordSystem,         Config::MUSTSET,  "IGS20", "reference frame used in file");
    readConfig(config, "comment",             comments,            Config::OPTIONAL, "",      "");
    if(isCreateSchema(config)) return;

    class Record
    {
    public:
      std::string identifier;
      Vector      quaternion;

      Record(const std::string &identifier, const_MatrixSliceRef quaternion) : identifier(identifier), quaternion(quaternion) {}
    };

    std::map<Time, std::vector<Record>> times2Records;
    std::vector<Time> times = timeSeriesPtr->times();
    for(auto &satellite : satellites)
    {
      satellite.identifier.resize(6, ' '); // max. 6 characters

      try
      {
        StarCameraArc arc = InstrumentFile::read(satellite.fileNameAttitude);
        for(UInt idEpoch=0; idEpoch<arc.size(); idEpoch++)
          arc.at(idEpoch).rotary = inverse(arc.at(idEpoch).rotary);
        if(earthRotation)
          for(UInt idEpoch=0; idEpoch<arc.size(); idEpoch++)
          arc.at(idEpoch).rotary *= inverse(earthRotation->rotaryMatrix(arc.at(idEpoch).time));
        Matrix quaternion = arc.matrix().column(1, 4);

        // (optionally) resample attitude
        std::vector<Time> satTimes = arc.times();
        if(times.size())
        {
          Polynomial polynomial(satTimes, interpolationDegree);
          quaternion = polynomial.interpolate(times, quaternion);
          satTimes = times;
        }

        for(UInt idEpoch=0; idEpoch<quaternion.rows(); idEpoch++)
          times2Records[satTimes.at(idEpoch)].push_back(Record(satellite.identifier, (quaternion.row(idEpoch)/norm(quaternion.row(idEpoch))).trans())); // normalize quaternions
      }
      catch(std::exception &e)
      {
        logWarning<<e.what()<<" continue..."<< Log::endl;
        satellite.identifier = "";
        continue;
      }
    }

    if(times2Records.empty())
      throw(Exception("No star camera input data found."));

    times.clear();
    for(const auto &epoch : times2Records)
      times.push_back(epoch.first);

    logStatus<<"writing ORBEX file to <"<<fileNameOrbex<<">"<<Log::endl;
    OutFile outfile(fileNameOrbex);
    outfile<<"%=ORBEX  0.09"<<std::endl;
    outfile<<"%%"<<std::endl;
    outfile<<"*------------------------------------------------------------------------------------------------------"<<std::endl;
    outfile<<"+FILE/DESCRIPTION"<<std::endl;
    outfile<<" DESCRIPTION         "<<description<<std::endl;
    outfile<<" CREATED_BY          "<<createdBy<<std::endl;
    outfile<<" CREATION_DATE       "<<System::now()%"%y %m %d %H %M %S"s<<std::endl;
    outfile<<" INPUT_DATA          "<<inputData<<std::endl;
    outfile<<" CONTACT             "<<contact<<std::endl;
    outfile<<" TIME_SYSTEM         GPS"<<std::endl;
    outfile<<" START_TIME          "<<times.front()%"%y %m %d %H %M %012.9S000"s<<std::endl; // ATTENTION: "fake" 15.12 precision due to max Time precision
    outfile<<" END_TIME            "<<times.back() %"%y %m %d %H %M %012.9S000"s<<std::endl; // ATTENTION: "fake" 15.12 precision due to max Time precision
    outfile<<" EPOCH_INTERVAL      "<<medianSampling(times).seconds()%"%9.3f"s<<std::endl;
    outfile<<" COORD_SYSTEM        "<<coordSystem<<std::endl;
    outfile<<" FRAME_TYPE          "<<(earthRotation ? "ECEF" : "ECI")<<std::endl;
    outfile<<" LIST_OF_REC_TYPES   ATT"<<std::endl;
    outfile<<"-FILE/DESCRIPTION"<<std::endl;
    if(comments.size())
    {
      outfile<<"*------------------------------------------------------------------------------------------------------"<<std::endl;
      for(auto &comment : comments)
        outfile<<"* "<<comment<<std::endl;
    }
    outfile<<"*------------------------------------------------------------------------------------------------------"<<std::endl;
    outfile<<"+SATELLITE/ID_AND_DESCRIPTION"<<std::endl;
    for(const auto &satellite : satellites)
      if(!satellite.identifier.empty())
        outfile<<" "<<satellite.identifier<<satellite.description<<std::endl;
    outfile<<"-SATELLITE/ID_AND_DESCRIPTION"<<std::endl;
    outfile<<"*------------------------------------------------------------------------------------------------------"<<std::endl;
    outfile<<"+EPHEMERIS/DATA"<<std::endl;
    outfile<<"*ATT RECORDS: TRANSFORMATION FROM "<<(earthRotation ? "TERRESTRIAL" : "CELESTIAL")<<" FRAME COORDINATES (T) TO SAT. BODY FRAME ONES (B) SUCH AS"<<std::endl;
    outfile<<"*                                 (0,B) = q.(0,T).trans(q)"<<std::endl;
    outfile<<"*REC ID_              N ___q0_(scalar)_____ ____q1__x__________ ____q2__y__________ ____q3__z__________"<<std::endl;
    for(const auto &epoch : times2Records)
    {
      outfile<<epoch.first%"## %y %m %d %H %M %012.9S000 "s<<epoch.second.size()%"%3i"s<<std::endl; // ATTENTION: "fake" 15.12 precision due to max Time precision
      for(const auto &record : epoch.second)
      {
        outfile<<" ATT "<<record.identifier<<std::string(11, ' ')<<record.quaternion.size();
        for(UInt i=0; i<record.quaternion.size(); i++)
          outfile<<record.quaternion(i)%" %19.16f"s;
        outfile<<std::endl;
      }
    }
    outfile<<"-EPHEMERIS/DATA"<<std::endl;
    outfile<<"*------------------------------------------------------------------------------------------------------"<<std::endl;
    outfile<<"%END_ORBEX"<<std::endl;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
