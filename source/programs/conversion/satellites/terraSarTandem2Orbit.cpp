/***********************************************/
/**
* @file terraSarTandem2Orbit.cpp
*
* @brief read TerraSar-X or Tandem-X orbits from the special CHORB format
*
* @author Norbert Zehentner
* @date 2014-03-26
*
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program reads in TerraSar-X or Tandem-X orbits in the special CHORB format and takes the appropriate
time frame as stated in the document header.
A description of the format can be found under: \url{http://op.gfz-potsdam.de/champ/docs_CHAMP/CH-GFZ-FD-002.pdf}
)";

/***********************************************/

#include "programs/program.h"
#include "base/string.h"
#include "files/fileInstrument.h"
#include "classes/earthRotation/earthRotation.h"


/***** CLASS ***********************************/

/** @brief read TerraSar-X or Tandem-X orbits from the special CHORB format
* @ingroup programsConversionGroup */
class TerraSarTandem2Orbit
{
  EarthRotationPtr earthRotation;

  void readOrbit(const FileName &fileName, std::vector<Time> &times, std::vector<Vector3d> &position, std::vector<Vector3d> &velocity, std::vector<Time> &timeSpan) const;

public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(TerraSarTandem2Orbit, SINGLEPROCESS, "read TerraSar-X or Tandem-X orbits from the special CHORB format", Conversion, Orbit, Instrument)
GROOPS_RENAMED_PROGRAM(Tsxtdx2Orbit, Tsxtdx2Starcamera, date2time(2020, 8, 4))

/***********************************************/

void TerraSarTandem2Orbit::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName              orbitName;
    std::vector<FileName> inOrbitName;

    readConfig(config, "outputfileOrbit", orbitName,     Config::MUSTSET, "", "");
    readConfig(config, "earthRotation",   earthRotation, Config::MUSTSET, "", "");
    readConfig(config, "inputfile",       inOrbitName,   Config::MUSTSET, "", "orbits in CHORB format");
    if(isCreateSchema(config)) return;

    // ==============================

    OrbitArc  orbitArc;

    for(UInt k=0; k<inOrbitName.size(); k++)
    {
      std::vector<Time>     times;
      std::vector<Vector3d> position;
      std::vector<Vector3d> velocity;
      std::vector<Time>     timeSpan;
      UInt count=0;

      logStatus<<"read file <"<<inOrbitName.at(k)<<">"<<Log::endl;
      readOrbit(inOrbitName.at(k), times, position, velocity, timeSpan);
      logInfo<<"  pos:  "<<times.size()<<Log::endl;

      for(UInt i=0; i<times.size(); i++)
      {
        if(!timeSpan.empty())
          if((times.at(i)<timeSpan.at(0)) || (times.at(i)>timeSpan.at(1)))
            continue;
        Double r = position.at(i).r();
        if((r!=r)||(r<6380e3)||(r>100000e3))
        {
          logWarning<<"strange position r="<<r<<Log::endl;
          continue;
        }

        OrbitEpoch orbitEpoch;
        orbitEpoch.time     = times.at(i);
        orbitEpoch.position = position.at(i);
        orbitEpoch.velocity = velocity.at(i);
        orbitArc.push_back(orbitEpoch);
        count++;
      } // for(epoch i)
      logInfo<<"  used: "<<count<<Log::endl;
    }

    // Daten speichern
    // ---------------
    orbitArc.sort();
    orbitArc.removeDuplicateEpochs(TRUE);

    if(!orbitName.empty())
    {
      logStatus<<"write orbit data to file <"<<orbitName<<">"<<Log::endl;
      InstrumentFile::write(orbitName, orbitArc);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void TerraSarTandem2Orbit::readOrbit(const FileName &fileName, std::vector<Time> &times, std::vector<Vector3d> &position, std::vector<Vector3d> &velocity, std::vector<Time> &timeSpan) const
{
  try
  {
    times.clear();
    position.clear();
    velocity.clear();
    timeSpan.clear();
    std::string rFrame;

    // open file
    // -------------
    InFile file(fileName);
    file.exceptions(std::ios::badbit|std::ios::failbit);

    // read data
    // ---------
    Bool header = FALSE;
    Bool isEOF = FALSE;
    while(!isEOF)
    {
      if(!file.good())
        break;
      std::string line;
      try
      {
        getline(file, line);
      }
      catch(std::exception &/*e*/)
      {
//        logWarning<<std::endl<<e.what()<<" continue..."<<Log::endl;
        isEOF = TRUE;
        break;
      }
      if(line.size()<5)
      {
        isEOF = TRUE;
        break;
      }

      if(!header)
      {
        std::string lineID = line.substr(0,5);
        if(lineID=="RFRAM")
          rFrame = line.substr(7,3);
        else if(lineID=="COMEN" && line.substr(7,8)=="For best")
        {
          getline(file, line);
          timeSpan.push_back(seconds2time( String::toDouble(line.substr( 7, 2))*3600 + String::toDouble(line.substr(10, 2))*60 +  String::toDouble(line.substr(13, 2))));
          timeSpan.push_back(seconds2time( String::toDouble(line.substr(18, 2))*3600 + String::toDouble(line.substr(21, 2))*60 +  String::toDouble(line.substr(24, 2))));
        }
        else if(lineID=="ORBIT")
          header = TRUE;
      }
      else
      {
      // Epochenzeit
      // -----------
        UInt   day;
        Double sec, x, y, z, vx, vy, vz;
        day = String::toInt(line.substr(0, 6));
        sec = String::toDouble(line.substr(6, 11));
        x = String::toDouble(line.substr(17, 12));
        y = String::toDouble(line.substr(29, 12));
        z = String::toDouble(line.substr(41, 12));
        vx = String::toDouble(line.substr(53, 12));
        vy = String::toDouble(line.substr(65, 12));
        vz = String::toDouble(line.substr(77, 12));
       //Epochezeitpunkt ist in TT seit J2000 gegeben -> Umrechnung in GPS Zeit
        Time time = mjd2time(J2000);
        time += seconds2time((day)*pow(10,-1)*86400);
        time += seconds2time(sec*pow(10,-6));
        time = timeTT2GPS(time);
//         if((!times.empty()) && (time<=times.back()))
//         {
//           logWarning<<"new  time: "<<time.dateTimeStr()<<Log::endl;
//           logWarning<<"last time: "<<times.back().dateTimeStr()<<Log::endl;
//           throw(Exception("times not in increasing order!"));
//         }
        if(times.empty())
        {
          if(!timeSpan.empty())
            timeSpan.at(0) += seconds2time(static_cast<Double>(time.mjdInt())*86400);
          logInfo<<"  First epoch: "<<time.dateTimeStr()<<Log::endl;
        }

        times.push_back(time);
        if(rFrame=="CTS")
        {
          Rotary3d rotation  = inverse(earthRotation->rotaryMatrix(time));
          position.push_back(rotation.rotate(0.001*Vector3d(x,y,z))); // mm -> m
          velocity.push_back(rotation.rotate(0.0000001*Vector3d(vx,vy,vz)));  // 10^-7 m/s -> 1 m/s
        }
        else
        {
        position.push_back(0.001*Vector3d(x,y,z)); // mm -> m
        velocity.push_back(0.0000001*Vector3d(vx,vy,vz));  // 10^-7 m/s -> 1 m/s
        }
      }
    } // while(!EOF)
    if(!timeSpan.empty() && times.size()>15)
    {
      timeSpan.at(1) += seconds2time(static_cast<Double>(times.at(times.size()-10).mjdInt())*86400);
      logInfo<<"  Cut orbit to: "<<timeSpan.at(0).dateTimeStr()<<" < "<<timeSpan.at(1).dateTimeStr()<<Log::endl;
    }

    if(!isEOF)
      logWarning<<"do not reached the end of file"<<Log::endl;
  }
  catch(std::exception &/*e*/)
  {
    logWarning<<"cannot open file: "<<fileName<<", continue..."<<Log::endl;
    return;
//     GROOPS_RETHROW(e)
  }
}

/***********************************************/
