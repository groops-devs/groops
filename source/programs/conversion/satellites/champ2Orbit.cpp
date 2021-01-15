/***********************************************/
/**
* @file champ2Orbit.cpp
*
* @brief read champ PSO orbits from the special CHORB format
*
* @author Norbert Zehentner
* @date 2011-10-04
*
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program reads in CHAMP precise science orbits in the special CHORB format.
A description of the format can be found under: \url{http://op.gfz-potsdam.de/champ/docs_CHAMP/CH-GFZ-FD-002.pdf}
)";

/***********************************************/

#include "programs/program.h"
#include "base/string.h"
#include "files/fileInstrument.h"
#include "classes/earthRotation/earthRotation.h"
#include "classes/timeSeries/timeSeries.h"

/***** CLASS ***********************************/

/** @brief read champ PSO orbits from the special CHORB format
* @ingroup programsConversionGroup */
class Champ2Orbit
{
  std::string satID;

  void readOrbit     (const FileName &fileName, std::vector<Time> &times, std::vector<Vector3d> &position, std::vector<Vector3d> &velocity) const;

public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(Champ2Orbit, SINGLEPROCESS, "read champ PSO orbits from the special CHORB format", Conversion, Orbit, Instrument)

/***********************************************/

void Champ2Orbit::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName              orbitName;
    std::vector<FileName> inOrbitName;
    EarthRotationPtr      earthRotation;
    TimeSeriesPtr         timeSeries;

    readConfig(config, "outputfileOrbit",           orbitName,     Config::MUSTSET,  "", "");
    readConfig(config, "earthRotation",             earthRotation, Config::MUSTSET,  "", "");
    readConfig(config, "timeSeries",                timeSeries,    Config::DEFAULT,  "", "");
    if(readConfigSequence(config, "inputOrbit", Config::MUSTSET, "", ""))
    {
      readConfig(config, "inputfile", inOrbitName, Config::MUSTSET, "", "orbits in SP3 format");
      endSequence(config);
    }
    if(isCreateSchema(config)) return;

    // ==============================
    // File testing
    std::vector<Time> timeSpan;
    timeSpan = timeSeries->times();

    // compare timeSeries with number of data files, if no timeSeries is given section is skipped
    if(timeSpan.size()!=0)
    {
      if(timeSpan.size() != inOrbitName.size()+1)
      {
        logError<<"fileCount("<<inOrbitName.size()<<") +1 != timeCount("<<timeSpan.size()<<")"<<Log::endl;
        throw(Exception("fileCount != timeCount-1"));
      }
    }
    logError<<"size timeSpan: "<<timeSpan.size()<<Log::endl;


    // ==============================

    OrbitArc    orbitArc;
    std::string rFrame;

    for(UInt k=0; k<inOrbitName.size(); k++)
    {
      std::vector<Time>     times;
      std::vector<Vector3d> position;
      std::vector<Vector3d> velocity;

      logStatus<<"read file <"<<inOrbitName.at(k)<<">"<<Log::endl;
      readOrbit(inOrbitName.at(k), times, position, velocity);
      logInfo<<"  pos:  "<<times.size()<<Log::endl;

      std::string rFram = inOrbitName.at(k).baseName();
      rFrame = rFram.substr(12,3);
      Time tmpTimes;
      UInt count  = 0;
      for(UInt i=0; i<times.size(); i++)
      {
        if((!timeSpan.empty()) && ((times.at(i)<timeSpan.at(k)) || (times.at(i)>=timeSpan.at(k+1))))
          continue;
        Double r = position.at(i).r();
        if((r!=r)||(r<6380e3)||(r>100000e3))
        {
          logWarning<<"strange position r="<<r<<Log::endl;
          continue;
        }

        OrbitEpoch orbitEpoch;
        if((orbitArc.size()!=0) && (orbitArc.at(orbitArc.size()-1).time >= times.at(i)))
        {
          tmpTimes = times.at(i);
          logWarning<<"times not in increasing order!"<<Log::endl;
          logWarning<<"new  time: "<<times.at(i).dateTimeStr()<<Log::endl;
          logWarning<<"last time: "<<orbitArc.at(orbitArc.size()-1).time.dateTimeStr()<<Log::endl;
        }
        else
        {
          orbitEpoch.time     = times.at(i);
          orbitEpoch.position = position.at(i);
          orbitEpoch.velocity = velocity.at(i);
          orbitArc.push_back(orbitEpoch);
          count++;
        }

      } // for(epoch i)
      if(count==0)
      {
        logWarning<<"Files may be in wrong order: Actual file is not used!!!!!!!!!!!!"<<Log::endl;
        logWarning<<"last time in actual file: "<<tmpTimes.dateTimeStr()<<Log::endl;
        logWarning<<"first time in orbit     : "<<orbitArc.at(0).time.dateTimeStr()<<Log::endl;
        logWarning<<"last time in orbit      : "<<orbitArc.at(orbitArc.size()-1).time.dateTimeStr()<<Log::endl;
      }
      logInfo<<"  used: "<<count<<Log::endl;
    }

    // ==============================

    // Drehung ins Raumfeste System nur wenn im erdfesten gegeben (TDS (True of date system) oder CTS (conventional terrestrial system) steht im Header der datei)
    // ----------------------------
    if(rFrame=="CTS")
    {
      logStatus<<"rotation from TRF to CRF"<<Log::endl;
      Single::forEach(orbitArc.size(), [&](UInt i)
      {
        const Rotary3d rotation = inverse(earthRotation->rotaryMatrix(orbitArc.at(i).time));
        orbitArc.at(i).position = rotation.rotate(orbitArc.at(i).position);
        if(orbitArc.at(i).velocity.r() > 0)
          orbitArc.at(i).velocity = rotation.rotate(orbitArc.at(i).velocity) + crossProduct(earthRotation->rotaryAxis(orbitArc.at(i).time), orbitArc.at(i).position);
      });
    }
    // ==============================

    // Daten speichern
    // ---------------
    if(!orbitName.empty())
    {
      logStatus<<"write orbit data to file <"<<orbitName<<">"<<Log::endl;
      std::list<Arc> arcList; arcList.push_back(orbitArc);
      InstrumentFile::write(orbitName, arcList);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void Champ2Orbit::readOrbit(const FileName &fileName, std::vector<Time> &times, std::vector<Vector3d> &position, std::vector<Vector3d> &velocity) const
{
  try
  {
    times.clear();
    position.clear();
    velocity.clear();

    // open file
    // -------------
    std::ifstream file(fileName.c_str());
    if(!file.good())
    {
      logWarning<<"cannot open file: "<<fileName<<", continue..."<<Log::endl;
      return;
    }
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
        if(lineID=="RFRAM") {} // wird aus dem Dateinamen gelesen
        else if(lineID=="ORBIT") header = TRUE;
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
        if((!times.empty()) && (time<=times.back()))
        {
          logWarning<<"new  time: "<<time.dateTimeStr()<<Log::endl;
          logWarning<<"last time: "<<times.back().dateTimeStr()<<Log::endl;
          throw(Exception("times not in increasing order!"));
        }
        if(times.empty())
          logInfo<<"  First epoch: "<<time.dateTimeStr()<<Log::endl;

        times.push_back(time);
        position.resize(times.size());
        velocity.resize(times.size());
        position.back() = 0.001*Vector3d(x,y,z); // mm -> m
        velocity.back() = 0.0000001*Vector3d(vx,vy,vz);  // 10^-7 m/s -> 1 m/s
      }
    } // while(!EOF)

    if(!isEOF)
      logWarning<<"do not reached the end of file"<<Log::endl;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
