/***********************************************/
/**
* @file cosmic2OrbitStar.cpp
*
* @brief read cosmic orbit star camera data.
*
* @author Norbert Zehentner
* @date 2014-07-28
*
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program reads in cosmic orbit and star camera data given in the CHAMP format.
In case of cosmic orbit and star camera data is stored in one file.
A description of the format can be found under: \url{http://op.gfz-potsdam.de/champ/docs_CHAMP/CH-GFZ-FD-001.pdf}
)";

/***********************************************/

#include "programs/program.h"
#include "base/string.h"
#include "files/fileInstrument.h"

/***** CLASS ***********************************/

/** @brief read cosmic orbit and star camera data.
* @ingroup programsConversionGroup */
class Cosmic2OrbitStar
{
  void readFileCosmic(const FileName &fileName, StarCameraArc &starArc, OrbitArc &orbArc);

public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(Cosmic2OrbitStar, SINGLEPROCESS, "read COSMIC orbit and star camera data", Conversion, Orbit, Instrument)

/***********************************************/

void Cosmic2OrbitStar::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName outNameOrb;
    FileName outNameStar;
    std::vector<FileName> fileName;

    readConfig(config, "outputfileOrbit",         outNameOrb,  Config::MUSTSET,  "", "");
    readConfig(config, "outputfileStarCamera",    outNameStar, Config::MUSTSET,  "", "");
    readConfig(config, "inputfile",               fileName,    Config::MUSTSET,  "", "");
    if(isCreateSchema(config)) return;

    logStatus<<"read input files"<<Log::endl;
    OrbitArc orbArc;
    StarCameraArc starArc;
    for(UInt i=0; i<fileName.size(); i++)
    {
      logStatus<<"read file <"<<fileName.at(i)<<">"<<Log::endl;
      readFileCosmic(fileName.at(i), starArc, orbArc);
    }

    // Daten speichern
    // ---------------
    logStatus<<"write orbit to file <"<<outNameOrb<<">"<<Log::endl;
    std::list<Arc> arcList3; arcList3.push_back(orbArc);
    InstrumentFile::write(outNameOrb, arcList3);
    logStatus<<"write star camera data to file <"<<outNameStar<<">"<<Log::endl;
    std::list<Arc> arcList2; arcList2.push_back(starArc);
    InstrumentFile::write(outNameStar, arcList2);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void Cosmic2OrbitStar::readFileCosmic(const FileName &fileName, StarCameraArc &starArc, OrbitArc &orbArc)
{
  try
  {
    OrbitEpoch orbEpoch;
    StarCameraEpoch starEpoch;
    Bool     attFlag   = FALSE;
    Bool     orbFlag   = FALSE;

    std::ifstream file(fileName.c_str());
    if(!file.good())
    {
      logWarning<<"cannot open file: "<<fileName.str()<<", continue..."<<Log::endl;
      return;
    }
    file.exceptions(std::ios::badbit|std::ios::failbit);

    //Daten einlesen, Headerzeilen werden hier direkt behandelt
    for(UInt i=0; ; i++)
    {
      std::string line;
      try
      {
        getline(file, line);
      }
      catch(std::exception &/*e*/)
      {
        //logWarning<<std::endl<<e.what()<<" continue..."<<Log::endl;
        break;
      }

      if(line.empty())
        break;

      std::string lineID1 = line.substr(0,3);

      if(lineID1 == "tim")
      {
        // vorherige Epoche Korrekturen anbringen falls vorhanden und anschließen zurückgeben

        if(attFlag)
        {
          if(starArc.size() && (starEpoch.time <= starArc.at(starArc.size()-1).time))
          {
            logWarning<<"(starEpoch.time <= starArc.at(starArc.size()-1).time)"<<Log::endl;
            continue;
          }
          starArc.push_back(starEpoch);
        }
        if(orbFlag)
        {
          if(orbArc.size() && (orbEpoch.time <= orbArc.at(orbArc.size()-1).time))
          {
            logWarning<<"(orbEpoch.time <= orbArc.at(orbArc.size()-1).time)"<<Log::endl;
            continue;
          }
          orbArc.push_back(orbEpoch);
        }
       // Flags für Korrekturen setzen, da diese nicht in allen epochen vorhanden sind
        orbFlag   = FALSE;
        attFlag   = FALSE;

        // epoch time der neuen Epoche einlesen
        UInt year, month, day, hour, minute;
        Double second;
        year = String::toInt(line.substr(4, 4));
        month = String::toInt(line.substr(9, 2));
        day = String::toInt(line.substr(12, 2));
        hour = String::toInt(line.substr(15, 2));
        minute = String::toInt(line.substr(18, 2));
        second = String::toDouble(line.substr(21, 10));
        orbEpoch.time = date2time(year, month, day, hour, minute, second);
        starEpoch.time = date2time(year, month, day, hour, minute, second);
      }
      else if(lineID1 == "pvi")  // orbit position and velocity
      {
        if(line.find("nan")==std::string::npos)
        {
          orbEpoch.position.x() = String::toDouble(line.substr(3, 14)) *1000;   // km => m
          orbEpoch.position.y() = String::toDouble(line.substr(17, 14)) *1000;   // km => m
          orbEpoch.position.z() = String::toDouble(line.substr(31, 14)) *1000;   // km => m
          orbEpoch.velocity.x() = String::toDouble(line.substr(45, 14)) *1000;   // km => m
          orbEpoch.velocity.y() = String::toDouble(line.substr(59, 14)) *1000;   // km => m
          orbEpoch.velocity.z() = String::toDouble(line.substr(73, 14)) *1000;   // km => m
          if(orbEpoch.position.r()<DEFAULT_R)
          {
  //           logWarning<<"strange position = "<<orbEpoch.position.r()<<Log::endl;
            continue;
          }
          orbFlag = TRUE;
        }
      }
      else if(lineID1 == "att")  // Attitude einlesen (als Quaternionen gegeben) Reihenfolge: Vektorkomponenten 1,2,3 unnd sklare Komponente
      {
        if(line.find("nan")==std::string::npos)
        {
          Vector q(4);
          q(1) = String::toDouble(line.substr(8, 14));
          q(2) = String::toDouble(line.substr(22, 14));
          q(3) = String::toDouble(line.substr(36, 14));
          q(0) = String::toDouble(line.substr(50, 14));
          if(fabs(norm(q)-1)>1e-5)
          {
  //           logWarning<<"strange norm = "<<norm(q)<<Log::endl;
            continue;
          }
          attFlag = TRUE;
          starEpoch.rotary = Rotary3d(q);
        }
      }
    }  //for(UInt i=0; ; i++) Schleife über die zeilen des Files
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}
/***********************************************/

