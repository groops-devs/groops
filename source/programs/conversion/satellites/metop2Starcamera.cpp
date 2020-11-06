/***********************************************/
/**
* @file metop2Starcamera.cpp
*
* @brief read MetOp star camera data.
*
* @author Torsten Mayer-Guerr
* @date 2010-07-26
* @author Norbert Zehentner
* @date 2014-05-12
*
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program reads in star camera data from MetOp satellites given in the special CHAMP format.
A description of the format can be found under: \url{http://op.gfz-potsdam.de/champ/docs_CHAMP/CH-GFZ-FD-001.pdf}
)";

/***********************************************/

#include "programs/program.h"
#include "base/string.h"
#include "files/fileMatrix.h"
#include "files/fileInstrument.h"

/***** CLASS ***********************************/

/** @brief read MetOp star camera data.
* @ingroup programsConversionGroup */
class Metop2Starcamera
{
  void readFileMetop(const FileName &fileName, StarCameraArc &starArc);
  void fillStarCamera(StarCameraArc &starArc);

public:
  void run(Config &config);
};

GROOPS_REGISTER_PROGRAM(Metop2Starcamera, SINGLEPROCESS, "read MetOp star camera data", Conversion, Instrument)

/***********************************************/

void Metop2Starcamera::run(Config &config)
{
  try
  {
    FileName outNameStar;
    std::vector<FileName> fileName;

    readConfig(config, "outputfileStarCamera", outNameStar, Config::MUSTSET,  "", "");
    readConfig(config, "inputfile", fileName, Config::MUSTSET,  "", "");
    if(isCreateSchema(config)) return;

    logStatus<<"read input files"<<Log::endl;
    StarCameraArc starArc;
    for(UInt i=0; i<fileName.size(); i++)
    {
      logStatus<<"read file <"<<fileName.at(i)<<">"<<Log::endl;
      readFileMetop(fileName.at(i), starArc);
    }

//     fillStarCamera(starArc);

    // Daten speichern
    // ---------------
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

void Metop2Starcamera::readFileMetop(const FileName &fileName, StarCameraArc &starArc)
{
  try
  {
    StarCameraEpoch starEpoch;
    Bool     attFlag   = FALSE;

    std::ifstream file(fileName.c_str());
    if(!file.good())
    {
      logWarning<<"cannot open file: "<<fileName.str()<<", continue..."<<Log::endl;
      return;
    }
    file.exceptions(std::ios::badbit|std::ios::failbit);

    //Daten einlesen, Headerzeilen werden hier direkt behandelt
//    logTimerStart;
    for(UInt i=0; ; i++)
    {
//      logTimerLoop(i, 259200);
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
       // Flag setzen, da diese nicht in allen epochen vorhanden sind
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
        starEpoch.time = date2time(year, month, day, hour, minute, second);
      }
      else if(lineID1 == "att")  // Attitude einlesen (als Quaternionen gegeben) Reihenfolge: Vektorkomponenten 1,2,3 unnd sklare Komponente
      {
        Vector q(4);
        q(1) = String::toDouble(line.substr(8, 14));
        q(2) = String::toDouble(line.substr(22, 14));
        q(3) = String::toDouble(line.substr(36, 14));
        q(0) = String::toDouble(line.substr(50, 14));
        attFlag = TRUE;
        if(fabs(norm(q)-1)>1e-5)
        {
          logWarning<<"strange norm = "<<norm(q)<<Log::endl;
          continue;
        }
        starEpoch.rotary = Rotary3d(q);
      }
      else if(lineID1 == "sca") {}  // Quaternion to convert Spacecraft reference frame
      else if(lineID1 == "ang") {}  // Euler angles
      else if(lineID1 == "pve") {}  // Position and velocity in Earth Centered Fixed coordinates
      else if(lineID1 == "pvi") {}  // Position and velocity in Inertial True-of-Epoch coordinates
      else if(lineID1 == "%eo")
        starArc.push_back(starEpoch);
    }  //for(UInt i=0; ; i++) Schleife über die zeilen des Files
//    logTimerLoopEnd(259200);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}
/***********************************************/

void Metop2Starcamera::fillStarCamera(StarCameraArc &starArc)
{
  try
  {
    Matrix quat(starArc.size(),5);
    for(UInt idStar=0; idStar<starArc.size(); idStar++)
    {
      Vector tmp = starArc.at(idStar).rotary.quaternion();
      quat(idStar,0)=starArc.at(idStar).time.mjd();
      quat(idStar,1)=tmp(0);
      quat(idStar,2)=tmp(1);
      quat(idStar,3)=tmp(2);
      quat(idStar,4)=tmp(3);
    }
    writeFileMatrix(FileName("quaternions.txt"), quat);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/
