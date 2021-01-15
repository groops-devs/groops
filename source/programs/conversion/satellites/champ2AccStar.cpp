/***********************************************/
/**
* @file champ2AccStar.cpp
*
* @brief read CHAMP accelerometer and star camera data.
*
* @author Torsten Mayer-Guerr
* @author Norbert Zehentner
* @date 2010-07-26
*
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program reads in CHAMP accelerometer and star camera data given in the special CHAMP format.
In case of CHAMP accelerometer and star camera data is both stored in one file.
A description of the format can be found under: \url{http://op.gfz-potsdam.de/champ/docs_CHAMP/CH-GFZ-FD-001.pdf}.
)";

/***********************************************/

#include "programs/program.h"
#include "base/string.h"
#include "inputOutput/file.h"
#include "files/fileInstrument.h"

/***** CLASS ***********************************/

/** @brief read CHAMP accelerometer and star camera data.
* @ingroup programsConversionGroup */
class Champ2AccStar
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(Champ2AccStar, SINGLEPROCESS, "read CHAMP accelerometer and star camera data", Conversion, Instrument)

/***********************************************/

void Champ2AccStar::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName              outNameAcc, outNameAca, outNameStar;
    std::vector<FileName> fileNameIn;

    readConfig(config, "outputfileAccelerometer",       outNameAcc,  Config::OPTIONAL, "", "");
    readConfig(config, "outputfileAngularAcceleration", outNameAca,  Config::OPTIONAL, "", "");
    readConfig(config, "outputfileStarCamera",          outNameStar, Config::OPTIONAL, "", "");
    readConfig(config, "inputfile",                     fileNameIn,  Config::MUSTSET,  "", "");
    if(isCreateSchema(config)) return;

    logStatus<<"read input files"<<Log::endl;
    AccelerometerArc arcAcc;
    AccelerometerArc arcAca;
    StarCameraArc    arcSca;
    for(UInt idFile=0; idFile<fileNameIn.size(); idFile++)
    {
      try
      {
        logStatus<<"read file <"<<fileNameIn.at(idFile)<<">"<<Log::endl;
        InFile file(fileNameIn.at(idFile));

        AccelerometerEpoch epochAcc;
        AccelerometerEpoch epochAca;
        StarCameraEpoch    epochSca;
        Vector3d           acl_k0, acl_k1, acc1, acc2;
        Bool biasFlag  = FALSE;
        Bool scaleFlag = FALSE;
        Bool aclFlag   = FALSE;
        Bool acaFlag   = FALSE;
        Bool attFlag   = FALSE;
        Bool acc1Flag  = FALSE;
        Bool acc2Flag  = FALSE;

        std::string line;
        while(std::getline(file, line))
        {
          if(line.empty())
            continue;
          std::stringstream ss(line);
          ss.exceptions(std::ios::badbit | std::ios::failbit);

          // Header Information (corrections for the accelerometer)
          if(line.find("+acl_k0") == 0)
          {
            acl_k0.x() = String::toDouble(line.substr(12, 14));
            acl_k0.y() = String::toDouble(line.substr(28, 14));
            acl_k0.z() = String::toDouble(line.substr(44, 14));
            biasFlag = TRUE;
          }
          else if(line.find("+acl_k1") == 0)
          {
            acl_k1.x() = String::toDouble(line.substr(12, 14));
            acl_k1.y() = String::toDouble(line.substr(28, 14));
            acl_k1.z() = String::toDouble(line.substr(44, 14));
            scaleFlag = TRUE;
          }

          // new epoch? => store old epoch
          if((line.find("tim") == 0) || (line.find("%eof") == 0))
          {
            if(aclFlag)
            {
              if(acc1Flag) epochAcc.acceleration += acc1;   // Lorentz correction
              if(acc2Flag) epochAcc.acceleration += acc2;   // correction for the radial component
              if(biasFlag) epochAcc.acceleration -= acl_k0; // bias correction
              if(scaleFlag)                                 // scale factor
              {
                epochAcc.acceleration.x() *= acl_k1.x();
                epochAcc.acceleration.y() *= acl_k1.y();
                epochAcc.acceleration.z() *= acl_k1.z();
              }
              // Umrechnung von mm/s^2  zu m/s^2 und Übergang von Akzelerometer System auf Satelliten Body system
              Vector3d tmp(epochAcc.acceleration.x()*0.001, epochAcc.acceleration.y()*0.001, epochAcc.acceleration.z()*0.001);
              epochAcc.acceleration.x() =  tmp.y();
              epochAcc.acceleration.y() = -tmp.z();
              epochAcc.acceleration.z() = -tmp.x();
              arcAcc.push_back(epochAcc);
            }

            if(attFlag)
            {
              arcSca.push_back(epochSca);
            }
            if(acaFlag)
            {
              Vector3d tmp(epochAca.acceleration.x(), epochAca.acceleration.y(), epochAca.acceleration.z());
              epochAca.acceleration.x() =  tmp.y();
              epochAca.acceleration.y() = -tmp.z();
              epochAca.acceleration.z() = -tmp.x();
              arcAca.push_back(epochAca);
            }
          }

          if(line.find("tim") == 0) // epoch time of the new epoch
          {
            const UInt   year   = String::toInt(line.substr(4, 4));
            const UInt   month  = String::toInt(line.substr(9, 2));
            const UInt   day    = String::toInt(line.substr(12, 2));
            const UInt   hour   = String::toInt(line.substr(15, 2));
            const UInt   minute = String::toInt(line.substr(18, 2));
            const Double second = String::toDouble(line.substr(21, 10));
            epochAcc.time = epochAca.time = epochSca.time = date2time(year, month, day, hour, minute, second);
            // start new epoch without data yet
            acc1Flag = acc2Flag = aclFlag = acaFlag = attFlag = FALSE;
          }
          else if(line.find("acl") == 0)   // accelerometer
          {
            epochAcc.acceleration.x() = String::toDouble(line.substr(8, 14));    // mm/s^2
            epochAcc.acceleration.y() = String::toDouble(line.substr(22, 14));   // mm/s^2
            epochAcc.acceleration.z() = String::toDouble(line.substr(36, 14));   // mm/s^2
            aclFlag = TRUE;
          }
          else if((line.find("aca") == 0) && (line.find("*******") == std::string::npos)) // angular accelerations
          {
            epochAca.acceleration.x() = String::toDouble(line.substr(8, 14))/1000;  // mrad/s^2 => rad/s^2
            epochAca.acceleration.y() = String::toDouble(line.substr(22, 14))/1000;  // mrad/s^2 => rad/s^2
            epochAca.acceleration.z() = String::toDouble(line.substr(36, 14))/1000;  // mrad/s^2 => rad/s^2
            acaFlag = TRUE;
          }
          else if(line.find("att") == 0)   // attitude
          {
            Vector q(4);
            q(1) = String::toDouble(line.substr(8, 14));
            q(2) = String::toDouble(line.substr(22, 14));
            q(3) = String::toDouble(line.substr(36, 14));
            q(0) = String::toDouble(line.substr(50, 14));
            attFlag = TRUE;
            if(fabs(norm(q)-1)>1e-5)
            {
              logWarning<<"quaternion strange norm = "<<norm(q)<<Log::endl;
              continue;
            }
            epochSca.rotary = Rotary3d(q);
          }
          else if(line.find("thr") == 0) {}  // Thruster events
          else if(line.find("hka") == 0) {}  // housekeeping data
          else if(line.find("acc") == 0)     // accelerometer corrections
          {
            UInt tmp1 = String::toInt(line.substr(4, 2));
            UInt tmp2 = String::toInt(line.substr(7, 1));
            if(tmp1 == 1 && tmp2 == 0)
            {
              acc1.x() = String::toDouble(line.substr(8, 14));
              acc1.y() = String::toDouble(line.substr(22, 14));
              acc1.z() = String::toDouble(line.substr(36, 14));
              acc1Flag = TRUE;
            }
            if(tmp1 == 2 && tmp2 == 0)
            {
              acc2.x() = String::toDouble(line.substr(8, 14));
              acc2.y() = String::toDouble(line.substr(22, 14));
              acc2.z() = String::toDouble(line.substr(36, 14));
              acc2Flag = TRUE;
            }
          }
        } // while(getline)
      }
      catch(std::exception &e)
      {
        logError<<e.what()<<": continue..."<<Log::endl;
      }
    } // for(idFile)

    if(!outNameAcc.empty())
    {
      logStatus<<"write accelerometer data to file <"<<outNameAcc<<">"<<Log::endl;
      arcAcc.sort();
      Arc::printStatistics(arcAcc);
      InstrumentFile::write(outNameAcc, arcAcc);
    }

    if(!outNameStar.empty())
    {
      logStatus<<"write star camera data to file <"<<outNameStar<<">"<<Log::endl;
      arcSca.sort();
      Arc::printStatistics(arcSca);
      InstrumentFile::write(outNameStar, arcSca);
    }

    if(!outNameAca.empty())
    {
      logStatus<<"write angular accelerations to file <"<<outNameAca<<">"<<Log::endl;
      arcAca.sort();
      Arc::printStatistics(arcAca);
      InstrumentFile::write(outNameAca, arcAca);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
