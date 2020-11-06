/***********************************************/
/**
* @file hw2TideGeneratingPotential.cpp
*
* @brief Read tide generating potential from Hartmann and Wenzel 1995.
* http://bowie.gsfc.nasa.gov/hw95/
*
* @author Torsten Mayer-Guerr
* @date 2011-09-17
*
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
Write \file{tide generating potential}{tideGeneratingPotential}
from Hartmann and Wenzel 1995 file, \url{https://doi.org/10.1029/95GL03324}.
)";

/***********************************************/

#include "programs/program.h"
#include "base/string.h"
#include "base/doodson.h"
#include "base/tideGeneratingPotential.h"
#include "inputOutput/file.h"
#include "files/fileTideGeneratingPotential.h"

/***** CLASS ***********************************/

/** @brief Read tide generating potential from Hartmann and Wenzel 1995.
* @ingroup programsConversionGroup */
class Hw2TideGeneratingPotential
{
public:
  void run(Config &config);
};

GROOPS_REGISTER_PROGRAM(Hw2TideGeneratingPotential, SINGLEPROCESS, "Read tide generating potential from Hartmann and Wenzel 1995", Conversion)

/***********************************************/

void Hw2TideGeneratingPotential::run(Config &config)
{
  try
  {
    FileName outName, inName;
    UInt     headerLines;
    Time     time0;

    readConfig(config, "outputfileTideGeneratingPotential", outName,     Config::MUSTSET, "",           "");
    readConfig(config, "inputfile",                         inName,      Config::MUSTSET, "hw95.dat",   "");
    readConfig(config, "headerLines",                       headerLines, Config::DEFAULT, "205",        "skip number of header lines");
    readConfig(config, "referenceTime",                     time0,       Config::DEFAULT, STRING_J2000, "reference time");
    if(isCreateSchema(config)) return;

    // =====================================================

    TideGeneratingPotential tgp;
    tgp.reserve(13000);

    UInt countPermanent = 0;
    UInt countNDegree   = 0;
    UInt countNDoodson  = 0;
    UInt countNSunMoon  = 0;
    Double maxAmp1      = 0.;
    Double maxAmp2      = 0.;
    Double maxAmp3      = 0.;

    Double T = timeGPS2JC(time0); // julian centuries

    // =====================================================

    InFile file(inName);

    // skip header
    std::string line;
    for(UInt i=0; i<headerLines; i++)
      getline(file, line);

    while(file.good())
    {
      std::string line;
      getline(file, line);
      if(line.empty())
        continue;

      // line number
      UInt number = String::toInt(line.substr(0, 6));
      if(number == 999999)
        break;

      std::string generatingBody = line.substr(7, 2);
      const UInt n = String::toInt(line.substr(9, 2));

      // doodson multipliers
      std::vector<Int> kn(11);
      for(UInt i=0; i<kn.size(); i++)
        kn.at(i) = String::toInt(line.substr(11+3*i, 3));

      // frequency
      const Double frequency = DEG2RAD/3600*String::toDouble(line.substr(44, 12)); // degree/hour -> rad/s

      // amplitudes
      Double c0 = 1e-10 * String::toDouble(line.substr(56, 12));
      Double s0 = 1e-10 * String::toDouble(line.substr(68, 12));
      c0 += 1e-10 * T * String::toDouble(line.substr(80, 10));
      s0 += 1e-10 * T * String::toDouble(line.substr(90, 10));

      // name
      const std::string name = String::trim(line.substr(101,4));

      // find permanent tides
      // --------------------
      {
        Int sum = 0;
        for(UInt i=0; i<kn.size(); i++)
          sum += static_cast<Int>(std::fabs(kn.at(i)));
        if(sum==0)
        {
          countPermanent++;
          continue;
        }
      }

      // use only degree n=2
      // -------------------
      if(n != 2)
      {
        countNDegree++;
        maxAmp1 = std::max(maxAmp1, sqrt(c0*c0+s0*s0));
        continue;
      }

      // use only tides from sun & moon
      // ------------------------------
      {
        Int sum = 0;
        for(UInt i=6; i<kn.size(); i++)
          sum += static_cast<Int>(std::fabs(kn.at(i)));
        if((sum!=0)||((generatingBody!="SU")&&(generatingBody!="MO")))
        {
          countNSunMoon++;
          maxAmp2 = std::max(maxAmp2, sqrt(c0*c0+s0*s0));
          continue;
        }
      }

      // skip all frequencies which cannot be doodson coded
      // --------------------------------------------------
      {
        Bool flag = FALSE;
        for(UInt i=0; i<6; i++)
          flag |= kn.at(i)<-5 || kn.at(i)>10;   //--> corresponding multipliers ranging from 0 to f
        if(flag)
        {
          countNDoodson++;
          maxAmp3 = std::max(maxAmp3, sqrt(c0*c0+s0*s0));
          continue;
        }
      }

      // test frequency
      if(fabs(Doodson(kn).frequency()-frequency*86400)>1e-8)
        logWarning<<Doodson(kn).name()<<" frequency difference: "<<Doodson(kn).frequency()<<" - "<<frequency*86400<<" = "<<Doodson(kn).frequency()-frequency*86400<<Log::endl;

      tgp.push_back(TideGeneratingConstituent(Doodson(kn), c0, s0));
    } // while(file.good())

    // find duplicates
    // ---------------
    UInt countDuplicate = 0;
    for(UInt i=1; i<tgp.size(); i++)
    {
      if(tgp.at(i-1)==tgp.at(i))
      {
        tgp.at(i-1).c += tgp.at(i).c;
        tgp.at(i-1).s += tgp.at(i).s;
        tgp.erase(tgp.begin()+i);
        i--;
        countDuplicate++;
      }
    }

    logInfo<<"constituents used:     "<<tgp.size()<<Log::endl;
    logInfo<<"constituents skiped:"<<Log::endl;
    logInfo<<"  permanent:           "<<countPermanent<<Log::endl;
    logInfo<<"  not degree 2:        "<<countNDegree<<Log::endl;
    logInfo<<"  not doodson codable: "<<countNDoodson<<Log::endl;
    logInfo<<"  not sun or moon:     "<<countNSunMoon<<Log::endl;
    logInfo<<"  duplicates:          "<<countDuplicate<<Log::endl;
    logInfo<<"constituents total:    "<<tgp.size()+countNDegree+countNDoodson+countNSunMoon+countPermanent+countDuplicate<<Log::endl;
    logInfo<<"  max. amplitude (degree>2):            "<<maxAmp1<<" m**2/s**2"<<Log::endl;
    logInfo<<"  max. amplitude (not sun or moon):     "<<maxAmp2<<" m**2/s**2"<<Log::endl;
    logInfo<<"  max. amplitude (not doodson codable): "<<maxAmp3<<" m**2/s**2"<<Log::endl;

    logStatus<<"save TGP <"<<outName<<">"<<Log::endl;
    writeFileTideGeneratingPotential(outName, tgp);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
