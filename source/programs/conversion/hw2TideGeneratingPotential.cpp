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
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(Hw2TideGeneratingPotential, SINGLEPROCESS, "Read tide generating potential from Hartmann and Wenzel 1995", Conversion)

/***********************************************/

void Hw2TideGeneratingPotential::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
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

    Double T = timeGPS2JC(time0); // julian centuries

    UInt countPermanent = 0;
    UInt countDegree1  = 0; Double ampDegree1  = 0.;
    UInt countDegree4  = 0; Double ampDegree4  = 0.;
    UInt countPlanet   = 0; Double ampPlanet   = 0.;
    UInt countNDoodson = 0; Double ampNDoodson = 0.;

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
      const UInt degree = String::toInt(line.substr(9, 2));

      // doodson multipliers
      std::vector<Int> kn(11);
      for(UInt i=0; i<kn.size(); i++)
        kn.at(i) = String::toInt(line.substr(11+3*i, 3));

      // frequency
      const Double frequency = 86400/3600*DEG2RAD*String::toDouble(line.substr(44, 12)); // degree/hour -> rad/day

      // amplitudes
      Double c0 = 1e-10 * String::toDouble(line.substr(56, 12));
      Double s0 = 1e-10 * String::toDouble(line.substr(68, 12));
      Double ct = 1e-10 * String::toDouble(line.substr(80, 10));
      Double st = 1e-10 * String::toDouble(line.substr(90, 10));

      // =======================================================

      c0 += ct * T;
      s0 += st * T;
      Double ampl = std::sqrt(c0*c0+s0*s0);

      // ignore permanent tides
      // ----------------------
      if(std::all_of(kn.begin(), kn.end(), [](Int n){return n==0;}))
      {
        countPermanent++;
        continue;
      }

      // ignore degree 1
      // ---------------
      if(degree == 1)
      {
        countDegree1++;
        ampDegree1 = std::max(ampDegree1, ampl);
        continue;
      }

      // use only tides from sun & moon
      // ------------------------------
      if(std::any_of(kn.begin()+6, kn.end(), [](Int n){return n!=0;}) || ((generatingBody != "SU") && (generatingBody != "MO")))
      {
        countPlanet++;
        ampPlanet = std::max(ampPlanet, ampl);
        continue;
      }

      // use only degree n=2, n=3
      // -------------------
      if(degree > 3)
      {
        countDegree4++;
        ampDegree4 = std::max(ampDegree4, ampl);
        continue;
      }

      // skip all frequencies which cannot be doodson coded
      // --------------------------------------------------
      {
        Bool flag = FALSE;
        for(UInt i=0; i<6; i++)
          flag |= (kn.at(i)+5 < -13) || (kn.at(i)+5 > 23);   //--> corresponding multipliers ranging from 0 to f
        if(flag)
        {
          countNDoodson++;
          ampNDoodson = std::max(ampNDoodson, ampl);
          continue;
        }
      }

      // test frequency
      if(std::fabs(Doodson(kn).frequency()-frequency) > 1e-8)
        logWarning<<Doodson(kn).name()<<" frequency difference: "<<Doodson(kn).frequency()<<" - "<<frequency<<" = "<<Doodson(kn).frequency()-frequency<<Log::endl;

      tgp.push_back(TideGeneratingConstituent(Doodson(kn), degree, c0, s0));
    } // while(file.good())

    // find duplicates
    // ---------------
    std::stable_sort(tgp.begin(), tgp.end());
    UInt countDuplicate = 0;
    for(UInt i=1; i<tgp.size(); i++)
      if(tgp.at(i-1) == tgp.at(i))
      {
        tgp.at(i-1).c += tgp.at(i).c;
        tgp.at(i-1).s += tgp.at(i).s;
        tgp.erase(tgp.begin()+i);
        i--;
        countDuplicate++;
      }

    logInfo<<"constituents used:     "<<tgp.size()<<Log::endl;
    logInfo<<"constituents skipped:"<<Log::endl;
    logInfo<<"  permanent:           "<<countPermanent<<Log::endl;
    logInfo<<"  degree 1:            "<<countDegree1<<Log::endl;
    logInfo<<"  degree >3:           "<<countDegree4<<Log::endl;
    logInfo<<"  not doodson codable: "<<countNDoodson<<Log::endl;
    logInfo<<"  not sun or moon:     "<<countPlanet<<Log::endl;
    logInfo<<"  duplicates:          "<<countDuplicate<<Log::endl;
    logInfo<<"constituents total:    "<<tgp.size()+countDegree1+countDegree4+countNDoodson+countPlanet+countPermanent+countDuplicate<<Log::endl;
    logInfo<<"  max. amplitude (degree 1):            "<<ampDegree1 <<" m^2/s^2"<<Log::endl;
    logInfo<<"  max. amplitude (degree >3):           "<<ampDegree4 <<" m^2/s^2"<<Log::endl;
    logInfo<<"  max. amplitude (not sun or moon):     "<<ampPlanet  <<" m^2/s^2"<<Log::endl;
    logInfo<<"  max. amplitude (not doodson codable): "<<ampNDoodson<<" m^2/s^2"<<Log::endl;

    logStatus<<"save TGP <"<<outName<<">"<<Log::endl;
    writeFileTideGeneratingPotential(outName, tgp);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
