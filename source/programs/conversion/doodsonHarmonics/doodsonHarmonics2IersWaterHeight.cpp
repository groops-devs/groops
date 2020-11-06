/***********************************************/
/**
* @file doodsonHarmonics2IersWaterHeight.cpp
*
* @brief Convert doodson harmonics to IERS conventions according to FES2004.
* cf. ftp://tai.bipm.org/iers/conv2010/chapter6/tidemodels/fes2004.dat
*
* @author Daniel Rieser
* @date 2012-03-05
*
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
Convert doodson harmonics to IERS conventions according to FES2004.
cf. \url{ftp://tai.bipm.org/iers/conv2010/chapter6/tidemodels/fes2004.dat}.
)";

/***********************************************/

#include "programs/program.h"
#include "base/doodson.h"
#include "inputOutput/file.h"
#include "files/fileDoodsonHarmonic.h"
#include "files/fileTideGeneratingPotential.h"
#include "classes/kernel/kernel.h"

/***** CLASS ***********************************/

/** @brief Convert doodson harmonics to IERS conventions according to FES2004.
* cf. ftp://tai.bipm.org/iers/conv2010/chapter6/tidemodels/fes2004.dat
* @ingroup programsConversionGroup */
class DoodsonHarmonics2IersWaterHeight
{
public:
  void run(Config &config);
};

GROOPS_REGISTER_PROGRAM(DoodsonHarmonics2IersWaterHeight, SINGLEPROCESS, "convert doodson harmonics file to IERS", Conversion, DoodsonHarmonics)

/***********************************************/

void DoodsonHarmonics2IersWaterHeight::run(Config &config)
{
  try
  {
    FileName    fileNameOut, fileNameIn, fileNameTGP;
    KernelPtr   kernel;
    Double      factor;
    UInt        minDegree, maxDegree = INFINITYDEGREE;
    std::vector<std::string> headers;

    readConfig(config, "outputfile",                       fileNameOut, Config::MUSTSET,  "","according to IERS2010, chapter 6.3.2, footnote 7");
    readConfig(config, "inputfileDoodsonHarmoncis",        fileNameIn,  Config::MUSTSET,  "", "");
    readConfig(config, "inputfileTideGeneratingPotential", fileNameTGP, Config::MUSTSET,  "{groopsDataDir}/tides/generatingTide_HW95.txt", "to compute Xi phase correction");
    readConfig(config, "header",                           headers,     Config::MUSTSET,  "Ocean tide model: FES2014b up to (100,100) in cm", "info for output header");
    readConfig(config, "kernel",                           kernel,      Config::MUSTSET,  "waterHeight", "data type of output values");
    readConfig(config, "factor",                           factor,      Config::DEFAULT,  "100", "e.g. from [m] to [cm]");
    readConfig(config, "minDegree",                        minDegree,   Config::DEFAULT,  "0",   "");
    readConfig(config, "maxDegree",                        maxDegree,   Config::OPTIONAL, "",    "");
    if(isCreateSchema(config)) return;

    // ==============================

    // read tides
    // ----------
    logStatus<<"read tides from <"<<fileNameIn<<">"<<Log::endl;
    DoodsonHarmonic d;
    readFileDoodsonHarmonic(fileNameIn, d);
    maxDegree = std::min(maxDegree, d.cnmCos.at(0).rows()-1);

    // read tide generating potential (TGP)
    // ------------------------------------
    logStatus<<"read tide generating potential <"<<fileNameTGP<<">"<<Log::endl;
    TideGeneratingPotential tgp;
    readFileTideGeneratingPotential(fileNameTGP, tgp);

    // Ocean tide model: FES2004 normalized model (fev. 2004) up to (100,100) in cm
    // (long period from FES2002 up to (50,50) + equilibrium Om1/Om2, atmospheric tide NOT included)
    // Doodson Darw  n   m    Csin+     Ccos+       Csin-     Ccos-       C+   eps+      C-   eps-
    // 255.555 M2    2   1 -0.991591 -0.360545   -0.253804  0.123315   1.0551 250.019 0.2822 295.914
    logStatus<<"write tides to <"<<fileNameOut<<">"<<Log::endl;
    OutFile file(fileNameOut);
    for(const auto &header : headers)
      file<<header<<std::endl;
    file<<"Doodson Darw  n   m    Ccos+     Csin+       Ccos-     Csin-       C+   eps+      C-   eps-"<<std::endl;
    const Vector kn = factor * d.GM/d.R * kernel->inverseCoefficients(Vector3d(d.R,0,0), maxDegree); // convert units
    for(UInt i=0; i<d.doodson.size(); i++)
    {
      const Double xi = tgp.xi(d.doodson.at(i));
      logStatus<<d.doodson.at(i).code()<<" "<<d.doodson.at(i).name()<<std::string(6ul-std::min(static_cast<std::size_t>(6), d.doodson.at(i).name().size()), ' ')<<"xi = "<<RAD2DEG*xi<<Log::endl;
      for(UInt m=0; m<=maxDegree; m++)
        for(UInt n=std::max(minDegree,m); n<=maxDegree; n++)
        {
          // iers conventions 2010, eq. (6.15)
          const Double cPlus  = 0.5 * (d.cnmCos.at(i)(n,m) - d.snmSin.at(i)(n,m));
          const Double cMinus = 0.5 * (d.cnmCos.at(i)(n,m) + d.snmSin.at(i)(n,m));
          const Double sPlus  = 0.5 * (d.cnmSin.at(i)(n,m) + d.snmCos.at(i)(n,m));
          const Double sMinus = 0.5 * (d.cnmSin.at(i)(n,m) - d.snmCos.at(i)(n,m));

          // iers conventions 2010, eq. (6.20) and (6.21)
          const Double cCosPlus  = kn(n) * (cPlus  * std::cos(xi) - sPlus  * std::sin(xi));
          const Double cCosMinus = kn(n) * (cMinus * std::cos(xi) - sMinus * std::sin(xi));
          const Double cSinPlus  = kn(n) * (cPlus  * std::sin(xi) + sPlus  * std::cos(xi));
          const Double cSinMinus = kn(n) * (cMinus * std::sin(xi) + sMinus * std::cos(xi));

          // iers conventions 2010, eq. (6.20)
          const Double cBarPlus  = std::sqrt(cCosPlus  * cCosPlus  + cSinPlus  * cSinPlus);
          const Double cBarMinus = std::sqrt(cCosMinus * cCosMinus + cSinMinus * cSinMinus);
          const Double epsPlus   = RAD2DEG * std::fmod(std::atan2(cCosPlus,  cSinPlus)  + 2*PI, 2*PI);
          const Double epsMinus  = RAD2DEG * std::fmod(std::atan2(cCosMinus, cSinMinus) + 2*PI, 2*PI);

          file<<d.doodson.at(i).code()<<" "<<d.doodson.at(i).name()<<std::string(4ul-std::min(static_cast<std::size_t>(3), d.doodson.at(i).name().size()), ' ');
          file<<n%"%3i "s<<m%"%3i "s;
          file<<cCosPlus%"%9.6f "s<<cSinPlus%"%9.6f   "s<<cCosMinus%"%9.6f "s<<cSinMinus%"%9.6f   "s;
          file<<cBarPlus%"%6.4f "s<<epsPlus%"%7.3f "s<<cBarMinus%"%6.4f "s<<epsMinus%"%7.3f"s<<std::endl;
        }
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
