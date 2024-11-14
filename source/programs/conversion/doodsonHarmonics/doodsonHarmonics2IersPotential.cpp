/***********************************************/
/**
* @file doodsonHarmonics2IersPotential.cpp
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
#include "classes/kernel/kernel.h"

/***** CLASS ***********************************/

/** @brief Convert doodson harmonics to IERS conventions according to FES2004.
* cf. ftp://tai.bipm.org/iers/conv2010/chapter6/tidemodels/fes2004.dat
* @ingroup programsConversionGroup */
class DoodsonHarmonics2IersPotential
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(DoodsonHarmonics2IersPotential, SINGLEPROCESS, "convert doodson harmonics file to IERS", Conversion, DoodsonHarmonics)

/***********************************************/

void DoodsonHarmonics2IersPotential::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName    fileNameOut, fileNameIn;
    Double      factor;
    UInt        minDegree, maxDegree = INFINITYDEGREE;
    std::vector<std::string> headers;

    readConfig(config, "outputfile",                       fileNameOut, Config::MUSTSET,  "","according to IERS2010, chapter 6.3.2, footnote 7");
    readConfig(config, "inputfileDoodsonHarmoncis",        fileNameIn,  Config::MUSTSET,  "", "");
    readConfig(config, "header",                           headers,     Config::MUSTSET,  R"|(["Coefficients to compute variations in normalized Stokes coefficients (unit = 10^-11)", "Ocean tide model: FES2014b up to (100,100) in cm"])|", "info for output header");
    readConfig(config, "factor",                           factor,      Config::DEFAULT,  "1e11", "");
    readConfig(config, "minDegree",                        minDegree,   Config::DEFAULT,  "1",    "");
    readConfig(config, "maxDegree",                        maxDegree,   Config::OPTIONAL, "",     "");
    if(isCreateSchema(config)) return;

    // ==============================

    // read tides
    // ----------
    logStatus<<"read tides from <"<<fileNameIn<<">"<<Log::endl;
    DoodsonHarmonic d;
    readFileDoodsonHarmonic(fileNameIn, d);
    maxDegree = std::min(maxDegree, d.cnmCos.at(0).rows()-1);

    // Coefficients to compute variations in normalized Stokes coefficients (unit = 10^-11)
    // Ocean tide model: FES2004 normalized model (fev. 2004) up to (100,100)
    // (long period from FES2002 up to (50,50) + equilibrium Om1/Om2, atmospheric tide NOT included)
    // Doodson Darw  l   m    DelC+     DelS+       DelC-     DelS-
    // 255.555 M2    2   1 -12.07164  -4.38919    -3.09008   1.50139
    logStatus<<"write tides to <"<<fileNameOut<<">"<<Log::endl;
    OutFile file(fileNameOut);
    for(const auto &header : headers)
      file<<header<<std::endl;
    file<<"Doodson Darw  l   m    DelC+     DelS+       DelC-     DelS-"<<std::endl;
    for(UInt i=0; i<d.doodson.size(); i++)
    {
      for(UInt m=0; m<=maxDegree; m++)
        for(UInt n=std::max(minDegree,m); n<=maxDegree; n++)
        {
          const Double cPlus  = 0.5 * factor * (d.cnmCos.at(i)(n,m) - d.snmSin.at(i)(n,m));
          const Double cMinus = 0.5 * factor * (d.cnmCos.at(i)(n,m) + d.snmSin.at(i)(n,m));
          const Double sPlus  = 0.5 * factor * (d.cnmSin.at(i)(n,m) + d.snmCos.at(i)(n,m));
          const Double sMinus = 0.5 * factor * (d.cnmSin.at(i)(n,m) - d.snmCos.at(i)(n,m));

          file<<d.doodson.at(i).code()<<" "<<d.doodson.at(i).name()<<std::string(4ul-std::min(static_cast<std::size_t>(3), d.doodson.at(i).name().size()), ' ');
          file<<n%"%3i "s<<m%"%3i "s;
          file<<cPlus%"%9.5f "s<<sPlus%"%9.5f   "s<<cMinus%"%9.5f "s<<sMinus%"%9.5f   "s<<std::endl;
        }
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
