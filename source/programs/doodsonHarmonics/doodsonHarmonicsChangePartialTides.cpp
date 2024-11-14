/***********************************************/
/**
* @file doodsonHarmonicsChangePartialTides.cpp
*
* @brief Change partial tides in doodsonHarmonics file.
*
* @author Torsten Mayer-Guerr
* @date 2024-11-07
*
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
Reads a file \configFile{inputfileDoodsonHarmonic}{doodsonHarmonic} and write it to
\configFile{outputfileDoodsonHarmonics}{doodsonHarmonic}. If set the spherical harmonics
expansion is limited in the range between \config{minDegree} and \config{maxDegree} inclusivly.
The \configClass{useDoodson}{doodson} and \configClass{ignoreDoodson}{doodson} can be used
to filter the partial types that will be exported.
Additional partial tides can be interpolated using the file \configFile{inputfileAdmittance}{admittance}.
)";

/***********************************************/

#include "programs/program.h"
#include "base/doodson.h"
#include "base/sphericalHarmonics.h"
#include "files/fileAdmittance.h"
#include "files/fileDoodsonHarmonic.h"

/***** CLASS ***********************************/

/** @brief Change partial tides in doodsonHarmonics file.
* @ingroup programsGroup */
class DoodsonHarmonicsChangePartialTides
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(DoodsonHarmonicsChangePartialTides, SINGLEPROCESS, "change partial tides in doodsonHarmonics file", DoodsonHarmonics)

/***********************************************/

void DoodsonHarmonicsChangePartialTides::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName fileNameOut, fileNameIn;
    FileName fileNameAdmit;
    std::vector<Doodson> useDoodson, ignoreDoodson;
    UInt     minDegree, maxDegree = INFINITYDEGREE;
    Double   GM = 0, R = 0;

    readConfig(config, "outputfileDoodsonHarmonics", fileNameOut,   Config::MUSTSET,  "", "");
    readConfig(config, "inputfileDoodsonHarmonics",  fileNameIn,    Config::MUSTSET,  "", "");
    readConfig(config, "inputfileAdmittance",        fileNameAdmit, Config::OPTIONAL, "", "interpolation of minor constituents");
    readConfig(config, "useDoodson",                 useDoodson,    Config::OPTIONAL, "", "use only these partial tides (additional tides will be interpolated)");
    readConfig(config, "ignoreDoodson",              ignoreDoodson, Config::OPTIONAL, "", "ignore these partial tides");
    readConfig(config, "minDegree",                  minDegree,     Config::DEFAULT,  "0", "");
    readConfig(config, "maxDegree",                  maxDegree,     Config::OPTIONAL, "",  "");
    readConfig(config, "GM",                         GM,            Config::OPTIONAL, STRING_DEFAULT_GM, "Geocentric gravitational constant");
    readConfig(config, "R",                          R,             Config::OPTIONAL, STRING_DEFAULT_R,  "reference radius");
    if(isCreateSchema(config)) return;

    logStatus<<"read doodson harmonics file <"<<fileNameIn<<">"<<Log::endl;
    DoodsonHarmonic d;
    readFileDoodsonHarmonic(fileNameIn, d);
    for(UInt i=0; i<d.doodson.size(); i++)
    {
      SphericalHarmonics harmCos = SphericalHarmonics(d.GM, d.R, d.cnmCos.at(i), d.snmCos.at(i)).get(maxDegree, minDegree, GM, R);
      SphericalHarmonics harmSin = SphericalHarmonics(d.GM, d.R, d.cnmSin.at(i), d.snmSin.at(i)).get(maxDegree, minDegree, GM, R);
      d.cnmCos.at(i) = harmCos.cnm();
      d.snmCos.at(i) = harmCos.snm();
      d.cnmSin.at(i) = harmSin.cnm();
      d.snmSin.at(i) = harmSin.snm();
    }

    Admittance admit;
    admit.doodsonMajor = admit.doodsonMinor = d.doodson;
    admit.admittance   = identityMatrix(d.doodson.size());
    if(!fileNameAdmit.empty())
    {
      logStatus<<"read admittance file <"<<fileNameAdmit<<">"<<Log::endl;
      readFileAdmittance(fileNameAdmit, admit);
    }

    logStatus<<"calculate partial tides"<<Log::endl;
    if(!useDoodson.size())
      useDoodson = d.doodson;
    std::sort(useDoodson.begin(), useDoodson.end());
    std::vector<Doodson> doodsonNew;
    std::vector<Matrix>  cnmCos, snmCos;
    std::vector<Matrix>  cnmSin, snmSin;
    for(auto &dood : useDoodson)
      if(std::find(ignoreDoodson.begin(), ignoreDoodson.end(), dood) == ignoreDoodson.end())
      {
        UInt idx = std::distance(admit.doodsonMinor.begin(), std::find(admit.doodsonMinor.begin(), admit.doodsonMinor.end(), dood));
        if(idx >= admit.doodsonMinor.size())
        {
          logWarning<<"Unable to interpolate "<<dood.code()<<" ("<<dood.name()<<"), skipped."<<Log::endl;
          continue;
        }
        cnmCos.push_back(Matrix(d.cnmCos.front().rows(), Matrix::SYMMETRIC, Matrix::LOWER));
        snmCos.push_back(Matrix(d.cnmCos.front().rows(), Matrix::SYMMETRIC, Matrix::LOWER));
        cnmSin.push_back(Matrix(d.cnmCos.front().rows(), Matrix::SYMMETRIC, Matrix::LOWER));
        snmSin.push_back(Matrix(d.cnmCos.front().rows(), Matrix::SYMMETRIC, Matrix::LOWER));
        for(UInt k=0; k<admit.doodsonMajor.size(); k++)
        {
          cnmCos.back() += admit.admittance(k, idx) * d.cnmCos.at(k);
          snmCos.back() += admit.admittance(k, idx) * d.snmCos.at(k);
          cnmSin.back() += admit.admittance(k, idx) * d.cnmSin.at(k);
          snmSin.back() += admit.admittance(k, idx) * d.snmSin.at(k);
        }
        doodsonNew.push_back(dood);
      }

    logStatus<<"writing doodson harmonics <"<<fileNameOut<<">"<<Log::endl;
    writeFileDoodsonHarmonic(fileNameOut, DoodsonHarmonic(d.GM, d.R, doodsonNew, cnmCos, snmCos, cnmSin, snmSin));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
