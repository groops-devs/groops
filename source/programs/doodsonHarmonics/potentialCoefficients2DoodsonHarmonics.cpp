/***********************************************/
/**
* @file potentialCoefficients2DoodsonHarmonics.cpp
*
* @brief DoodsonHarmonic file from a list of cos/sin potentialCoefficients.
*
* @author Torsten Mayer-Guerr
* @date 2007-10-24
*
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
Create a \file{DoodsonHarmonic file}{doodsonHarmonic} from a list of
cos/sin \file{potentialCoefficients}{potentialCoefficients} for given \config{doodson}
(Doodson number or Darwin´s name, e.g. 255.555 or M2) tidal constituents.
If \config{applyXi} the Doodson-Warburg phase correction (see IERS conventions) is applied before.
)";

/***********************************************/

#include "programs/program.h"
#include "base/doodson.h"
#include "files/fileSphericalHarmonics.h"
#include "files/fileDoodsonHarmonic.h"
#include "files/fileTideGeneratingPotential.h"

/***** CLASS ***********************************/

/** @brief DoodsonHarmonic file from a list of cos/sin potentialCoefficients.
* @ingroup programsGroup */
class PotentialCoefficients2DoodsonHarmonics
{
public:
  void run(Config &config);

  class Constituent
  {
    public:
    Doodson  doodson;
    FileName cosName, sinName;
  };
};

GROOPS_REGISTER_PROGRAM(PotentialCoefficients2DoodsonHarmonics, SINGLEPROCESS, "doodsonHarmonic file from cos/sin potentialCoefficients.", DoodsonHarmonics, PotentialCoefficients)

/***********************************************/

template<> Bool readConfig(Config &config, const std::string &name, PotentialCoefficients2DoodsonHarmonics::Constituent &var, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation)
{
  if(!readConfigSequence(config, name, mustSet, defaultValue, annotation))
    return FALSE;
  readConfig(config, "doodson",                           var.doodson, Config::MUSTSET, "", "");
  readConfig(config, "inputfileCosPotentialCoefficients", var.cosName, Config::MUSTSET, "", "");
  readConfig(config, "inputfileSinPotentialCoefficients", var.sinName, Config::MUSTSET, "", "");
  endSequence(config);
  return TRUE;
}

/***********************************************/

void PotentialCoefficients2DoodsonHarmonics::run(Config &config)
{
  try
  {
    FileName fileNameOut;
    FileName fileNameTGP;
    std::vector<Constituent> constituent;
    UInt     minDegree, maxDegree = INFINITYDEGREE;
    Double   GM, R;
    Bool     applyXi;

    readConfig(config, "outputfileDoodsonHarmonics",       fileNameOut, Config::MUSTSET,   "",  "");
    readConfig(config, "inputfileTideGeneratingPotential", fileNameTGP, Config::OPTIONAL, "{groopsDataDir}/tides/generatingTide_HW95.txt", "to compute Xi phase correction");
    readConfig(config, "constituent",                      constituent, Config::MUSTSET,   "",  "");
    readConfig(config, "minDegree",                        minDegree,   Config::DEFAULT,   "0", "");
    readConfig(config, "maxDegree",                        maxDegree,   Config::OPTIONAL, "",  "");
    readConfig(config, "GM",                               GM,          Config::DEFAULT,  STRING_DEFAULT_GM, "Geocentric gravitational constant");
    readConfig(config, "R",                                R,           Config::DEFAULT,  STRING_DEFAULT_R,  "reference radius");
    readConfig(config, "applyXi",                          applyXi,     Config::DEFAULT,  "1",               "apply Doodson-Warburg phase correction (see IERS conventions)");
    if(isCreateSchema(config)) return;

    // read tide generating potential (TGP)
    // ------------------------------------
    TideGeneratingPotential tgp;
    if(applyXi)
    {
      if(fileNameTGP.empty())
        throw(Exception("Need TideGeneratingPotential to compute xi phase correction"));
      logStatus<<"read tide generating potential <"<<fileNameTGP<<">"<<Log::endl;
      readFileTideGeneratingPotential(fileNameTGP, tgp);
    }

    // read coefficients
    // -----------------
    logStatus<<"read potential coefficients"<<Log::endl;
    std::vector<Doodson> doodson(constituent.size());
    std::vector<Matrix>  cnmCos(constituent.size());
    std::vector<Matrix>  snmCos(constituent.size());
    std::vector<Matrix>  cnmSin(constituent.size());
    std::vector<Matrix>  snmSin(constituent.size());
    for(UInt i=0; i<constituent.size(); i++)
    {
      SphericalHarmonics harmCos, harmSin;
      readFileSphericalHarmonics(constituent.at(i).cosName, harmCos);
      readFileSphericalHarmonics(constituent.at(i).sinName, harmSin);
      harmCos = harmCos.get(maxDegree, minDegree, GM, R);
      harmSin = harmSin.get(maxDegree, minDegree, GM, R);

      Double xi = 0.;
      if(applyXi)
        xi = tgp.xi(constituent.at(i).doodson);

      doodson.at(i) =  constituent.at(i).doodson;
      cnmCos.at(i)  =  harmCos.cnm() * cos(xi) + harmSin.cnm() * sin(xi);
      snmCos.at(i)  =  harmCos.snm() * cos(xi) + harmSin.snm() * sin(xi);
      cnmSin.at(i)  = -harmCos.cnm() * sin(xi) + harmSin.cnm() * cos(xi);
      snmSin.at(i)  = -harmCos.snm() * sin(xi) + harmSin.snm() * cos(xi);
    }

    // write doodson harmonics
    // -----------------------
    logStatus<<"writing doodson harmonics <"<<fileNameOut<<">"<<Log::endl;
    writeFileDoodsonHarmonic(fileNameOut, DoodsonHarmonic(GM, R, doodson, cnmCos, snmCos, cnmSin, snmSin));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
