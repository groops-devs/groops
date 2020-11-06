/***********************************************/
/**
* @file gravityfieldReplacePotentialCoefficients.cpp
*
* @brief Replace single potential coefficients in a gravity field
*
* @author Andreas Kvas
* @author Torsten Mayer-Guerr
* @date 2015-11-17
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
Replaces single potential coefficients in a gravity field.
Both \configClass{gravityfield}{gravityfieldType}
and \configClass{gravityfieldReplacement}{gravityfieldType} are evalutated
at \config{time} and converted to spherical harmonic coefficients.
Single \config{coefficients} are then replaced in \configClass{gravityfield}{gravityfieldType}
by the values from \configClass{gravityfieldReplacement}{gravityfieldType}
and the result is written to \configFile{outputfilePotentialCoefficients}{potentialCoefficients}
from \config{minDegree} to \config{maxDegree},
)";

/***********************************************/

#include "programs/program.h"
#include "files/fileSphericalHarmonics.h"
#include "classes/gravityfield/gravityfield.h"

/***** CLASS ***********************************/

/** @brief Replace single potential coefficients in a gravity field.
* @ingroup programsGroup */
class GravityfieldReplacePotentialCoefficients
{
public:
  class Coeff
  {
  public:
    Bool isCnm;
    UInt degree, order;
  };

  void run(Config &config);
};

GROOPS_REGISTER_PROGRAM(GravityfieldReplacePotentialCoefficients, SINGLEPROCESS, "Replace single potential coefficients in a gravity field", Gravityfield)
GROOPS_RENAMED_PROGRAM(GravityfieldReplaceC20, GravityfieldReplacePotentialCoefficients, date2time(2020, 5, 24))

/***********************************************/

template<> Bool readConfig(Config &config, const std::string &name, GravityfieldReplacePotentialCoefficients::Coeff &var, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation)
{
  std::string choice;
  if(!readConfigChoice(config, name, choice, mustSet, defaultValue, annotation))
    return FALSE;
  if(readConfigChoiceElement(config, "cnm", choice, ""))
  {
    var.isCnm = TRUE;
    readConfig(config, "degree", var.degree, Config::MUSTSET, "", "");
    readConfig(config, "order",  var.order,  Config::MUSTSET, "", "");
  }
  if(readConfigChoiceElement(config, "snm", choice, ""))
  {
    var.isCnm = FALSE;
    readConfig(config, "degree", var.degree, Config::MUSTSET, "", "");
    readConfig(config, "order",  var.order,  Config::MUSTSET, "", "");
  }
  endChoice(config);
  return TRUE;
}

/***********************************************/

void GravityfieldReplacePotentialCoefficients::run(Config &config)
{
  try
  {
    FileName           fileNameCoeff;
    GravityfieldPtr    gravityfield1, gravityfield2;
    std::vector<Coeff> coefficients;
    UInt               minDegree, maxDegree = INFINITYDEGREE;
    Double             GM, R;
    Time               time;

    readConfig(config, "outputfilePotentialCoefficients", fileNameCoeff, Config::MUSTSET,  "", "");
    readConfig(config, "gravityfield",                    gravityfield1, Config::MUSTSET,  "", "single coefficients are replaced by the other gravityfield");
    readConfig(config, "gravityfieldReplacement",         gravityfield2, Config::MUSTSET,  "", "contains the coefficients for replacement");
    readConfig(config, "coefficients",                    coefficients,  Config::MUSTSET,  "", "");
    readConfig(config, "minDegree",                       minDegree,     Config::DEFAULT,  "0", "");
    readConfig(config, "maxDegree",                       maxDegree,     Config::OPTIONAL, "",  "");
    readConfig(config, "GM",                              GM,            Config::DEFAULT,  STRING_DEFAULT_GM, "Geocentric gravitational constant");
    readConfig(config, "R",                               R,             Config::DEFAULT,  STRING_DEFAULT_R,  "reference radius");
    readConfig(config, "time",                            time,          Config::OPTIONAL, "", "at this time the gravity field will be evaluated");
    if(isCreateSchema(config)) return;

    logStatus<<"create spherical harmonics"<<Log::endl;
    SphericalHarmonics harm1 = gravityfield1->sphericalHarmonics(time, maxDegree, minDegree, GM, R);
    SphericalHarmonics harm2 = gravityfield2->sphericalHarmonics(time, maxDegree, minDegree, GM, R);

    for(auto &coeff : coefficients)
    {
      if(coeff.isCnm)
      {
        harm1.cnm()(coeff.degree, coeff.order) = harm2.cnm()(coeff.degree, coeff.order);
        if(harm1.sigma2cnm().size())
          harm1.sigma2cnm()(coeff.degree, coeff.order) = 0;
        if(harm1.sigma2cnm().size() && harm2.sigma2cnm().size())
          harm1.sigma2cnm()(coeff.degree, coeff.order) = harm2.sigma2cnm()(coeff.degree, coeff.order);
      }
      else
      {
        harm1.snm()(coeff.degree, coeff.order) = harm2.snm()(coeff.degree, coeff.order);
        if(harm1.sigma2snm().size())
          harm1.sigma2snm()(coeff.degree, coeff.order) = 0;
        if(harm1.sigma2snm().size() && harm2.sigma2snm().size())
          harm1.sigma2snm()(coeff.degree, coeff.order) = harm2.sigma2snm()(coeff.degree, coeff.order);
      }
    }

    logStatus<<"writing potential coefficients to file <"<<fileNameCoeff<<">"<<Log::endl;
    writeFileSphericalHarmonics(fileNameCoeff, harm1);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
