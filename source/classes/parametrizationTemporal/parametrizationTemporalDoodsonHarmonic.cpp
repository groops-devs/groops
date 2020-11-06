/***********************************************/
/**
* @file parametrizationTemporalDoodsonHarmonic.cpp
*
* @brief Tidal variations.
* @see ParametrizationTemporal
*
* @author Torsten Mayer-Guerr
* @date 2011-12-28
*
*/
/***********************************************/

#include "base/import.h"
#include "config/config.h"
#include "base/doodson.h"
#include "files/fileAdmittance.h"
#include "classes/parametrizationTemporal/parametrizationTemporal.h"
#include "parametrizationTemporalDoodsonHarmonic.h"

/***********************************************/

ParametrizationTemporalDoodsonHarmonic::ParametrizationTemporalDoodsonHarmonic(Config &config)
{
  try
  {
    FileName admittanceName;

    readConfig(config, "doodson",             majorDoodson,   Config::MUSTSET,  "", "code number (e.g. 255.555) or darwin name (e.g. M2)");
    readConfig(config, "inputfileAdmittance", admittanceName, Config::OPTIONAL, "", "interpolation of minor constituents");
    if(isCreateSchema(config)) return;

    // read admittace file
    // -------------------
    if(!admittanceName.empty())
    {
      Admittance admit;
      readFileAdmittance(admittanceName, admit);

      admittance = Matrix(majorDoodson.size(), admit.doodsonMinor.size());
      for(UInt i=0; i<majorDoodson.size(); i++)
      {
        Bool found = FALSE;
        for(UInt k=0; k<admit.doodsonMajor.size(); k++)
          if(majorDoodson.at(i) == admit.doodsonMajor.at(k))
          {
            found = TRUE;
            copy(admit.admittance.row(k), admittance.row(i));
            break;
          }
        if(!found)
          throw(Exception(majorDoodson.at(i).name()+" not found in admittance file"));
      }

      // Matrix with Doodson multiplicators
      doodsonMatrix = Doodson::matrix(admit.doodsonMinor);
    }
    else
    {
      // Matrix with Doodson multiplicators
      doodsonMatrix = Doodson::matrix(majorDoodson);
    }

    _parameterCount = 2 * majorDoodson.size();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void ParametrizationTemporalDoodsonHarmonic::parameterName(std::vector<ParameterName> &name) const
{
  try
  {
    for(UInt i=0; i<majorDoodson.size(); i++)
    {
      name.push_back(ParameterName("", "", "doodson.cos("+majorDoodson.at(i).name()+")"));
      name.push_back(ParameterName("", "", "doodson.sin("+majorDoodson.at(i).name()+")"));
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void ParametrizationTemporalDoodsonHarmonic::factors(const Time &time, UInt startIndex, std::vector<UInt> &index, std::vector<Double> &factor) const
{
  try
  {
    Vector thetaf = doodsonMatrix * Doodson::arguments(time);

    Matrix csMinor(thetaf.rows(), 2);
    for(UInt i=0; i<thetaf.rows(); i++)
    {
      csMinor(i,0) = cos(thetaf(i));
      csMinor(i,1) = sin(thetaf(i));
    }
    Matrix csMajor = ((admittance.size()) ? (admittance * csMinor) : csMinor);

    for(UInt i=0; i<csMajor.rows(); i++)
    {
      index.push_back(startIndex++);
      index.push_back(startIndex++);
      factor.push_back(csMajor(i,0));
      factor.push_back(csMajor(i,1));
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
