/***********************************************/
/**
* @file parametrizationGnssAntenna.cpp
*
* @brief Parametrization of antenna center variations.
*
* @author Torsten Mayer-Guerr
* @date 2012-11-08
*
*/
/***********************************************/

#define DOCSTRING_ParametrizationGnssAntenna

#include "base/import.h"
#include "base/gnssType.h"
#include "config/configRegister.h"
#include "parametrizationGnssAntennaSphericalHarmonics.h"
#include "parametrizationGnssAntennaRadialBasis.h"
#include "parametrizationGnssAntennaCenter.h"
#include "parametrizationGnssAntenna.h"

/***********************************************/

GROOPS_REGISTER_CLASS(ParametrizationGnssAntenna, "parametrizationGnssAntennaType",
                      ParametrizationGnssAntennaCenter,
                      ParametrizationGnssAntennaSphericalHarmonics,
                      ParametrizationGnssAntennaRadialBasis)

GROOPS_RENAMED_CLASS(gnssAntennaRepresentationType, parametrizationGnssAntennaType, date2time(2020, 6, 3))

GROOPS_READCONFIG_UNBOUNDED_CLASS(ParametrizationGnssAntenna, "parametrizationGnssAntennaType")

/***********************************************/

ParametrizationGnssAntenna::ParametrizationGnssAntenna(Config &config, const std::string &name)
{
  try
  {
    std::string type;
    while(readConfigChoice(config, name, type, Config::OPTIONAL, "", "antenna center variations"))
    {
      if(readConfigChoiceElement(config, "center",  type, ""))
        base.push_back(new ParametrizationGnssAntennaCenter(config));
      if(readConfigChoiceElement(config, "sphericalHarmonics",  type, ""))
        base.push_back(new ParametrizationGnssAntennaSphericalHarmonics(config));
      if(readConfigChoiceElement(config, "radialBasis",  type, ""))
        base.push_back(new ParametrizationGnssAntennaRadialBasis(config));

      endChoice(config);
      if(isCreateSchema(config))
        return;
    }

    parameterCount_ = 0;
    index.resize(base.size());
    for(UInt i=0; i<base.size(); i++)
    {
      index.at(i) = parameterCount_;
      parameterCount_ += base.at(i)->parameterCount();
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

ParametrizationGnssAntenna::~ParametrizationGnssAntenna()
{
  for(UInt i=0; i<base.size(); i++)
    delete base.at(i);
}

/***********************************************/

void ParametrizationGnssAntenna::parameterName(std::vector<ParameterName> &name) const
{
  for(UInt i=0; i<base.size(); i++)
    base.at(i)->parameterName(name);
}

/***********************************************/

Matrix ParametrizationGnssAntenna::designMatrix(Angle azimut, Angle elevation)
{
  try
  {
    Matrix A(1, parameterCount_);
    for(UInt i=0; i<base.size(); i++)
      base.at(i)->designMatrix(azimut, elevation, A.column(index.at(i), base.at(i)->parameterCount()));
    return A;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
