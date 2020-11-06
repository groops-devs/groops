/***********************************************/
/**
* @file parametrizationTemporal.cpp
*
* @brief Parametrization of temporal variations.
*
* @author Torsten Mayer-Guerr
* @date 2011-12-28
*
*/
/***********************************************/

#define DOCSTRING_ParametrizationTemporal

#include "base/import.h"
#include "config/configRegister.h"
#include "classes/parametrizationTemporal/parametrizationTemporalConstant.h"
#include "classes/parametrizationTemporal/parametrizationTemporalTrend.h"
#include "classes/parametrizationTemporal/parametrizationTemporalSplines.h"
#include "classes/parametrizationTemporal/parametrizationTemporalPolynomial.h"
#include "classes/parametrizationTemporal/parametrizationTemporalOscillation.h"
#include "classes/parametrizationTemporal/parametrizationTemporalFourier.h"
#include "classes/parametrizationTemporal/parametrizationTemporalDoodsonHarmonic.h"
#include "classes/parametrizationTemporal/parametrizationTemporal.h"

/***********************************************/

GROOPS_REGISTER_CLASS(ParametrizationTemporal, "parametrizationTemporalType",
                      ParametrizationTemporalConstant,
                      ParametrizationTemporalTrend,
                      ParametrizationTemporalSplines,
                      ParametrizationTemporalPolynomial,
                      ParametrizationTemporalOscillation,
                      ParametrizationTemporalFourier,
                      ParametrizationTemporalDoodsonHarmonic)

GROOPS_RENAMED_CLASS(temporalRepresentationType, parametrizationTemporalType, date2time(2020, 6, 3))

GROOPS_READCONFIG_UNBOUNDED_CLASS(ParametrizationTemporal, "parametrizationTemporalType")

/***********************************************/

ParametrizationTemporal::ParametrizationTemporal(Config &config, const std::string &name)
{
  try
  {
    std::string type;
    while(readConfigChoice(config, name, type, Config::OPTIONAL, "", "Parametrization of temporal variations"))
    {
      if(readConfigChoiceElement(config, "constant",        type, "Constant in time"))
        representation.push_back(new ParametrizationTemporalConstant(config));
      if(readConfigChoiceElement(config, "trend",           type, "Linear trend"))
        representation.push_back(new ParametrizationTemporalTrend(config));
      if(readConfigChoiceElement(config, "splines",         type, "Splines"))
        representation.push_back(new ParametrizationTemporalSplines(config));
      if(readConfigChoiceElement(config, "polynomial",      type, "Polynomial"))
        representation.push_back(new ParametrizationTemporalPolynomial(config));
      if(readConfigChoiceElement(config, "oscillation",     type, "Oscillations"))
        representation.push_back(new ParametrizationTemporalOscillation(config));
      if(readConfigChoiceElement(config, "fourier",         type, "Fourier series"))
        representation.push_back(new ParametrizationTemporalFourier(config));
      if(readConfigChoiceElement(config, "doodsonHarmonic", type, "Tidal variations"))
        representation.push_back(new ParametrizationTemporalDoodsonHarmonic(config));
      endChoice(config);
      if(isCreateSchema(config))
        return;
    }

    parameterCount_ = 0;
    index.resize(representation.size());
    for(UInt i=0; i<representation.size(); i++)
    {
      index.at(i) = parameterCount_;
      parameterCount_ += representation.at(i)->parameterCount();
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

ParametrizationTemporal::~ParametrizationTemporal()
{
  for(UInt i=0; i<representation.size(); i++)
    delete representation.at(i);
}

/***********************************************/

void ParametrizationTemporal::setInterval(const Time &timeStart, const Time &timeEnd, Bool estimatePerArc)
{
  try
  {
    if(timeStart>=timeEnd)
      throw(Exception("wrong time interval: "+timeStart.dateTimeStr()+" - "+timeEnd.dateTimeStr()));

    for(UInt i=0; i<representation.size(); i++)
      representation.at(i)->setInterval(timeStart, timeEnd, estimatePerArc);

    // count parameter
    parameterCount_ = 0;
    index.resize(representation.size());
    for(UInt i=0; i<representation.size(); i++)
    {
      index.at(i) = parameterCount_;
      parameterCount_ += representation.at(i)->parameterCount();
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void ParametrizationTemporal::factors(const Time &time, std::vector<UInt> &index, std::vector<Double> &factor) const
{
  index.clear();
  factor.clear();
  for(UInt i=0; i<representation.size(); i++)
    representation.at(i)->factors(time, this->index.at(i), index, factor);
}

/***********************************************/

Vector ParametrizationTemporal::factors(const Time &time) const
{
  std::vector<UInt>   index;
  std::vector<Double> factor;
  factors(time, index, factor);
  Vector A(parameterCount());
  for(UInt i=0; i<index.size(); i++)
    A(index.at(i)) = factor.at(i);
  return A;
}

/***********************************************/

void ParametrizationTemporal::designMatrix(const Time &time, const const_MatrixSlice &B, MatrixSliceRef A) const
{
  std::vector<UInt>   index;
  std::vector<Double> factor;
  factors(time, index, factor);
  A.slice(0, 0, B.rows(), B.columns()*parameterCount_).setNull();
  for(UInt i=0; i<index.size(); i++)
    axpy(factor.at(i), B, A.slice(0, index.at(i)*B.columns(), B.rows(), B.columns()));
}

/***********************************************/

void ParametrizationTemporal::parameterName(std::vector<ParameterName> &name) const
{
  for(UInt i=0; i<representation.size(); i++)
    representation.at(i)->parameterName(name);
}

/***********************************************/

void ParametrizationTemporal::parameterName(const std::vector<ParameterName> &baseName, std::vector<ParameterName> &name) const
{
  std::vector<ParameterName> temporalName;
  parameterName(temporalName);
  for(UInt i=0; i<temporalName.size(); i++)
    for(UInt k=0; k<baseName.size(); k++)
    {
      ParameterName param = baseName.at(k);
      param.interval = temporalName.at(i).interval;
      param.temporal = temporalName.at(i).temporal;
      name.push_back(param);
    }
}

/***********************************************/
/***********************************************/
