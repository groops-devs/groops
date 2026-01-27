/***********************************************/
/**
* @file slrParametrizationTimeBias.cpp
*
* @brief Time biases.
* @see SlrParametrization
*
* @author Torsten Mayer-Guerr
* @date 2022-04-28
*
*/
/***********************************************/

#include "base/import.h"
#include "config/config.h"
#include "classes/parametrizationTemporal/parametrizationTemporal.h"
#include "classes/platformSelector/platformSelector.h"
#include "slr/slr.h"
#include "slr/slrParametrization/slrParametrizationTimeBias.h"

/***********************************************/

SlrParametrizationTimeBias::SlrParametrizationTimeBias(Config &config)
{
  try
  {
    readConfig(config, "name",             name,                    Config::OPTIONAL, "parameter.timeBias", "used for parameter selection");
    readConfig(config, "selectStations",   selectorStations,        Config::MUSTSET,  R"(["all"])", "");
    readConfig(config, "estimateTimeBias", parametrizationTemporal, Config::MUSTSET,  "", "[ms]");
    if(isCreateSchema(config)) return;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

SlrParametrizationTimeBias::~SlrParametrizationTimeBias()
{
  for(auto para : paraStations)
    delete para;
}

/***********************************************/

void SlrParametrizationTimeBias::init(Slr *slr, const std::vector<const SlrParametrizationGravityField*> &/*paramGravityField*/)
{
  try
  {
    this->slr = slr;

    auto selectedStations = slr->selectStations(selectorStations);
    paraStations.resize(slr->stations.size(), nullptr);
    for(UInt idStat=0; idStat<slr->stations.size(); idStat++)
      if(selectedStations.at(idStat) && slr->stations.at(idStat)->useable())
      {
        auto para = new Parameter();
        paraStations.at(idStat) = para;
        para->x       = Vector(parametrizationTemporal->parameterCount());
        para->station = slr->stations.at(idStat);
        if(!para->station->timeBiases.size())
          para->station->timeBiases = Vector(para->station->times.size());
      }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void SlrParametrizationTimeBias::initParameter(SlrNormalEquationInfo &normalEquationInfo)
{
  try
  {
    for(auto para : paraStations)
      if(para)
        para->index = SlrParameterIndex();
    if(!isEnabled(normalEquationInfo, name))
      return;

    UInt countParaStat = 0;
    for(auto para : paraStations)
      if(para && para->station->useable() && normalEquationInfo.estimateStation.at(para->station->idStat()))
      {
        std::vector<ParameterName> parameterNames;
        parametrizationTemporal->parameterName({ParameterName(para->station->name(), "timeBias")}, parameterNames);
        para->index = normalEquationInfo.parameterNamesStation(para->station->idStat(), parameterNames);
        countParaStat += para->x.rows();
      }
    if(countParaStat)
      logInfo<<countParaStat%"%9i time bias parameters"s<<Log::endl;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void SlrParametrizationTimeBias::aprioriParameter(const SlrNormalEquationInfo &normalEquationInfo, MatrixSliceRef x0) const
{
  try
  {
    for(auto para : paraStations)
      if(para && para->index)
        copy(para->x, x0.row(normalEquationInfo.index(para->index), para->x.rows()));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void SlrParametrizationTimeBias::designMatrix(const SlrNormalEquationInfo &/*normalEquationInfo*/, const SlrObservationEquation &eqn, SlrDesignMatrix &A) const
{
  try
  {
    auto paraStations = this->paraStations.at(eqn.station->idStat());
    if(paraStations && paraStations->index)
    {
      MatrixSlice Design(A.column(paraStations->index));
      for(UInt idEpoch=0; idEpoch<eqn.timesStat.size(); idEpoch++)
      {
        const const_MatrixSlice B(eqn.A.slice(eqn.index.at(idEpoch), SlrObservationEquation::idxTime, eqn.count.at(idEpoch), 1));
        std::vector<UInt>   idx;
        std::vector<Double> factor;
        parametrizationTemporal->factors(eqn.timesStat.at(idEpoch), idx, factor);
        for(UInt i=0; i<factor.size(); i++)
          axpy(1e-3*factor.at(i), B, Design.slice(eqn.index.at(idEpoch), idx.at(i), eqn.count.at(idEpoch), 1)); // [ms]
      }
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Double SlrParametrizationTimeBias::updateParameter(const SlrNormalEquationInfo &normalEquationInfo, const_MatrixSliceRef x, const_MatrixSliceRef /*Wz*/)
{
  try
  {
    Double maxChange = 0;
    Slr::InfoParameterChange info("ms");
    for(auto para : paraStations)
      if(para && para->index)
      {
        const Vector dx = 1e-3 * x.row(normalEquationInfo.index(para->index), para->x.rows());
        para->x += dx;

        for(UInt idEpoch=0; idEpoch<para->station->times.size(); idEpoch++)
        {
          Double dt = inner(parametrizationTemporal->factors(para->station->times.at(idEpoch)), dx);
          para->station->timeBiases(idEpoch) += dt;
          if(info.update(1e3*dt))
            info.info = "time bias station ("+para->station->name()+")";
        }
      }
    info.print(0., maxChange);

    return maxChange;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
