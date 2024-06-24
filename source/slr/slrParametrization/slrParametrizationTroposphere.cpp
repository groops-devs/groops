/***********************************************/
/**
* @file slrParametrizationTroposphere.cpp
*
* @brief Troposphere.
* @see SlrParametrization
*
* @author Torsten Mayer-Guerr
* @date 2022-04-28
*
*/
/***********************************************/

#include "base/import.h"
#include "config/config.h"
#include "files/fileMatrix.h"
#include "classes/troposphere/troposphere.h"
#include "classes/parametrizationTemporal/parametrizationTemporal.h"
#include "classes/platformSelector/platformSelector.h"
#include "slr/slrParametrization/slrParametrizationTroposphere.h"

/***********************************************/

SlrParametrizationTroposphere::SlrParametrizationTroposphere(Config &config)
{
  try
  {
    readConfig(config, "name",                  name,             Config::OPTIONAL, "parameter.troposphere", "used for parameter selection");
    readConfig(config, "selectStations",        selectorStations, Config::MUSTSET,  "",  "");
    readConfig(config, "outputfileTroposphere", fileNameTropo,    Config::OPTIONAL, "output/troposphere_{loopTime:%D}.{station}.txt", "columns: MJD, ZHD, ZWD, dry north gradient, wet north gradient, dry east gradient, wet east gradient");
    readConfig(config, "troposphere",           troposphere,      Config::MUSTSET,  R"({"mendesAndPavlis":{}})",  "a priori troposphere model");
    readConfig(config, "troposphereEstimation", parametrization,  Config::DEFAULT,  "",  "[m] parametrization of zenith delays");
    if(isCreateSchema(config)) return;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

SlrParametrizationTroposphere::~SlrParametrizationTroposphere()
{
  for(Parameter *para : parameters)
    delete para;
}

/***********************************************/

void SlrParametrizationTroposphere::init(Slr *slr, const std::vector<const SlrParametrizationGravityField*> &/*paramGravityField*/)
{
  try
  {
    this->slr = slr;
    auto selectedStations = slr->selectStations(selectorStations);
    parameters.resize(slr->stations.size(), nullptr);
    UInt idTropo = 0;
    std::vector<std::string> names;
    std::vector<Vector3d>    positions;
    for(UInt idStat=0; idStat<slr->stations.size(); idStat++)
      if(selectedStations.at(idStat) && slr->stations.at(idStat)->useable())
      {
        names.push_back(slr->stations.at(idStat)->name());
        positions.push_back(slr->stations.at(idStat)->position(slr->times.front()));
        auto para = new Parameter();
        parameters.at(idStat) = para;
        para->idStat = idStat;
        para->idTropo   = idTropo++;
        para->x         = Vector(parametrization->parameterCount());
        para->zenitDelay.resize(slr->times.size(), 0);
      }

    troposphere->init(names, positions);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void SlrParametrizationTroposphere::initParameter(SlrNormalEquationInfo &normalEquationInfo)
{
  try
  {
    for(auto para : parameters)
      if(para)
        para->index = SlrParameterIndex();
    if(!isEnabled(normalEquationInfo, name))
      return;

    // wet troposphere
    UInt countPara = 0;
    if(parametrization->parameterCount())
      for(auto para : parameters)
        if(para && normalEquationInfo.estimateStation.at(para->idStat))
        {
          std::vector<ParameterName> parameterNames;
          parametrization->parameterName({ParameterName(slr->stations.at(para->idStat)->name(), "troposphere")}, parameterNames);
          para->index = normalEquationInfo.parameterNamesStation(para->idStat, parameterNames);
          countPara += parameterNames.size();
        }
    if(countPara)
      logInfo<<countPara%"%9i troposphere parameters"s<<Log::endl;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void SlrParametrizationTroposphere::observationCorrections(SlrObservationEquation &eqn) const
{
  try
  {
    auto para = parameters.at(eqn.station->idStat());
    if(!para)
      return;

    for(UInt idEpoch=0; idEpoch<eqn.timesTrans.size(); idEpoch++)
    {
      // apriori value
      Double delay = troposphere->slantDelay(para->idTropo, eqn.timesTrans.at(idEpoch), LIGHT_VELOCITY/eqn.laserWavelength(idEpoch),
                                             eqn.azimutStat.at(idEpoch), eqn.elevationStat.at(idEpoch));

      // estimated wet effect
      if(parametrization->parameterCount())
        delay += troposphere->mappingFunctionWet(para->idTropo, eqn.timesTrans.at(idEpoch), LIGHT_VELOCITY/eqn.laserWavelength(idEpoch),
                                                 eqn.azimutStat.at(idEpoch), eqn.elevationStat.at(idEpoch))
               * inner(parametrization->factors(eqn.timesTrans.at(idEpoch)), para->x);

      eqn.l(idEpoch) -= delay;
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void SlrParametrizationTroposphere::aprioriParameter(const SlrNormalEquationInfo &normalEquationInfo, MatrixSliceRef x0) const
{
  try
  {
    // update wet troposphere
    for(auto para : parameters)
      if(para && para->index)
        copy(para->x, x0.row(normalEquationInfo.index(para->index), para->x.rows()));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void SlrParametrizationTroposphere::designMatrix(const SlrNormalEquationInfo &/*normalEquationInfo*/, const SlrObservationEquation &eqn, SlrDesignMatrix &A) const
{
  try
  {
    auto para = parameters.at(eqn.station->idStat());
    if(!para || !para->index)
      return;

    for(UInt idEpoch=0; idEpoch<eqn.timesTrans.size(); idEpoch++)
    {
      const Double mappingFunction = troposphere->mappingFunctionWet(para->idTropo, eqn.timesTrans.at(idEpoch), LIGHT_VELOCITY/eqn.laserWavelength(idEpoch),
                                                                        eqn.azimutStat.at(idEpoch), eqn.elevationStat.at(idEpoch));
      const Double B = mappingFunction * eqn.A(idEpoch, SlrObservationEquation::idxRange);
      std::vector<UInt>   idx;
      std::vector<Double> factor;
      parametrization->factors(eqn.timesTrans.at(idEpoch), idx, factor);
      MatrixSlice Design(A.column(para->index));
      for(UInt i=0; i<factor.size(); i++)
        Design(idEpoch, idx.at(i)) = factor.at(i) * B;
    } // for(idEpoch)
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Double SlrParametrizationTroposphere::updateParameter(const SlrNormalEquationInfo &normalEquationInfo, const_MatrixSliceRef x, const_MatrixSliceRef /*Wz*/)
{
  try
  {
    // update wet troposphere
    Double maxChange = 0;
    Slr::InfoParameterChange info("mm");
    for(auto para : parameters)
      if(para && para->index)
      {
        auto recv = slr->stations.at(para->idStat);
        para->x += x.row(normalEquationInfo.index(para->index), parametrization->parameterCount());
        std::vector<UInt>   index;
        std::vector<Double> factor;
        for(UInt idEpoch=0; idEpoch<slr->times.size(); idEpoch++)
        {
          parametrization->factors(slr->times.at(idEpoch), index, factor);
          Double z = 0;
          for(UInt k=0; k<factor.size(); k++)
            z += factor.at(k) * para->x(index.at(k));
          const Double zOld = para->zenitDelay.at(idEpoch);
          para->zenitDelay.at(idEpoch) = z;
          if(info.update(1e3*(z-zOld)))
            info.info = "troposphere wet ("+recv->name()+", "+slr->times.at(idEpoch).dateTimeStr()+")";
        }
      }
    info.print(1e-3, maxChange);

    return maxChange;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void SlrParametrizationTroposphere::writeResults(const SlrNormalEquationInfo &normalEquationInfo, const std::string &suffix) const
{
  try
  {
    if(!isEnabled(normalEquationInfo, name) || fileNameTropo.empty())
      return;

    VariableList varList;
    varList.setVariable("station", "****");
    for(auto para : parameters)
      if(para && normalEquationInfo.estimateStation.at(para->idStat))
      {
        Double laserWavelength = 0;
        for(const auto &obsSat : slr->stations.at(para->idStat)->observations)
          if(obsSat.size())
          {
            laserWavelength = obsSat.front()->laserWavelength(0);
            break;
          }

        Matrix A(slr->times.size(), 12);
        std::vector<Time> times;
        for(UInt idEpoch=0; idEpoch<slr->times.size(); idEpoch++)
        {
          Double zenithWetDelay, zenithDryDelay, gradientWetNorth, gradientDryNorth, gradientWetEast, gradientDryEast, aDry, aWet;
          troposphere->getAprioriValues(para->idTropo, slr->times.at(idEpoch), LIGHT_VELOCITY/laserWavelength,
                                        zenithDryDelay, zenithWetDelay,
                                        gradientDryNorth, gradientWetNorth, gradientDryEast, gradientWetEast, aDry, aWet);
          times.push_back(slr->times.at(idEpoch));
          A(idEpoch,  0) = slr->times.at(idEpoch).mjd();
          A(idEpoch,  1) = zenithDryDelay;                                // tropospheric zenith dry delay [m] (only from model)
          A(idEpoch,  2) = zenithWetDelay + para->zenitDelay.at(idEpoch); // tropospheric zenith wet delay [m] (model + delta estimate)
          A(idEpoch,  3) = gradientDryNorth;                              // tropospheric dry gradient - north direction [m] (model + delta estimate, due to same mapping function)
          A(idEpoch,  4) = gradientWetNorth;                              // tropospheric wet gradient - north direction [m] (only from model)
          A(idEpoch,  5) = gradientDryEast;                               // tropospheric dry gradient - east component [m] (model + delta estimate, due to same mapping function)
          A(idEpoch,  6) = gradientWetEast;                               // tropospheric wet gradient - east component [m] (only from model)
          A(idEpoch,  7) = para->zenitDelay.at(idEpoch);                  // tropospheric zenith wet delay [m] (delta estimate)
          A(idEpoch,  8) = 0;                                             // tropospheric gradient - north [m] (delta estimate)
          A(idEpoch,  9) = 0;                                             // tropospheric gradient - east  [m] (delta estimate)
          A(idEpoch, 10) = aDry;                                          // dry mapping function coefficient a []
          A(idEpoch, 11) = aWet;                                          // wet mapping function coefficient a []
        }

        varList.setVariable("station", slr->stations.at(para->idStat)->name());
        InstrumentFile::write(fileNameTropo(varList).appendBaseName(suffix), Arc(times, A));
      }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
