/***********************************************/
/**
* @file slrParametrizationRangeBiasStationSatellite.h
*
* @brief Range biases.
* @see SlrParametrization
*
* @author Torsten Mayer-Guerr
* @date 2022-04-28
*
*/
/***********************************************/

#ifndef __GROOPS_SLRPARAMETRIZATIONRANGEBIASSTATIONSATELLITE__
#define __GROOPS_SLRPARAMETRIZATIONRANGEBIASSTATIONSATELLITE__

// Latex documentation
#ifdef DOCSTRING_SlrParametrization
static const char *docstringSlrParametrizationRangeBiasStationSatellite = R"(
\subsection{RangeBiasStationSatellite}\label{slrParametrizationType:rangeBiasStationSatellite}
Estimates the station-satellite range bias in $[m]$ between all
\configClass{selectStations}{platformSelectorType} -
\configClass{selectSatellites}{platformSelectorType} pairs.

The \file{parameter names}{parameterName} are \verb|<station>.<satellite>:rangeBias::|.
)";
#endif

/***********************************************/

#include "classes/platformSelector/platformSelector.h"
#include "slr/slrParametrization/slrParametrization.h"

/***** CLASS ***********************************/

/** @brief Range biases.
* @ingroup slrParametrizationGroup
* @see SlrParametrization */
class SlrParametrizationRangeBiasStationSatellite : public SlrParametrizationBase
{
  class Parameter
  {
  public:
    SlrStationPtr     station;
    SlrSatellitePtr   satellite;
    SlrParameterIndex index;
    Double            range;
  };

  Slr                 *slr;
  std::string          name;
  PlatformSelectorPtr  selectorStations, selectorSatellites;
  std::vector<std::vector<Parameter*>> parameters; // for each station, satellite
  FileName             fileNameRangeBias;

public:
  SlrParametrizationRangeBiasStationSatellite(Config &config);
 ~SlrParametrizationRangeBiasStationSatellite();

  void   init(Slr *slr, const std::vector<const SlrParametrizationGravityField*> &paramGravityField) override;
  void   observationCorrections(SlrObservationEquation &eqn) const override;
  void   initParameter(SlrNormalEquationInfo &normalEquationInfo) override;
  void   aprioriParameter(const SlrNormalEquationInfo &normalEquationInfo, MatrixSliceRef x0) const override;
  void   designMatrix(const SlrNormalEquationInfo &normalEquationInfo, const SlrObservationEquation &eqn, SlrDesignMatrix &A) const override;
  Double updateParameter(const SlrNormalEquationInfo &normalEquationInfo, const_MatrixSliceRef x, const_MatrixSliceRef Wz) override;
  void   writeResults(const SlrNormalEquationInfo &normalEquationInfo, const std::string &suffix) const override;
};

/***********************************************/
/***********************************************/

inline SlrParametrizationRangeBiasStationSatellite::SlrParametrizationRangeBiasStationSatellite(Config &config)
{
  try
  {
    readConfig(config, "name",                name,               Config::OPTIONAL, "parameter.rangeBiasSatellite", "used for parameter selection");
    readConfig(config, "selectStations",      selectorStations,   Config::MUSTSET,  R"(["all"])", "");
    readConfig(config, "selectSatellites",    selectorSatellites, Config::MUSTSET,  R"(["all"])", "");
    readConfig(config, "outputfileRangeBias", fileNameRangeBias,  Config::MUSTSET,  "rangeBias_{station}_{satellite}.txt", "variable {station} and {satellite} available");
    if(isCreateSchema(config)) return;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline SlrParametrizationRangeBiasStationSatellite::~SlrParametrizationRangeBiasStationSatellite()
{
  for(auto &parametersStation : parameters)
    for(auto para : parametersStation)
      delete para;
}

/***********************************************/

inline void SlrParametrizationRangeBiasStationSatellite::init(Slr *slr, const std::vector<const SlrParametrizationGravityField*> &/*paramGravityField*/)
{
  try
  {
    this->slr = slr;

    auto selectedStations   = slr->selectStations(selectorStations);
    auto selectedSatellites = slr->selectSatellites(selectorSatellites);
    parameters.resize(slr->stations.size());
    for(UInt idStat=0; idStat<slr->stations.size(); idStat++)
      if(selectedStations.at(idStat))
        for(UInt idSat=0; idSat<slr->satellites.size(); idSat++)
          if(selectedSatellites.at(idSat))
          {
            auto para = new Parameter();
            parameters.at(idStat).resize(slr->satellites.size(), nullptr);
            parameters.at(idStat).at(idSat) = para;
            para->station   = slr->stations.at(idStat);
            para->satellite = slr->satellites.at(idSat);
            para->range     = 0.;
          }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline void SlrParametrizationRangeBiasStationSatellite::observationCorrections(SlrObservationEquation &eqn) const
{
  try
  {
    const UInt idStat = eqn.station->idStat();
    const UInt idSat  = eqn.satellite->idSat();
    if((eqn.type == SlrObservationEquation::RANGE) && (parameters.size() > idStat) && (parameters.at(idStat).size() > idSat) && parameters.at(idStat).at(idSat))
      eqn.l -= parameters.at(idStat).at(idSat)->range;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline void SlrParametrizationRangeBiasStationSatellite::initParameter(SlrNormalEquationInfo &normalEquationInfo)
{
  try
  {
    for(auto parametersStation : parameters)
      for(auto para : parametersStation)
        if(para)
          para->index = SlrParameterIndex();
    if(!isEnabled(normalEquationInfo, name))
      return;

    UInt countPara = 0;
    for(auto &parametersStation : parameters)
      for(auto para : parametersStation)
        if(para && para->satellite->useable() && normalEquationInfo.estimateSatellite.at(para->satellite->idSat())
                && para->station->useable()   && normalEquationInfo.estimateStation.at(para->station->idStat()))
        {
          para->index = normalEquationInfo.parameterNamesStation(para->station->idStat(),
                                                                 {ParameterName(para->station->name()+"."+para->satellite->name(), "rangeBias")});
          countPara++;
        }
    if(countPara)
      logInfo<<countPara%"%9i satellite range bias parameters"s<<Log::endl;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline void SlrParametrizationRangeBiasStationSatellite::aprioriParameter(const SlrNormalEquationInfo &normalEquationInfo, MatrixSliceRef x0) const
{
  try
  {
    for(auto &parametersStation : parameters)
      for(auto para : parametersStation)
        x0(normalEquationInfo.index(para->index), 0) = para->range;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline void SlrParametrizationRangeBiasStationSatellite::designMatrix(const SlrNormalEquationInfo &/*normalEquationInfo*/, const SlrObservationEquation &eqn, SlrDesignMatrix &A) const
{
  try
  {
    auto &parametersStation = this->parameters.at(eqn.station->idStat());
    if(parametersStation.size())
    {
      auto para = parametersStation.at(eqn.satellite->idSat());
      if((eqn.type == SlrObservationEquation::RANGE) && para && para->index)
        axpy(1., eqn.A.column(SlrObservationEquation::idxRange), A.column(para->index));
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline Double SlrParametrizationRangeBiasStationSatellite::updateParameter(const SlrNormalEquationInfo &normalEquationInfo, const_MatrixSliceRef x, const_MatrixSliceRef /*Wz*/)
{
  try
  {
    Double maxChange = 0;
    Slr::InfoParameterChange infoStat("mm");
    for(auto &parametersStation : parameters)
      for(auto para : parametersStation)
        if(para && para->index)
        {
          const Double dBias = x(normalEquationInfo.index(para->index), 0);
          para->range += dBias;
          if(infoStat.update(1e3*dBias))
            infoStat.info = "range bias station - satellite ("+para->station->name()+"."+para->satellite->name()+")";
        }
    infoStat.print(1e-3, maxChange);

    return maxChange;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline void SlrParametrizationRangeBiasStationSatellite::writeResults(const SlrNormalEquationInfo &normalEquationInfo, const std::string &suffix) const
{
  try
  {
    if(!isEnabled(normalEquationInfo, name))
      return;

    if(!fileNameRangeBias.empty())
    {
      VariableList varList;
      varList.setVariable("station",   "****");
      varList.setVariable("satellite", "****");
      logStatus<<"write range bias to files <"<<fileNameRangeBias(varList).appendBaseName(suffix)<<">"<<Log::endl;
      for(auto &parametersStation : parameters)
        for(auto para : parametersStation)
          if(para && para->index)
          {
            MiscValueEpoch epoch;
            epoch.time  = slr->times.front();
            epoch.value = para->range;
            MiscValueArc arc;
            arc.push_back(epoch);
            varList.setVariable("station",   para->station->name());
            varList.setVariable("satellite", para->satellite->name());
            InstrumentFile::write(fileNameRangeBias(varList).appendBaseName(suffix), arc);
          }
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

#endif
