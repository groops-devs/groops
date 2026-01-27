/***********************************************/
/**
* @file slrParametrizationRangeBiasStation.h
*
* @brief Range biases.
* @see SlrParametrization
*
* @author Torsten Mayer-Guerr
* @date 2022-04-28
*
*/
/***********************************************/

#ifndef __GROOPS_SLRPARAMETRIZATIONRANGEBIASSTATION__
#define __GROOPS_SLRPARAMETRIZATIONRANGEBIASSTATION__

// Latex documentation
#ifdef DOCSTRING_SlrParametrization
static const char *docstringSlrParametrizationRangeBiasStation = R"(
\subsection{RangeBiasStation}\label{slrParametrizationType:rangeBiasStation}
Estimates a constant station range bias in $[m]$ for
\configClass{selectStations}{platformSelectorType}.

The \file{parameter names}{parameterName} are \verb|<station>:rangeBias::|.
)";
#endif

/***********************************************/

#include "classes/platformSelector/platformSelector.h"
#include "slr/slrParametrization/slrParametrization.h"

/***** CLASS ***********************************/

/** @brief Range biases.
* @ingroup slrParametrizationGroup
* @see SlrParametrization */
class SlrParametrizationRangeBiasStation : public SlrParametrizationBase
{
  class Parameter
  {
  public:
    SlrStationPtr     station;
    SlrParameterIndex index;
    Double            range;
  };

  Slr                    *slr;
  std::string             name;
  PlatformSelectorPtr     selectorStations;
  std::vector<Parameter*> parameters;
  FileName                fileNameRangeBias;

public:
  SlrParametrizationRangeBiasStation(Config &config);
 ~SlrParametrizationRangeBiasStation();

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

inline SlrParametrizationRangeBiasStation::SlrParametrizationRangeBiasStation(Config &config)
{
  try
  {
    readConfig(config, "name",                name,              Config::OPTIONAL, "parameter.rangeBiasStation", "used for parameter selection");
    readConfig(config, "selectStations",      selectorStations,  Config::MUSTSET,  R"(["all"])", "");
    readConfig(config, "outputfileRangeBias", fileNameRangeBias, Config::MUSTSET,  "rangeBias_{station}.txt", "variable {station} available");
    if(isCreateSchema(config)) return;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline SlrParametrizationRangeBiasStation::~SlrParametrizationRangeBiasStation()
{
  for(auto para : parameters)
    delete para;
}

/***********************************************/

inline void SlrParametrizationRangeBiasStation::init(Slr *slr, const std::vector<const SlrParametrizationGravityField*> &/*paramGravityField*/)
{
  try
  {
    this->slr = slr;

    auto selectedStations = slr->selectStations(selectorStations);
    parameters.resize(slr->stations.size(), nullptr);
    for(UInt idStat=0; idStat<slr->stations.size(); idStat++)
      if(selectedStations.at(idStat))
      {
        auto para = new Parameter();
        parameters.at(idStat) = para;
        para->station = slr->stations.at(idStat);
        para->range   = 0.;
      }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline void SlrParametrizationRangeBiasStation::observationCorrections(SlrObservationEquation &eqn) const
{
  try
  {
    const UInt idStat = eqn.station->idStat();
    if((eqn.type == SlrObservationEquation::RANGE) && (parameters.size() > idStat) && parameters.at(idStat))
      eqn.l -= parameters.at(idStat)->range;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline void SlrParametrizationRangeBiasStation::initParameter(SlrNormalEquationInfo &normalEquationInfo)
{
  try
  {
    for(auto para : parameters)
      if(para)
        para->index = SlrParameterIndex();
    if(!isEnabled(normalEquationInfo, name))
      return;

    UInt countPara = 0;
    for(auto para : parameters)
      if(para && para->station->useable() && normalEquationInfo.estimateStation.at(para->station->idStat()))
      {
        para->index = normalEquationInfo.parameterNamesStation(para->station->idStat(), {ParameterName(para->station->name(), "rangeBias")});
        countPara++;
      }
    if(countPara)
      logInfo<<countPara%"%9i station range bias parameters"s<<Log::endl;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline void SlrParametrizationRangeBiasStation::aprioriParameter(const SlrNormalEquationInfo &normalEquationInfo, MatrixSliceRef x0) const
{
  try
  {
    for(auto para : parameters)
      if(para && para->index)
        x0(normalEquationInfo.index(para->index), 0) = para->range;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline void SlrParametrizationRangeBiasStation::designMatrix(const SlrNormalEquationInfo &/*normalEquationInfo*/, const SlrObservationEquation &eqn, SlrDesignMatrix &A) const
{
  try
  {
    auto para = this->parameters.at(eqn.station->idStat());
    if((eqn.type == SlrObservationEquation::RANGE) && para && para->index)
      axpy(1., eqn.A.column(SlrObservationEquation::idxRange), A.column(para->index));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline Double SlrParametrizationRangeBiasStation::updateParameter(const SlrNormalEquationInfo &normalEquationInfo, const_MatrixSliceRef x, const_MatrixSliceRef /*Wz*/)
{
  try
  {
    Double maxChange = 0;
    Slr::InfoParameterChange infoStat("mm");
    for(auto para : parameters)
      if(para && para->index)
      {
        const Double dBias = x(normalEquationInfo.index(para->index), 0);
        para->range += dBias;
        if(infoStat.update(1e3*dBias))
          infoStat.info = "range bias station ("+para->station->name()+")";
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

inline void SlrParametrizationRangeBiasStation::writeResults(const SlrNormalEquationInfo &normalEquationInfo, const std::string &suffix) const
{
  try
  {
    if(!isEnabled(normalEquationInfo, name))
      return;

    if(!fileNameRangeBias.empty())
    {
      VariableList varList;
      varList.setVariable("station", "****");
      logStatus<<"write range bias to files <"<<fileNameRangeBias(varList).appendBaseName(suffix)<<">"<<Log::endl;
      for(auto para : parameters)
        if(para && para->index)
        {
          MiscValueEpoch epoch;
          epoch.time  = slr->times.front();
          epoch.value = para->range;
          MiscValueArc arc;
          arc.push_back(epoch);
          varList.setVariable("station", para->station->name());
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
