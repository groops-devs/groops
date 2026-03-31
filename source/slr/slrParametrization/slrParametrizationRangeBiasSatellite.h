/***********************************************/
/**
* @file slrParametrizationRangeBiasSatellite.h
*
* @brief Range biases.
* @see SlrParametrization
*
* @author Torsten Mayer-Guerr
* @date 2022-04-28
*
*/
/***********************************************/

#ifndef __GROOPS_SLRPARAMETRIZATIONRANGEBIASSATELLITE__
#define __GROOPS_SLRPARAMETRIZATIONRANGEBIASSATELLITE__

// Latex documentation
#ifdef DOCSTRING_SlrParametrization
static const char *docstringSlrParametrizationRangeBiasSatellite = R"(
\subsection{RangeBiasSatellite}\label{slrParametrizationType:rangeBiasSatellite}
Estimates a constant satellite range bias in $[m]$ for
\configClass{selectSatellites}{platformSelectorType}.

The \file{parameter names}{parameterName} a \verb|<satellite>:rangeBias::|.
)";
#endif

/***********************************************/

#include "classes/platformSelector/platformSelector.h"
#include "slr/slrParametrization/slrParametrization.h"

/***** CLASS ***********************************/

/** @brief Range biases.
* @ingroup slrParametrizationGroup
* @see SlrParametrization */
class SlrParametrizationRangeBiasSatellite : public SlrParametrizationBase
{
  class Parameter
  {
  public:
    SlrSatellitePtr   satellite;
    SlrParameterIndex index;
    Double            range;
  };

  Slr                    *slr;
  std::string             name;
  PlatformSelectorPtr     selectorSatellites;
  std::vector<Parameter*> parameters;
  FileName                fileNameRangeBias;

public:
  SlrParametrizationRangeBiasSatellite(Config &config);
 ~SlrParametrizationRangeBiasSatellite();

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

inline SlrParametrizationRangeBiasSatellite::SlrParametrizationRangeBiasSatellite(Config &config)
{
  try
  {
    readConfig(config, "name",                     name,                     Config::OPTIONAL, "parameter.rangeBiasSatellite", "used for parameter selection");
    readConfig(config, "selectSatellites",         selectorSatellites,       Config::MUSTSET,  R"(["all"])", "");
    readConfig(config, "outputfileRangeBias",      fileNameRangeBias,        Config::MUSTSET,  "rangeBias_{satellite}.dat", "variable {satellite} available");
    if(isCreateSchema(config)) return;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline SlrParametrizationRangeBiasSatellite::~SlrParametrizationRangeBiasSatellite()
{
  for(auto para : parameters)
    delete para;
}

/***********************************************/

inline void SlrParametrizationRangeBiasSatellite::init(Slr *slr, const std::vector<const SlrParametrizationGravityField*> &/*paramGravityField*/)
{
  try
  {
    this->slr = slr;
    std::vector<Byte> selectedSatellites = slr->selectSatellites(selectorSatellites);
    parameters.resize(slr->satellites.size(), nullptr);
    for(UInt idSat=0; idSat<slr->satellites.size(); idSat++)
      if(selectedSatellites.at(idSat))
      {
        auto para = new Parameter();
        parameters.at(idSat) = para;
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

inline void SlrParametrizationRangeBiasSatellite::observationCorrections(SlrObservationEquation &eqn) const
{
  try
  {
    const UInt idSat = eqn.satellite->idSat();
    if((eqn.type == SlrObservationEquation::RANGE) && (parameters.size() > idSat) && parameters.at(idSat))
      eqn.l -= parameters.at(idSat)->range;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline void SlrParametrizationRangeBiasSatellite::initParameter(SlrNormalEquationInfo &normalEquationInfo)
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
      if(para && para->satellite->useable() && normalEquationInfo.estimateSatellite.at(para->satellite->idSat()))
      {
        para->index = normalEquationInfo.parameterNamesSatellite(para->satellite->idSat(), {ParameterName(para->satellite->name(), "rangeBias")});
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

inline void SlrParametrizationRangeBiasSatellite::aprioriParameter(const SlrNormalEquationInfo &normalEquationInfo, MatrixSliceRef x0) const
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

inline void SlrParametrizationRangeBiasSatellite::designMatrix(const SlrNormalEquationInfo &/*normalEquationInfo*/, const SlrObservationEquation &eqn, SlrDesignMatrix &A) const
{
  try
  {
    auto para = this->parameters.at(eqn.satellite->idSat());
    if((eqn.type == SlrObservationEquation::RANGE) && para && para->index)
      axpy(1., eqn.A.column(SlrObservationEquation::idxRange), A.column(para->index));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline Double SlrParametrizationRangeBiasSatellite::updateParameter(const SlrNormalEquationInfo &normalEquationInfo, const_MatrixSliceRef x, const_MatrixSliceRef /*Wz*/)
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
          infoStat.info = "range bias satellite ("+para->satellite->name()+")";
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

inline void SlrParametrizationRangeBiasSatellite::writeResults(const SlrNormalEquationInfo &normalEquationInfo, const std::string &suffix) const
{
  try
  {
    if(!isEnabled(normalEquationInfo, name))
      return;

    if(!fileNameRangeBias.empty())
    {
      VariableList varList;
      varList.setVariable("satellite", "****");
      logStatus<<"write range bias to files <"<<fileNameRangeBias(varList).appendBaseName(suffix)<<">"<<Log::endl;
      for(auto para : parameters)
        if(para && para->index)
        {
          MiscValueEpoch epoch;
          epoch.time  = slr->times.front();
          epoch.value = para->range;
          MiscValueArc arc;
          arc.push_back(epoch);
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
