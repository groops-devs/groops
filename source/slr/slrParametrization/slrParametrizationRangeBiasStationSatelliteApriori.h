/***********************************************/
/**
* @file slrParametrizationRangeBiasStationSatelliteApriori.h
*
* @brief Range biases.
* @see SlrParametrization
*
* @author Torsten Mayer-Guerr
* @date 2022-04-28
*
*/
/***********************************************/

#ifndef __GROOPS_SLRPARAMETRIZATIONRANGEBIASSTATIONSATELLITEAPRIORI__
#define __GROOPS_SLRPARAMETRIZATIONRANGEBIASSTATIONSATELLITEAPRIORI__

// Latex documentation
#ifdef DOCSTRING_SlrParametrization
static const char *docstringSlrParametrizationRangeBiasStationSatelliteApriori = R"(
\subsection{RangeBiasStationSatelliteApriori}\label{slrParametrizationType:rangeBiasStationSatelliteApriori}
A priori station-satellite range bias value between all \configClass{selectStations}{platformSelectorType} -
\configClass{selectSatellites}{platformSelectorType} pairs.

For standard \href{https://ilrs.gsfc.nasa.gov/}{ILRS} processing this class should be setup twice.
Once for the model from José Rodríguez (see \program{SlrComModel2RangeBiasStationSatellite}) and additionally for
biases from the \href{https://ilrs.gsfc.nasa.gov/network/site_information/data_correction.html}{ILRS Data Handling File} converted with \program{SlrSinexDataHandling2Files}.
)";
#endif

/***********************************************/

#include "classes/platformSelector/platformSelector.h"
#include "slr/slrParametrization/slrParametrization.h"

/***** CLASS ***********************************/

/** @brief Range biases.
* @ingroup slrParametrizationGroup
* @see SlrParametrization */
class SlrParametrizationRangeBiasStationSatelliteApriori : public SlrParametrizationBase
{
  PlatformSelectorPtr selectorStations, selectorSatellites;
  FileName            fileNameRangeBias;
  std::vector<std::vector<MiscValueArc>> range; // per station, satellite

public:
  SlrParametrizationRangeBiasStationSatelliteApriori(Config &config);
 ~SlrParametrizationRangeBiasStationSatelliteApriori() {}

  void init(Slr *slr, const std::vector<const SlrParametrizationGravityField*> &paramGravityField) override;
  void observationCorrections(SlrObservationEquation &eqn) const override;
};

/***********************************************/

inline SlrParametrizationRangeBiasStationSatelliteApriori::SlrParametrizationRangeBiasStationSatelliteApriori(Config &config)
{
  try
  {
    readConfig(config, "selectStations",     selectorStations,   Config::MUSTSET, R"(["all"])", "");
    readConfig(config, "selectSatellites",   selectorSatellites, Config::MUSTSET, R"(["all"])", "");
    readConfig(config, "inputfileRangeBias", fileNameRangeBias,  Config::MUSTSET, "{groopsDataDir}/slr/rangeBiasStationSatellite/modelRodriguez/rangeBias.{station}.{satellite}.txt", "variable {station} and {satellite} available");
    if(isCreateSchema(config)) return;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline void SlrParametrizationRangeBiasStationSatelliteApriori::init(Slr *slr, const std::vector<const SlrParametrizationGravityField*> &/*paramGravityField*/)
{
  try
  {
    Bool found = FALSE;
    auto selectedStations   = slr->selectStations(selectorStations);
    auto selectedSatellites = slr->selectSatellites(selectorSatellites);
    VariableList varList;
    for(UInt idStat=0; idStat<slr->stations.size(); idStat++)
      if(selectedStations.at(idStat))
        for(UInt idSat=0; idSat<slr->satellites.size(); idSat++)
          if(selectedSatellites.at(idSat))
          {
            varList.setVariable("station",   slr->stations.at(idStat)->name());
            varList.setVariable("satellite", slr->satellites.at(idSat)->name());
            try
            {
              MiscValueArc arc = InstrumentFile::read(fileNameRangeBias(varList));
              if(range.size() <= idStat)
                range.resize(idStat+1);
              if(range.at(idStat).size() <= idSat)
                range.at(idStat).resize(idSat+1);
              range.at(idStat).at(idSat) = arc;
              found = TRUE;
            }
            catch(std::exception &/*e*/)
            {
              // logWarningOnce<<"Unable to read range bias file <"<<fileNameRangeBias(varList)<<">, ignoring."<<Log::endl;
            }
          }

    if(!found)
    {
      varList.setVariable("station",   "****");
      varList.setVariable("satellite", "****");
      logWarningOnce<<"Initialization of all range bias failed. Wrong file name <"<<fileNameRangeBias(varList)<<">?"<<Log::endl;
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline void SlrParametrizationRangeBiasStationSatelliteApriori::observationCorrections(SlrObservationEquation &eqn) const
{
  try
  {
    const UInt idStat = eqn.station->idStat();
    const UInt idSat  = eqn.satellite->idSat();
    if((eqn.type != SlrObservationEquation::RANGE) || (range.size() <= idStat) || (range.at(idStat).size() <= idSat) || !range.at(idStat).at(idSat).size())
      return;

    const MiscValueArc &arc = range.at(idStat).at(idSat);
    for(UInt i=0; (i<arc.size()) && (arc.at(i).time <= eqn.timesStat.back()); i++)
      if((i+1 >= arc.size()) || (arc.at(i+1).time >= eqn.timesStat.front()))
      {
        eqn.l -= arc.at(i).value;
        break;
      }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

#endif
