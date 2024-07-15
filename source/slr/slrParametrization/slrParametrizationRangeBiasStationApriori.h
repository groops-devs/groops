/***********************************************/
/**
* @file slrParametrizationRangeBiasStationApriori.h
*
* @brief Range biases.
* @see SlrParametrization
*
* @author Torsten Mayer-Guerr
* @date 2022-04-28
*
*/
/***********************************************/

#ifndef __GROOPS_SLRPARAMETRIZATIONRANGEBIASSTATIONAPRIORI__
#define __GROOPS_SLRPARAMETRIZATIONRANGEBIASSTATIONAPRIORI__

// Latex documentation
#ifdef DOCSTRING_SlrParametrization
static const char *docstringSlrParametrizationRangeBiasStationApriori = R"(
\subsection{RangeBiasStationApriori}\label{slrParametrizationType:rangeBiasStationApriori}
A priori station range bias value for all \configClass{selectStations}{platformSelectorType}.
The \href{https://ilrs.gsfc.nasa.gov/}{ILRS} provides the mean range biases \href{https://ilrs.gsfc.nasa.gov/network/site_information/data_correction.html}{ILRS Data Handling File},
but these have been determined using the passive satellites LAGEOS and Etalon and are therefore only suitable for passive
satellites and not for active ones.
Use \program{SlrSinexDataHandling2Files} to convert the range biases from
\href{https://ilrs.gsfc.nasa.gov/network/site_information/data_correction.html}{ILRS Data Handling File} to \file{instrument file}{instrument}.
)";
#endif

/***********************************************/

#include "classes/platformSelector/platformSelector.h"
#include "slr/slrParametrization/slrParametrization.h"

/***** CLASS ***********************************/

/** @brief Range biases.
* @ingroup slrParametrizationGroup
* @see SlrParametrization */
class SlrParametrizationRangeBiasStationApriori : public SlrParametrizationBase
{
  PlatformSelectorPtr       selectorStations;
  FileName                  fileNameRangeBias;
  std::vector<MiscValueArc> range; // per station

public:
  SlrParametrizationRangeBiasStationApriori(Config &config);
 ~SlrParametrizationRangeBiasStationApriori() {}

  void init(Slr *slr, const std::vector<const SlrParametrizationGravityField*> &paramGravityField) override;
  void observationCorrections(SlrObservationEquation &eqn) const override;
};

/***********************************************/

inline SlrParametrizationRangeBiasStationApriori::SlrParametrizationRangeBiasStationApriori(Config &config)
{
  try
  {
    readConfig(config, "selectStations",     selectorStations,  Config::MUSTSET, R"(["all"])", "");
    readConfig(config, "inputfileRangeBias", fileNameRangeBias, Config::MUSTSET, "{groopsDataDir}/slr/stations/rangeBias/rangeBias.{station}.txt", "variable {station} available");
    if(isCreateSchema(config)) return;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline void SlrParametrizationRangeBiasStationApriori::init(Slr *slr, const std::vector<const SlrParametrizationGravityField*> &/*paramGravityField*/)
{
  try
  {
    Bool found = FALSE;
    auto selectedStations   = slr->selectStations(selectorStations);
    VariableList varList;
    for(UInt idStat=0; idStat<slr->stations.size(); idStat++)
      if(selectedStations.at(idStat))
      {
        varList.setVariable("station", slr->stations.at(idStat)->name());
        try
        {
          MiscValueArc arc = InstrumentFile::read(fileNameRangeBias(varList));
          if(range.size() <= idStat)
            range.resize(idStat+1);
          range.at(idStat) = arc;
          found = TRUE;
        }
        catch(std::exception &/*e*/)
        {
          // logWarningOnce<<"Unable to read range bias file <"<<fileNameRangeBias(varList)<<">, ignoring."<<Log::endl;
        }
      }

    if(!found)
    {
      varList.setVariable("station", "****");
      logWarningOnce<<"Initialization of all range bias failed. Wrong file name <"<<fileNameRangeBias(varList)<<">?"<<Log::endl;
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline void SlrParametrizationRangeBiasStationApriori::observationCorrections(SlrObservationEquation &eqn) const
{
  try
  {
    const UInt idStat = eqn.station->idStat();
    if((range.size() <= idStat) ||  !range.at(idStat).size())
      return;

    const MiscValueArc &arc = range.at(idStat);
    for(UInt i=0; (i<arc.size()) && (arc.at(i).time <= eqn.timesTrans.back()); i++)
      if((i+1 >= arc.size()) || (arc.at(i+1).time >= eqn.timesTrans.front()))
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
