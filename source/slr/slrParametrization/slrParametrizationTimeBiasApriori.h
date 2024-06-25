/***********************************************/
/**
* @file slrParametrizationTimeBiasApriori.h
*
* @brief Time biases.
* @see SlrParametrization
*
* @author Torsten Mayer-Guerr
* @date 2022-04-28
*
*/
/***********************************************/

#ifndef __GROOPS_SLRPARAMETRIZATIONTIMEBIASAPRIORI__
#define __GROOPS_SLRPARAMETRIZATIONTIMEBIASAPRIORI__

// Latex documentation
#ifdef DOCSTRING_SlrParametrization
static const char *docstringSlrParametrizationTimeBiasApriori = R"(
\subsection{TimeBiasApriori}\label{slrParametrizationType:timeBiasApriori}
A priori time bias value for all \configClass{selectStations}{platformSelectorType}.
The \href{https://ilrs.gsfc.nasa.gov/}{ILRS} provides the mean time biases \href{https://ilrs.gsfc.nasa.gov/network/site_information/data_correction.html}{ILRS Data Handling File},
but these have been determined using the passive satellites LAGEOS and Etalon and are therefore only suitable for passive
satellites and not for active ones.
Use \program{SlrSinexDataHandling2Files} to convert the time biases from
\href{https://ilrs.gsfc.nasa.gov/network/site_information/data_correction.html}{ILRS Data Handling File} to \file{instrument file}{instrument}.
)";
#endif

/***********************************************/

#include "classes/parametrizationTemporal/parametrizationTemporal.h"
#include "classes/platformSelector/platformSelector.h"
#include "slr/slrParametrization/slrParametrization.h"

/***** CLASS ***********************************/

/** @brief Range biases.
* @ingroup slrParametrizationGroup
* @see SlrParametrization */
class SlrParametrizationTimeBiasApriori : public SlrParametrizationBase
{
  FileName            fileNameTimeBias;
  PlatformSelectorPtr selectorStations;

public:
  SlrParametrizationTimeBiasApriori(Config &config);
 ~SlrParametrizationTimeBiasApriori() {}

  void init(Slr *slr, const std::vector<const SlrParametrizationGravityField*> &paramGravityField) override;
};

/***********************************************/

inline SlrParametrizationTimeBiasApriori::SlrParametrizationTimeBiasApriori(Config &config)
{
  try
  {
    readConfig(config, "selectStations",    selectorStations, Config::MUSTSET,  R"(["all"])", "");
    readConfig(config, "inputfileTimeBias", fileNameTimeBias, Config::MUSTSET,  "{groopsDataDir}/slr/stations/timeBias/timeBias.{station}.txt", "variable {station} available");
    if(isCreateSchema(config)) return;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline void SlrParametrizationTimeBiasApriori::init(Slr *slr, const std::vector<const SlrParametrizationGravityField*> &/*paramGravityField*/)
{
  try
  {
    Bool found = FALSE;
    auto selectedStations = slr->selectStations(selectorStations);
    VariableList varList;
    for(UInt idStat=0; idStat<slr->stations.size(); idStat++)
      if(selectedStations.at(idStat))
      {
        auto station = slr->stations.at(idStat);
        varList.setVariable("station", station->name());
        try
        {
          MiscValuesArc arc = InstrumentFile::read(fileNameTimeBias(varList));
          Vector bias(station->times.size());
          UInt idArc = 0;
          for(UInt i=0; i<bias.rows(); i++)
          {
            while((idArc+1 < arc.size()) && (arc.at(idArc).time <= station->times.at(i)))
              idArc++;
            if(arc.at(idArc).time > station->times.at(i))
              break;
            bias(i) = arc.at(idArc).values(0) + arc.at(idArc).values(1) * (station->times.at(i)-arc.at(idArc).time).mjd();
          }
          if(quadsum(bias))
          {
            if(station->timeBiases.size())
              station->timeBiases += bias;
            else
              station->timeBiases = bias;
          }
          found = TRUE;
        }
        catch(std::exception &/*e*/)
        {
          // logWarningOnce<<"Unable to read time bias file <"<<fileNameTimeBias(varList)<<">, ignoring."<<Log::endl;
        }
      }

    if(!found)
    {
      varList.setVariable("station", "****");
      logWarningOnce<<"Initialization of all time bias failed. Wrong file name <"<<fileNameTimeBias(varList)<<">?"<<Log::endl;
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

#endif
