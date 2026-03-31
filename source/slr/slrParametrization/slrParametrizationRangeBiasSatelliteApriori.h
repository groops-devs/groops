/***********************************************/
/**
* @file slrParametrizationRangeBiasSatelliteApriori.h
*
* @brief Range biases.
* @see SlrParametrization
*
* @author Torsten Mayer-Guerr
* @date 2022-04-28
*
*/
/***********************************************/

#ifndef __GROOPS_SLRPARAMETRIZATIONRANGEBIASSATELLITEAPRIORI__
#define __GROOPS_SLRPARAMETRIZATIONRANGEBIASSATELLITEAPRIORI__

// Latex documentation
#ifdef DOCSTRING_SlrParametrization
static const char *docstringSlrParametrizationRangeBiasSatelliteApriori = R"(
\subsection{RangeBiasSatelliteApriori}\label{slrParametrizationType:rangeBiasSatelliteApriori}
A priori satellite range bias value for \configClass{selectSatellites}{platformSelectorType}.
)";
#endif

/***********************************************/

#include "classes/platformSelector/platformSelector.h"
#include "slr/slrParametrization/slrParametrization.h"

/***** CLASS ***********************************/

/** @brief Range biases.
* @ingroup slrParametrizationGroup
* @see SlrParametrization */
class SlrParametrizationRangeBiasSatelliteApriori : public SlrParametrizationBase
{
  PlatformSelectorPtr       selectorSatellites;
  FileName                  fileNameRangeBias;
  std::vector<MiscValueArc> range; // per satellite

public:
  SlrParametrizationRangeBiasSatelliteApriori(Config &config);
 ~SlrParametrizationRangeBiasSatelliteApriori() {}

  void init(Slr *slr, const std::vector<const SlrParametrizationGravityField*> &paramGravityField) override;
  void observationCorrections(SlrObservationEquation &eqn) const override;
};

/***********************************************/

inline SlrParametrizationRangeBiasSatelliteApriori::SlrParametrizationRangeBiasSatelliteApriori(Config &config)
{
  try
  {
    readConfig(config, "selectSatellites",   selectorSatellites, Config::MUSTSET,  R"(["all"])", "");
    readConfig(config, "inputfileRangeBias", fileNameRangeBias,  Config::MUSTSET,  "", "variable {satellite} available");
    if(isCreateSchema(config)) return;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline void SlrParametrizationRangeBiasSatelliteApriori::init(Slr *slr, const std::vector<const SlrParametrizationGravityField*> &/*paramGravityField*/)
{
  try
  {
    Bool found = FALSE;
    auto selectedSatellites = slr->selectSatellites(selectorSatellites);
    VariableList varList;
    for(UInt idSat=0; idSat<slr->satellites.size(); idSat++)
      if(selectedSatellites.at(idSat))
      {
        varList.setVariable("satellite", slr->satellites.at(idSat)->name());
        try
        {
          MiscValueArc arc = InstrumentFile::read(fileNameRangeBias(varList));
          if(range.size() <= idSat)
            range.resize(idSat+1);
          range.at(idSat) = arc;
          found = TRUE;
        }
        catch(std::exception &/*e*/)
        {
          // logWarningOnce<<"Unable to read range bias file <"<<fileNameRangeBias(varList)<<">, ignoring."<<Log::endl;
        }
      }

    if(!found)
    {
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

inline void SlrParametrizationRangeBiasSatelliteApriori::observationCorrections(SlrObservationEquation &eqn) const
{
  try
  {
    const UInt idSat = eqn.satellite->idSat();
    if((eqn.type != SlrObservationEquation::RANGE) || (range.size() <= idSat) ||  !range.at(idSat).size())
      return;

    const MiscValueArc &arc = range.at(idSat);
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
