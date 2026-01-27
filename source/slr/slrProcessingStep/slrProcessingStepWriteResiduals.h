/***********************************************/
/**
* @file slrProcessingStepWriteResiduals.h
*
* @brief SLR processing step: WriteResiduals.
*
* @author Torsten Mayer-Guerr
* @date 2022-04-28
*
*/
/***********************************************/

#ifndef __GROOPS_SLRPROCESSINGSTEPWRITERESIDUALS__
#define __GROOPS_SLRPROCESSINGSTEPWRITERESIDUALS__

// Latex documentation
#ifdef DOCSTRING_SlrProcessingStep
static const char *docstringSlrProcessingStepWriteResiduals = R"(
\subsection{WriteResiduals}\label{slrProcessingStepType:writeResiduals}
Writes the \file{observation residuals}{instrument} for all
\configClass{selectStations}{platformSelectorType}. For for each station-satellite
pair a file is written. The file name is interpreted as a template with
the variables \verb|{station}| and \verb|{satellite}| being replaced by the station name.
)";
#endif

/***********************************************/

#include "config/config.h"
#include "classes/platformSelector/platformSelector.h"
#include "slr/slrProcessingStep/slrProcessingStep.h"

/***** CLASS ***********************************/

/** @brief SLR processing step: WriteResiduals.
* @ingroup slrProcessingStepGroup
* @see SlrProcessingStep */
class SlrProcessingStepWriteResiduals : public SlrProcessingStepBase
{
  PlatformSelectorPtr selectorStations, selectorSatellites;
  FileName            fileNameResiduals;

public:
  SlrProcessingStepWriteResiduals(Config &config);
  void process(SlrProcessingStep::State &state) override;
};

/***********************************************/

inline SlrProcessingStepWriteResiduals::SlrProcessingStepWriteResiduals(Config &config)
{
  try
  {
    readConfig(config, "selectStations",      selectorStations,   Config::MUSTSET, "", "subset of used stations");
    readConfig(config, "selectSatellites",    selectorSatellites, Config::MUSTSET, "", "subset of used satellites");
    readConfig(config, "outputfileResiduals", fileNameResiduals,  Config::MUSTSET, "output/residuals_{loopTime:%D}.{station}.{satellite}.dat", "variable {station} available");
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline void SlrProcessingStepWriteResiduals::process(SlrProcessingStep::State &state)
{
  try
  {
    auto selectedStations   = state.slr->selectStations(selectorStations);
    auto selectedSatellites = state.slr->selectSatellites(selectorSatellites);

    VariableList varList;
    varList.setVariable("station",   "****");
    varList.setVariable("satellite", "****");
    logStatus<<"write residuals to file <"<<fileNameResiduals(varList)<<">"<<Log::endl;
    for(auto stat : state.slr->stations)
      if(selectedStations.at(stat->idStat()) && state.normalEquationInfo.estimateStation.at(stat->idStat()) && stat->useable())
        for(auto sat : state.slr->satellites)
          if(selectedSatellites.at(sat->idSat()) && state.normalEquationInfo.estimateSatellite.at(sat->idSat()) && sat->useable())
          {
            std::vector<SatelliteLaserRangingArc> arcList;
            for(auto &obs : stat->observations.at(sat->idSat()))
            {
              SlrObservationEquation eqn;
              eqn.compute(*obs, *stat, *sat, state.slr->funcRotationCrf2Trf, state.slr->funcReduceModels, FALSE);
              Matrix data(obs->observations.rows(), 1+Epoch::dataCount(Epoch::SATELLITELASERRANGING));
              copy(obs->residuals,    data.column(1)); // data1 range [seconds]
              copy(obs->redundancies, data.column(3)); // data3 redundancy
              for(UInt i=0; i<data.rows(); i++)
              {
                data(i, 2) = obs->sigmas(i)/obs->sigmas0(i); // data2 rms
                data(i, 5) = obs->laserWavelength(i);        // data5 laser wavelength
                data(i, 6) = eqn.azimutStat.at(i);           // data5 azmiuth
                data(i, 7) = eqn.elevationStat.at(i);        // data6 elevation
              }

              arcList.push_back(Arc(obs->timesTrans, data, Epoch::SATELLITELASERRANGING));
            }

            if(arcList.size())
            {
              VariableList varList;
              varList.setVariable("station",   stat->name());
              varList.setVariable("satellite", sat->name());
              InstrumentFile::write(fileNameResiduals(varList), arcList);
            }
          } // for(stat, sat)
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

#endif
