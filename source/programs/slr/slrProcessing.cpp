/***********************************************/
/**
* @file slrProcessing.cpp
*
* @brief Satellite Laser Ranging (SLR) processing.
*
* @author Torsten Mayer-Guerr
* @date 2022-04-26
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program processes SLR normal point or full rate observations. It calculates the linearized observation equations,
accumulates them into a system of normal equations and solves it.

To calculate observation equations from the passes, the model parameters or unknown parameters need to be
defined beforehand. These unknown parameters can be chosen arbitrarily by the user with an adequate list of defined
\configClass{parametrization}{slrParametrizationType}.
Some of the \configClass{parametrization}{slrParametrizationType} also include a priori models.

Lastly it is required to define the process flow of the SLR processing. This is accomplished
with a list of \configClass{processingSteps}{slrProcessingStepType}.
Each step is processed consecutively. Some steps allow the selection of parameters, station, or satellites,
which affects all subsequent steps.

The \configClass{timeSeries}{timeSeriesType} is used to precompute Earth rotation and station displacements
with a uniform sampling. In a second step these values are interpolated to the observation epochs.
A sampling of about 10 minutes should be adequate.

It should be noted that GROOPS uses GPS time format, but normal point/full rate data files and CPF files, provided by ILRS data centers
are given in UTC time format.
%See also \program{SlrSimulation}.
)";

/***********************************************/

#include "programs/program.h"
#include "classes/timeSeries/timeSeries.h"
#include "classes/earthRotation/earthRotation.h"
#include "slr/slr.h"
#include "slr/slrStationGenerator/slrStationGenerator.h"
#include "slr/slrSatelliteGenerator/slrSatelliteGenerator.h"
#include "slr/slrParametrization/slrParametrization.h"
#include "slr/slrProcessingStep/slrProcessingStep.h"

/***** CLASS ***********************************/

/** @brief Satellite Laser Ranging (SLR) processing.
* @ingroup programsGroup */
class SlrProcessing
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(SlrProcessing, SINGLEPROCESS, "Satellite Laser Ranging (SLR) processing", Slr)

/***********************************************/

void SlrProcessing::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    TimeSeriesPtr            timeSeries;
    SlrSatelliteGeneratorPtr satelliteGenerator;
    SlrStationGeneratorPtr   stationGenerator;
    SlrParametrizationPtr    slrParametrization;
    EarthRotationPtr         earthRotation;
    SlrProcessingStepPtr     processingSteps;

    readConfig(config, "timeSeries",      timeSeries,          Config::MUSTSET, "",  "defines station movements and earth rotation epochs");
    readConfig(config, "satellite",       satelliteGenerator,  Config::MUSTSET, "",  "satellites");
    readConfig(config, "station",         stationGenerator,    Config::MUSTSET, "",  "ground station network");
    readConfig(config, "earthRotation",   earthRotation,       Config::MUSTSET, "",  "apriori earth rotation");
    readConfig(config, "parametrization", slrParametrization,  Config::MUSTSET, "",  "models and parameters");
    readConfig(config, "processingStep",  processingSteps,     Config::MUSTSET, "",  "steps are processed consecutively");
    if(isCreateSchema(config)) return;

    // ============================

    // init the SLR system
    // --------------------
    logInfo<<"Init SLR"<<Log::endl;
    SlrPtr slr = std::make_shared<Slr>();
    slr->init(timeSeries->times(), satelliteGenerator, stationGenerator, earthRotation, slrParametrization);
//     stationGenerator->preprocessing(slr.get(), comm);
    satelliteGenerator = nullptr;
    stationGenerator   = nullptr;
    slrParametrization = nullptr;
    earthRotation      = nullptr;
    logInfo<<"summary:"<<Log::endl;
    logInfo<<"  satellites:   "<<std::count_if(slr->satellites.begin(), slr->satellites.end(), [](auto s) {return s->useable();})<<Log::endl;
    logInfo<<"  stations:     "<<std::count_if(slr->stations.begin(),   slr->stations.end(),   [](auto s) {return s->useable();})<<Log::endl;
    if(!std::any_of(slr->satellites.begin(), slr->satellites.end(), [](auto s){return s->useable();}))
    {
      logWarningOnce<<"no useable satellites"<<Log::endl;
      return;
    }
    if(!std::any_of(slr->stations.begin(), slr->stations.end(), [](auto s){return s->useable();}))
    {
      logWarningOnce<<"no useable stations"<<Log::endl;
      return;
    }

    // count observations
    // ------------------
    UInt countObs    = 0;
    UInt countPasses = 0;
    for(auto &station : slr->stations)
      for(auto &obsSat : station->observations)
        for(auto &obs : obsSat)
        {
          countObs += obs->observations.rows();
          countPasses++;
        }
    logInfo<<"  observations: "<<countObs<<Log::endl;
    logInfo<<"  passes:       "<<countPasses<<Log::endl;

    // Processing steps
    // ----------------
    SlrProcessingStep::State state(slr);
    processingSteps->process(state);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
