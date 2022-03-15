/***********************************************/
/**
* @file gnssProcessingStepDisableTransmitterShadowEpochs.h
*
* @brief GNSS processing step: DisableTransmitterShadowEpochs.
*
* @author Torsten Mayer-Guerr
* @date 2021-09-05
*
*/
/***********************************************/

#ifndef __GROOPS_GNSSPROCESSINGSTEPDISABLETRANSMITTERSHADOWEPOCHS__
#define __GROOPS_GNSSPROCESSINGSTEPDISABLETRANSMITTERSHADOWEPOCHS__

// Latex documentation
#ifdef DOCSTRING_GnssProcessingStep
static const char *docstringGnssProcessingStepDisableTransmitterShadowEpochs = R"(
\subsection{DisableTransmitterShadowEpochs}\label{gnssProcessingStepType:disableTransmitterShadowEpochs}
Disable transmitter epochs during eclipse.
With proper attitude modeling (see \program{SimulateStarCameraGnss}) this is usually not necessary.
)";
#endif

/***********************************************/

#include "config/config.h"
#include "base/kepler.h"
#include "classes/eclipse/eclipse.h"
#include "gnss/gnssTransceiverSelector/gnssTransceiverSelector.h"
#include "gnss/gnssProcessingStep/gnssProcessingStep.h"

/***** CLASS ***********************************/

/** @brief GNSS processing step: DisableTransmitterShadowEpochs.
* @ingroup gnssProcessingStepGroup
* @see GnssProcessingStep */
class GnssProcessingStepDisableTransmitterShadowEpochs : public GnssProcessingStepBase
{
  GnssTransceiverSelectorPtr selectTransmitters;
  Bool                       disableShadowEpochs, disablePostShadowEpochs;
  EphemeridesPtr             ephemerides;
  EclipsePtr                 eclipse;

public:
  GnssProcessingStepDisableTransmitterShadowEpochs(Config &config);
  void process(GnssProcessingStep::State &state) override;
  Bool expectInitializedParameters() const override {return FALSE;}
};

/***********************************************/

inline GnssProcessingStepDisableTransmitterShadowEpochs::GnssProcessingStepDisableTransmitterShadowEpochs(Config &config)
{
  try
  {
    readConfig(config, "selectTransmitters",              selectTransmitters,      Config::MUSTSET,  "",  "");
    readConfig(config, "disableShadowEpochs",             disableShadowEpochs,     Config::DEFAULT,  "1", "disable epochs if satellite is in Earth's/Moon's shadow");
    readConfig(config, "disablePostShadowRecoveryEpochs", disablePostShadowEpochs, Config::DEFAULT,  "1", "disable epochs if satellite is in post-shadow recovery maneuver for GPS block IIA");
    readConfig(config, "ephemerides",                     ephemerides,             Config::MUSTSET,  "",  "");
    readConfig(config, "eclipse",                         eclipse,                 Config::MUSTSET,  "",  "eclipse model used to determine if a satellite is in Earth's shadow");
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline void GnssProcessingStepDisableTransmitterShadowEpochs::process(GnssProcessingStep::State &state)
{
  try
  {
    logStatus<<"=== disable transmitter epochs in eclipse ==================="<<Log::endl;
    if(state.normalEquationInfo.isEachReceiverSeparately)
    {
      logWarning<<"DisableTransmitterShadowEpochs is not allowed in single receiver loop"<<Log::endl;
      return;
    }

    UInt countEpochs = 0;
    auto selectedTransmitters = selectTransmitters->select(state.gnss->transmitters);
    for(UInt idTrans=0; idTrans<state.gnss->transmitters.size(); idTrans++)
      if(selectedTransmitters.at(idTrans) && state.gnss->transmitters.at(idTrans)->useable())
      {
        auto trans = state.gnss->transmitters.at(idTrans);

        // 30 minute post-shadow recovery time for block IIA satellites
        Time recoveryTime;
        if(trans->info.antenna.at(trans->info.findAntenna(state.gnss->times.at(0))).name.find(std::string("BLOCK IIA")) != std::string::npos)
          recoveryTime = seconds2time(30*60);

        // Look for a potential shadow exit before the start of the interval that would result in the
        // satellite already performing a post-shadow recovery maneuver at the start of the interval
        Time timeShadowExit;
        if(disablePostShadowEpochs && recoveryTime.seconds() > 0)
        {
          Kepler kepler(state.gnss->times.at(0), trans->positionCoM(state.gnss->times.at(0)), trans->velocity(state.gnss->times.at(0)));
          for(UInt i=0; i<recoveryTime.seconds()/60; i++)
          {
            const Time time = state.gnss->times.at(0) - seconds2time(i*60.);
            if(eclipse->factor(time, kepler.position(time), ephemerides) < 0.5)
            {
              timeShadowExit = time;
              break;
            }
          }
        }

        Double factorPreviousEpoch = 1.0;
        for(UInt idEpoch=0; idEpoch<state.gnss->times.size(); idEpoch++)
          if(trans->useable(idEpoch))
          {
            const Double factor = eclipse->factor(state.gnss->times.at(idEpoch), trans->positionCoM(state.gnss->times.at(idEpoch)), ephemerides);
            if((factorPreviousEpoch < 0.5) && (factor >= 0.5))
              timeShadowExit = state.gnss->times.at(idEpoch);

            // set satellite unuseable during shadow crossing and post-shadow recovery maneuver
            if((disableShadowEpochs && factor < 0.5) || (disablePostShadowEpochs && state.gnss->times.at(idEpoch) < timeShadowExit+recoveryTime))
            {
              trans->disable(idEpoch, "during shadow crossing and post-shadow recovery maneuver");
              countEpochs++;
            }

            factorPreviousEpoch = factor;
          } // for(idEpoch)
      }

    if(countEpochs)
      state.changedNormalEquationInfo = TRUE;
    logInfo<<"  "<<countEpochs<<" disabled transmitter epochs"<<Log::endl;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

#endif
