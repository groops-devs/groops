/***********************************************/
/**
* @file gnssProcessingStepWriteResiduals.h
*
* @brief GNSS processing step: WriteResiduals.
*
* @author Torsten Mayer-Guerr
* @date 2021-09-05
*
*/
/***********************************************/

#ifndef __GROOPS_GNSSPROCESSINGSTEPWRITERESIDUALS__
#define __GROOPS_GNSSPROCESSINGSTEPWRITERESIDUALS__

// Latex documentation
#ifdef DOCSTRING_GnssProcessingStep
static const char *docstringGnssProcessingStepWriteResiduals = R"(
\subsection{WriteResiduals}\label{gnssProcessingStepType:writeResiduals}
Writes the \file{observation residuals}{instrument} for all
\configClass{selectReceivers}{gnssTransceiverSelectorType}.
For for each station a file is written. The file name is interpreted as
a template with the variable \verb|{station}| being replaced by the station name.
)";
#endif

/***********************************************/

#include "config/config.h"
#include "gnss/gnssTransceiverSelector/gnssTransceiverSelector.h"
#include "gnss/gnssProcessingStep/gnssProcessingStep.h"

/***** CLASS ***********************************/

/** @brief GNSS processing step: WriteResiduals.
* @ingroup gnssProcessingStepGroup
* @see GnssProcessingStep */
class GnssProcessingStepWriteResiduals : public GnssProcessingStepBase
{
  GnssTransceiverSelectorPtr selectReceivers;
  FileName                   fileNameResiduals;

public:
  GnssProcessingStepWriteResiduals(Config &config);
  void process(GnssProcessingStep::State &state) override;
};

/***********************************************/

inline GnssProcessingStepWriteResiduals::GnssProcessingStepWriteResiduals(Config &config)
{
  try
  {
    readConfig(config, "selectReceivers",     selectReceivers,   Config::MUSTSET, "", "subset of used stations");
    readConfig(config, "outputfileResiduals", fileNameResiduals, Config::MUSTSET, "output/residuals_{loopTime:%D}.{station}.dat", "variable {station} available");
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssProcessingStepWriteResiduals::process(GnssProcessingStep::State &state)
{
  try
  {
    auto selectedReceivers = selectReceivers->select(state.gnss->receivers);
    VariableList fileNameVariableList;
    addVariable("station", "****", fileNameVariableList);
    logStatus<<"write residuals to file <"<<fileNameResiduals(fileNameVariableList)<<">"<<Log::endl;
    for(auto recv : state.gnss->receivers)
      if(selectedReceivers.at(recv->idRecv()) && state.normalEquationInfo.estimateReceiver.at(recv->idRecv()) && recv->isMyRank())
      {
        GnssReceiverArc arc;
        for(UInt idEpoch : state.normalEquationInfo.idEpochs)
          if(recv->useable(idEpoch))
          {
            // get types
            std::vector<GnssType> types;
            for(UInt idTrans=0; idTrans<recv->idTransmitterSize(idEpoch); idTrans++)
              if(recv->observation(idTrans, idEpoch) && state.gnss->transmitters.at(idTrans)->useable(idEpoch))
                for(UInt idType=0; idType<recv->observation(idTrans, idEpoch)->size(); idType++)
                  if(!recv->observation(idTrans, idEpoch)->at(idType).type.isInList(types))
                    types.push_back(recv->observation(idTrans, idEpoch)->at(idType).type & ~(GnssType::PRN+GnssType::FREQ_NO));
            std::sort(types.begin(), types.end());
            if(!types.size())
              continue;

            GnssReceiverEpoch epoch;
            GnssType system = GnssType::SYSTEM;
            for(UInt idType=0; idType<types.size(); idType++)
            {
              if(types.at(idType) != system)
              {
                system = types.at(idType) & GnssType::SYSTEM;
                epoch.obsType.push_back( GnssType::AZIMUT    + GnssType::L1 + system );
                epoch.obsType.push_back( GnssType::ELEVATION + GnssType::L1 + system );
                epoch.obsType.push_back( GnssType::AZIMUT    + GnssType::L2 + system );
                epoch.obsType.push_back( GnssType::ELEVATION + GnssType::L2 + system );
                epoch.obsType.push_back( GnssType::IONODELAY + system );
              }
              // residuals, redundancy, sigma
              epoch.obsType.insert(epoch.obsType.end(), {types.at(idType), types.at(idType), types.at(idType)});
            }

            for(UInt idTrans=0; idTrans<recv->idTransmitterSize(idEpoch); idTrans++)
              if(recv->observation(idTrans, idEpoch) && state.gnss->transmitters.at(idTrans)->useable(idEpoch))
              {
                const GnssObservation &obs = *recv->observation(idTrans, idEpoch);
                const GnssObservationEquation eqn(obs, *recv, *state.gnss->transmitters.at(idTrans),
                                                  state.gnss->funcRotationCrf2Trf, state.gnss->funcReduceModels, idEpoch, FALSE, {});
                const GnssType prn = obs.at(0).type & (GnssType::SYSTEM + GnssType::PRN + GnssType::FREQ_NO);
                UInt idType = std::distance(epoch.obsType.begin(), std::find(epoch.obsType.begin(), epoch.obsType.end(), prn));
                if(idType >= epoch.obsType.size())
                  continue;

                epoch.time = eqn.timeRecv;
                epoch.satellite.push_back(prn);
                epoch.observation.insert(epoch.observation.end(), {eqn.azimutRecvAnt, eqn.elevationRecvAnt, eqn.azimutTrans, eqn.elevationTrans, obs.STEC});
                idType += 5;

                for(; (idType<epoch.obsType.size()) && (epoch.obsType.at(idType) == prn); idType+=3)
                {
                  epoch.observation.insert(epoch.observation.end(), {0., 0., 1.});
                  for(UInt i=0; i<obs.size(); i++)
                    if(obs.at(i).type == epoch.obsType.at(idType))
                    {
                      epoch.observation.at(epoch.observation.size()-3) = obs.at(i).residuals;
                      epoch.observation.at(epoch.observation.size()-2) = obs.at(i).redundancy;
                      epoch.observation.at(epoch.observation.size()-1) = obs.at(i).sigma/obs.at(i).sigma0;
                      break;
                    }
                }
              } // for(idTrans)

            if(epoch.satellite.size())
              arc.push_back(epoch);
          } // for(idEpoch)

        VariableList fileNameVariableList;
        addVariable("station", recv->name(), fileNameVariableList);
        InstrumentFile::write(fileNameResiduals(fileNameVariableList), arc);
      } // for(recv)
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

#endif
