/***********************************************/
/**
* @file gnssReceiver.cpp
*
* @brief GNSS receiver.
*
* @author Torsten Mayer-Guerr
* @author Norbert Zehentner
* @author Sebastian Strasser
* @date 2013-06-28
*
*/
/***********************************************/

#include "base/import.h"
#include "parser/dataVariables.h"
#include "files/fileMatrix.h"
#include "files/fileInstrument.h"
#include "inputOutput/logging.h"
#include "gnss/gnss.h"
#include "gnss/gnssObservation.h"
#include "gnss/gnssTransmitter.h"
#include "gnss/gnssReceiver.h"
#include "gnss/gnssParametrizationAmbiguities.h"
#include "misc/varianceComponentEstimation.h"

/***********************************************/

Gnss::Receiver::Receiver()  : _gnss(nullptr), _idRecv(NULLINDEX), useAtAll(FALSE), observationSampling(0), huber(2.5), huberPower(1.5)
{
}

/***********************************************/

Gnss::Receiver::~Receiver()
{
  free();
}

/***********************************************/

void Gnss::Receiver::free()
{
  try
  {
    // delete observations
//     for(UInt idTrans=0; idTrans<observations.size(); idTrans++)
//       for(UInt idEpoch=0; idEpoch<observations.at(idTrans).size(); idEpoch++)
//         delete observations.at(idTrans).at(idEpoch);
obsMem.clear();
    observations_.clear();
    track.clear();

    sigma0Type.clear();
    sigma0Factor.clear();

    clk.clear();
    use.clear();
    useAtAll = FALSE;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Gnss &Gnss::Receiver::gnss() const
{
  if(_gnss == nullptr)
    throw(Exception("Receiver is not registered in Gnss class"));
  return *_gnss;
}

/***********************************************/

void Gnss::Receiver::initInterval(const std::vector<Time> &times)
{
  try
  {
    this->useAtAll = TRUE;
    this->times = times;
    this->clk.clear(); this->clk.resize(times.size(), 0);
    this->use.clear(); this->use.resize(times.size(), TRUE);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void Gnss::Receiver::disable(UInt idEpoch)
{
  try
  {
    if(idEpoch == NULLINDEX)
    {
      useAtAll = FALSE;
obsMem.clear();
obsMem.shrink_to_fit();
      for(auto &obsEpoch : observations_)
        for(auto &obs : obsEpoch)
          obs = nullptr;
    }
    else
    {
      use.at(idEpoch) = FALSE;
      if(idEpoch < observations_.size())
        observations_.at(idEpoch).clear();
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Bool Gnss::Receiver::isReceiverEstimable(Double minEstimableEpochsRatio) const
{
  try
  {
    if(!useable() || observationSampling == 0.)
      return FALSE;

    UInt countEpoch = 0;
    for(UInt idEpoch=0; idEpoch<times.size(); idEpoch++)
      if(useable(idEpoch))
        countEpoch++;

    const UInt minEpochCount = static_cast<UInt>(std::round(minEstimableEpochsRatio * gnss().times.size() * medianSampling(gnss().times).seconds()/observationSampling));
    return countEpoch >= minEpochCount;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Gnss::Observation *Gnss::Receiver::observation(UInt idTrans, UInt idEpoch) const
{
  if((idEpoch < observations_.size()) && (idTrans < observations_.at(idEpoch).size()))
    return observations_[idEpoch][idTrans];
  return nullptr;
}

/***********************************************/

UInt Gnss::Receiver::countObservations(UInt idTrans, UInt idEpochStart, UInt idEpochEnd) const
{
  try
  {
    if(!useable())
     return 0;

    UInt count = 0;
    const UInt idTransStart = ((idTrans!=NULLINDEX) ? idTrans : 0);
    const UInt idTransEnd   = ((idTrans!=NULLINDEX) ? idTrans : MAX_UINT);
    if(idEpochStart==NULLINDEX) idEpochStart = 0;
    if(idEpochEnd  ==NULLINDEX) idEpochEnd   = MAX_UINT;

    for(UInt idEpoch=idEpochStart; (idEpoch<=idEpochEnd) && (idEpoch<idEpochSize()); idEpoch++)
      if(useable(idEpoch))
        for(UInt idTrans=idTransStart; (idTrans<=idTransEnd) && (idTrans<idTransmitterSize(idEpoch)); idTrans++)
          if(observation(idTrans, idEpoch))
            count++;

    return count;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void Gnss::Receiver::countObservationsPerType(UInt idTrans, UInt idEpochStart, UInt idEpochEnd, std::vector<GnssType> &types, std::vector<UInt> &count) const
{
  try
  {
    if(!useable())
     return;

    const UInt idTransStart = ((idTrans!=NULLINDEX) ? idTrans : 0);
    const UInt idTransEnd   = ((idTrans!=NULLINDEX) ? idTrans : MAX_UINT);
    if(idEpochStart==NULLINDEX) idEpochStart = 0;
    if(idEpochEnd  ==NULLINDEX) idEpochEnd   = MAX_UINT;

    for(UInt idEpoch=idEpochStart; (idEpoch<=idEpochEnd) && (idEpoch<idEpochSize()); idEpoch++)
      if(useable(idEpoch))
        for(UInt idTrans=idTransStart; (idTrans<=idTransEnd) && (idTrans<idTransmitterSize(idEpoch)); idTrans++)
        {
          auto obs = observation(idTrans, idEpoch);
          if(obs)
            for(UInt idType=0; idType<obs->size(); idType++)
            {
              UInt idx = GnssType::index(types, obs->at(idType).type);
              if(idx == NULLINDEX)
              {
                idx = types.size();
                types.push_back(obs->at(idType).type & ~GnssType::PRN);
                count.push_back(0);
              }
              count.at(idx)++;
            }
        }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

std::vector<GnssType> Gnss::Receiver::observationsTypes(UInt idTrans) const
{
  try
  {
    std::vector<GnssType> types;
    const UInt idTransStart = ((idTrans!=NULLINDEX) ? idTrans : 0);
    const UInt idTransEnd   = ((idTrans!=NULLINDEX) ? idTrans : MAX_UINT);
    if(useable())
      for(UInt idEpoch=0; idEpoch<idEpochSize(); idEpoch++)
        if(useable(idEpoch))
          for(UInt idTrans=idTransStart; (idTrans<=idTransEnd) && (idTrans<idTransmitterSize(idEpoch)); idTrans++)
          {
            auto obs = observation(idTrans, idEpoch);
            if(obs)
              for(UInt idType=0; idType<obs->size(); idType++)
                if(GnssType::index(types, obs->at(idType).type) == NULLINDEX)
                  types.push_back(obs->at(idType).type & ~GnssType::PRN);
          }
    std::sort(types.begin(), types.end());
    return types;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void Gnss::Receiver::residualsStatistics(UInt idTrans, const std::vector<UInt> &idEpochs,
                                         std::vector<GnssType> &types, std::vector<Double> &ePe, std::vector<Double> &redundancy,
                                         std::vector<UInt> &obsCount, std::vector<UInt> &outlierCount) const
{
  try
  {
    if(!useable())
     return;

    const UInt idTransStart = ((idTrans!=NULLINDEX) ? idTrans : 0);
    const UInt idTransEnd   = ((idTrans!=NULLINDEX) ? idTrans : MAX_UINT);

    for(UInt idEpoch : idEpochs)
      if(useable(idEpoch))
        for(UInt idTrans=idTransStart; (idTrans<=idTransEnd) && (idTrans<idTransmitterSize(idEpoch)); idTrans++)
          if(observation(idTrans, idEpoch))
          {
            const Observation &obs = *observation(idTrans, idEpoch);
            for(UInt idType=0; idType<obs.size(); idType++)
              if(obs.at(idType).sigma0>0)
              {
                UInt idx = GnssType::index(types, obs.at(idType).type);
                if(idx == NULLINDEX) // new type?
                {
                  idx = types.size();
                  types.push_back(obs.at(idType).type & GnssType::NOPRN); // without PRN
                  ePe.push_back(0);
                  redundancy.push_back(0);
                  obsCount.push_back(0);
                  outlierCount.push_back(0);
                }

                ePe.at(idx)        += std::pow(obs.at(idType).residuals/obs.at(idType).sigma, 2);
                redundancy.at(idx) += obs.at(idType).redundancy;
                obsCount.at(idx)++;
                if(obs.at(idType).sigma > obs.at(idType).sigma0)
                  outlierCount.at(idx)++;
              }
          } // for(idTrans, idEpoch)
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void Gnss::Receiver::computeWeightsFromResiduals(const std::vector<UInt> &idEpochs, Gnss::Receiver::WeightingType computeWeights, Gnss::Receiver::WeightingType adjustSigma0)
{
  try
  {
    if(adjustSigma0 != WeightingType::NONE)
    {
      std::vector<GnssType> types;
      std::vector<Double>   ePe, redundancy;
      std::vector<UInt>     obsCount, outlierCount;

      // same sigma for all phase types
      if((adjustSigma0 == WeightingType::GROUPED) || (adjustSigma0 == WeightingType::GROUPEDPHASES))
      {
        types.push_back(GnssType::PHASE);
        ePe.push_back(0);
        redundancy.push_back(0);
        obsCount.push_back(0);
        outlierCount.push_back(0);
      }
      // same sigma also for all code types
      if(adjustSigma0 == WeightingType::GROUPED)
      {
        types.push_back(GnssType::RANGE);
        ePe.push_back(0);
        redundancy.push_back(0);
        obsCount.push_back(0);
        outlierCount.push_back(0);
      }

      residualsStatistics(NULLINDEX, idEpochs, types, ePe, redundancy, obsCount, outlierCount);

      // compute accuracies
      std::vector<Double> factor(types.size(), 1.0);
      for(UInt i=0; i<types.size(); i++)
        if(obsCount.at(i))
        {
          factor.at(i) = Vce::standardDeviation(ePe.at(i), redundancy.at(i), huber, huberPower);
          if(std::isnan(factor.at(i)) || (factor.at(i) <= 0))
          {
            logWarning<<name()<<":"<<types.at(i).str()<<": cannot determine factor = "<<ePe.at(i)<<"/"<<redundancy.at(i)<<Log::endl;
            factor.at(i) = 10;
          }
        }

      // adjust sigma and sigma0
      for(UInt idEpoch=0; idEpoch<idEpochSize(); idEpoch++)
        for(UInt idTrans=0; idTrans<idTransmitterSize(idEpoch); idTrans++)
          if(observation(idTrans, idEpoch))
          {
            Observation &obs = *observation(idTrans, idEpoch);
            for(UInt idType=0; idType<obs.size(); idType++)
              if(obs.at(idType).sigma0>0)
              {
                const UInt idx = GnssType::index(types, obs.at(idType).type);
                if(idx == NULLINDEX)
                  continue;
                obs.at(idType).sigma  *= factor.at(idx);
                obs.at(idType).sigma0 *= factor.at(idx);
              } // for(idType)
          } // for(idTrans, idEpoch)


      // apply old factors to store total factor
      for(UInt idType=0; idType<types.size(); idType++)
      {
        const UInt idx = GnssType::index(sigma0Type, types.at(idType));
        if(idx != NULLINDEX)
          factor.at(idType) *= sigma0Factor.at(idx);
      }
      sigma0Type   = types;
      sigma0Factor = factor;
    } // if(adjustSigma0)

    // ===========================================================

    // adjust sigmas of observations (downweight outliers)
    if(computeWeights != WeightingType::NONE)
    {
      for(UInt idEpoch : idEpochs)
        if(useable(idEpoch))
          for(UInt idTrans=0; idTrans<idTransmitterSize(idEpoch); idTrans++)
            if(observation(idTrans, idEpoch))
            {
              Observation &obs = *observation(idTrans, idEpoch);
              std::vector<GnssType> types;      // ranges, phases, ...
              std::vector<Double>   ePe;        // quadsum of residuals of each group
              std::vector<Double>   redundancy; // redundancy of each group

              // same sigma for all phase types
              if((computeWeights == WeightingType::GROUPED) || (computeWeights == WeightingType::GROUPEDPHASES))
              {
                types.push_back(GnssType::PHASE);
                ePe.push_back(0);
                redundancy.push_back(0);
              }
              // same sigma also for all code types
              if(computeWeights == WeightingType::GROUPED)
              {
                types.push_back(GnssType::RANGE);
                ePe.push_back(0);
                redundancy.push_back(0);
              }

              // assign observations to groups
              for(UInt idType=0; idType<obs.size(); idType++)
              {
                UInt idx = GnssType::index(types, obs.at(idType).type);
                if(idx == NULLINDEX) // new type
                {
                  idx = types.size();
                  types.push_back(obs.at(idType).type);
                  ePe.push_back(0.);
                  redundancy.push_back(0.);
                }
                ePe.at(idx)        += std::pow(obs.at(idType).residuals/obs.at(idType).sigma0, 2);
                redundancy.at(idx) += obs.at(idType).redundancy;
              }

              // adjust sigmas
              for(UInt idType=0; idType<obs.size(); idType++)
              {
                const UInt idx = GnssType::index(types, obs.at(idType).type);
                obs.at(idType).sigma = obs.at(idType).sigma0;
                if(redundancy.at(idx) > 0.1) // redundancy: it is possible to estimate sigma?
                {
                  const Double s = std::sqrt(ePe.at(idx)/redundancy.at(idx));
                  if(s > huber) // outlier? ==> change weight
                    obs.at(idType).sigma *= std::pow(s/huber, huberPower);
                }
              } // for(idType)
          } // for(idTrans, idEpoch)
    } // if(computeWeights)
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Double Gnss::Receiver::sigmaFactor(GnssType type) const
{
  try
  {
    const UInt idx = GnssType::index(sigma0Type, type);
    if(idx != NULLINDEX)
      return sigma0Factor.at(idx);
    return 1;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

GnssReceiverArc Gnss::Receiver::residuals(const std::vector<UInt> &idEpochs) const
{
  try
  {
    GnssReceiverArc arc;
    for(UInt idEpoch : idEpochs)
      if(useable(idEpoch))
      {
        GnssReceiverEpoch epoch;

        std::vector<GnssType> types;
        std::vector<UInt>     count;
        countObservationsPerType(NULLINDEX, idEpoch, idEpoch, types, count);

        // remove GLONASS frequency number and consequential duplicate types
        for(auto &&type : types)
          type &= ~GnssType::FREQ_NO;
        std::sort(types.begin(), types.end());
        types.erase(std::unique(types.begin(), types.end()), types.end());

        if(types.size() == 0)
          continue;
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
            epoch.obsType.push_back( GnssType::IONODELAY + system );
            epoch.obsType.push_back( GnssType::IONODELAY + system );
          }
          epoch.obsType.push_back(types.at(idType)); // residuals
          epoch.obsType.push_back(types.at(idType)); // redundancy
          epoch.obsType.push_back(types.at(idType)); // sigma
        }

        for(UInt idTrans=0; idTrans<idTransmitterSize(idEpoch); idTrans++)
          if(gnss().transmitter.at(idTrans)->useable(idEpoch) && observation(idTrans, idEpoch))
          {
            const Observation &obs = *observation(idTrans, idEpoch);
            const ObservationEquation eqn(obs, *this, *gnss().transmitter.at(idTrans), gnss().ionosphere, idEpoch, FALSE, {});

            GnssType prn = obs.at(0).type & (GnssType::SYSTEM + GnssType::PRN + GnssType::FREQ_NO);
            UInt idType = 0;
            while((idType<epoch.obsType.size()) && (epoch.obsType.at(idType) != prn))
              idType++;
            if(idType>=epoch.obsType.size())
              continue;

            epoch.time = eqn.timeRecv;
            epoch.satellite.push_back( prn );
            epoch.observation.push_back( eqn.azimutRecvAnt );
            epoch.observation.push_back( eqn.elevationRecvAnt );
            epoch.observation.push_back( eqn.azimutTrans );
            epoch.observation.push_back( eqn.elevationTrans );
            idType += 4;

            epoch.observation.push_back( obs.dSTEC );
            epoch.observation.push_back( obs.dSTEC );
            epoch.observation.push_back( obs.dSTEC );
            idType += 3;

            for(; (idType<epoch.obsType.size()) && (epoch.obsType.at(idType) == prn); idType+=3)
            {
              Bool found = FALSE;
              for(UInt i=0; i<obs.size(); i++)
                if(obs.at(i).type == epoch.obsType.at(idType))
                {
                  epoch.observation.push_back( obs.at(i).residuals );
                  epoch.observation.push_back( obs.at(i).redundancy );
                  epoch.observation.push_back( obs.at(i).sigma/obs.at(i).sigma0 );
                  found = TRUE;
                  break;
                }
              if(!found)
              {
                epoch.observation.push_back( 0 );
                epoch.observation.push_back( 0 );
                epoch.observation.push_back( 1 );
              }
            }
          } // for(idTrans)

        if(epoch.satellite.size())
          arc.push_back(epoch);
      } // for(idEpoch)

    return arc;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void Gnss::Receiver::readObservations(InstrumentFile &fileReceiver,
                                      const std::vector<GnssType> &useType, const std::vector<GnssType> &ignoreType,
                                      const Time &timeMargin)
{
  try
  {
    std::vector<Time> observationTimes;
    UInt idEpoch = 0;
    for(UInt arcNo=0; arcNo<fileReceiver.arcCount(); arcNo++)
    {
      GnssReceiverArc arc = fileReceiver.readArc(arcNo);
      for(UInt arcEpoch=0; arcEpoch<arc.size(); arcEpoch++)
      {
        // search time slot
        while((idEpoch < times.size()) && (times.at(idEpoch)+timeMargin < arc.at(arcEpoch).time))
          idEpoch++;
        if(idEpoch >= times.size())
          break;
        if(arc.at(arcEpoch).time+timeMargin < times.at(idEpoch))
          continue;
        times.at(idEpoch) = arc.at(arcEpoch).time;
//         clk.at(idEpoch)   = arc.at(arcEpoch).clockError;
        observationTimes.push_back(arc.at(arcEpoch).time);

        // create observation class for each satellite
        UInt idObs  = 0;
        for(UInt k=0; k<arc.at(arcEpoch).satellite.size(); k++)
        {
          Observation *obs = new Observation();

          // find list of observation types for this satellite
          GnssType satType = arc.at(arcEpoch).satellite.at(k);
          UInt idType = 0;
          while(arc.at(arcEpoch).obsType.at(idType) != satType)
            idType++;

          for(; (idType<arc.at(arcEpoch).obsType.size()) && (arc.at(arcEpoch).obsType.at(idType)==satType); idType++, idObs++)
          {
            if(arc.at(arcEpoch).observation.at(idObs) != 0)
            {
              GnssType type = arc.at(arcEpoch).obsType.at(idType) + satType;
              // remove GLONASS frequency number
              if((type == GnssType::GLONASS) && !((type == GnssType::G1) || (type == GnssType::G2)))
                type.setFrequencyNumber(9999);

              // check completeness
              if(type.hasWildcard(GnssType::TYPE + GnssType::FREQUENCY + GnssType::SYSTEM + GnssType::PRN))
              {
                logWarning<<name()<<" "<<type.str()<<" is not complete"<<Log::endl;
                continue;
              }
              if((type.frequencyNumber() == 9999) && (type == GnssType::GLONASS) && ((type == GnssType::G1) || (type == GnssType::G2)))
              {
                logWarning<<name()<<" "<<type.str()<<": GLONASS frequency number not set"<<Log::endl;
                continue;
              }
              if((type.frequencyNumber() != 9999) && !((type == GnssType::GLONASS) && ((type == GnssType::G1) || (type == GnssType::G2))))
              {
                logWarning<<name()<<" "<<type.str()<<": GLONASS frequency number is set"<<Log::endl;
                continue;
              }

              Bool use = (useType.size()==0) ? TRUE : FALSE;
              if(GnssType::index(useType, type) != NULLINDEX)
                use = TRUE;
              if(GnssType::index(ignoreType, type) != NULLINDEX)
                use = FALSE;
              if(use)
                obs->push_back(SingleObservation(type, arc.at(arcEpoch).observation.at(idObs)));
            }
          }
          if(obs->size() == 0)
          {
            delete obs;
            continue;
          }

          // search transmitter index for satellite number (PRN)
          UInt idTrans = 0;
          while((idTrans<gnss().transmitter.size()) && (gnss().transmitter.at(idTrans)->PRN() != satType))
            idTrans++;
          if(idTrans >= gnss().transmitter.size()) // observation without transmitter
          {
            delete obs;
            continue;
          }

          if(observations_.size() <= idEpoch)
            observations_.resize(idEpoch+1);
          if(observations_.at(idEpoch).size() <= idTrans)
            observations_.at(idEpoch).resize(idTrans+1, nullptr);
          if(observations_[idEpoch][idTrans])
            logWarning<<"observation already exists"<<Log::endl;
          std::swap(observations_[idEpoch][idTrans], obs);
          delete obs;
        } // for(satellite)
        idEpoch++;
      } // for(arcEpoch)
      if(idEpoch >= times.size())
        break;
    } // for(arcNo)

    // median sampling
    // ---------------
    observationSampling = medianSampling(observationTimes).seconds();

    // Copy observations to a continuous memory block
    // ----------------------------------------------
    // count observations
    UInt count = 0;
    for(const auto &obsEpoch : observations_)
      for(const Observation *obs : obsEpoch)
        if(obs)
          count++;
    obsMem.resize(count);
    obsMem.shrink_to_fit();
    // copy
    count = 0;
    for(auto &obsEpoch : observations_)
      for(Observation *&obs : obsEpoch)
        if(obs)
        {
          obsMem.at(count) = *obs;
          obsMem.at(count).shrink_to_fit();
          delete obs;
          obs = &obsMem.at(count++);
        }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void Gnss::Receiver::initObservation(Gnss::AnalysisType analysisType)
{
  try
  {
    for(UInt idTrans=0; idTrans<gnss().transmitter.size(); idTrans++)
    {
      Double phaseWindup = 0;
      for(UInt idEpoch=0; idEpoch<idEpochSize(); idEpoch++)
      {
        Observation *obs = observation(idTrans, idEpoch);
        if(obs)
        {
          std::vector<GnssType> types;
          if(!(obs->init(*this, *gnss().transmitter.at(idTrans), gnss().ionosphere, idEpoch, phaseWindup) && obs->observationList(analysisType, types)))
            deleteObservation(idTrans, idEpoch);
        }
      } // for idEpoch
    } // for idTrans
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void Gnss::Receiver::deleteObservation(UInt idTrans, UInt idEpoch)
{
  try
  {
    if(observation(idTrans, idEpoch))
      observations_[idEpoch][idTrans] = nullptr;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Bool Gnss::Receiver::deleteUndefinedObservations()
{
  try
  {
    std::map<GnssType, UInt> removedObs;
    for(UInt idEpoch = 0; idEpoch < times.size(); idEpoch++)
      for(UInt idTrans = 0; idTrans < gnss().transmitter.size(); idTrans++)
      {
        Observation *obs = observation(idTrans, idEpoch);
        if(obs == nullptr)
          continue;

        const std::vector<GnssType> receiverTypes    = definedTypes(idEpoch);
        const std::vector<GnssType> transmitterTypes = gnss().transmitter.at(idTrans)->definedTypes(idEpoch);
        for(UInt idObs = obs->size(); idObs --> 0; )
        {
          // check against receiver types
          if(receiverTypes.size() && GnssType::index(receiverTypes, obs->at(idObs).type) == NULLINDEX)
          {
            removedObs[obs->at(idObs).type]++;
            obs->erase(idObs);
            continue;
          }

          // check against transmitter types
          if(transmitterTypes.size())
            for(const auto &type : gnss().replaceCompositeSignals({obs->at(idObs).type})) // ATTENTION: S*X and D*X are not replaced and will thus be removed
              if(GnssType::index(transmitterTypes, type) == NULLINDEX)
              {
                removedObs[obs->at(idObs).type]++;
                obs->erase(idObs);
                break;
              }
        }

        if(!obs->size())
          deleteObservation(idTrans, idEpoch);
      }

    if(removedObs.size())
    {
      std::stringstream ss;
      for(const auto &obs : removedObs)
        ss << " " << obs.first.str() << "=" << obs.second;
      logWarning<<name()<<": removed undefined observations"<<ss.str()<<Log::endl;
    }

    return removedObs.size();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Bool Gnss::Receiver::deleteObservationsOfInestimableEpochs(Gnss::AnalysisType analysisType)
{
  try
  {
    Bool epochRemoved = FALSE;
    for(UInt idEpoch=0; idEpoch<times.size(); idEpoch++)
      if(!isEpochEstimable(analysisType, idEpoch))
      {
        epochRemoved = countObservations(NULLINDEX/*idTrans*/, idEpoch, idEpoch);
        disable(idEpoch);
      }
    return epochRemoved;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void Gnss::Receiver::estimateInitialClockErrorFromCodeObservations(Double maxPosDiff, Bool estimateKinematicPosition)
{
  try
  {
    Vector posDiff(estimateKinematicPosition ? times.size() : 1);
    Vector clockDiff(times.size()); // in meters
    for(UInt iter=0; iter<10; iter++)
    {
      // setup observation equations: position, clock
      // --------------------------------------------
      std::vector<Matrix> listl, listA;
      std::vector<Matrix> listlFull, listAFull, listBFull;
      std::vector<UInt>   listEpoch, listObsCount;

      for(UInt idEpoch=0; idEpoch<times.size(); idEpoch++)
        if(isEpochEstimable(ANALYSIS_CODE, idEpoch))
        {
          // count observations and setup observation equations for each transmitter
          UInt obsCount = 0;
          std::vector<ObservationEquation> eqnList;
          std::vector<GnssType> types, systems;
          for(UInt idTrans=0; idTrans<idTransmitterSize(idEpoch); idTrans++)
          {
            Observation *obs = observation(idTrans, idEpoch);
            if(obs && obs->observationList(ANALYSIS_CODE, types))
            {
              eqnList.emplace_back(*obs, *this, *gnss().transmitter.at(idTrans), gnss().ionosphere, idEpoch, TRUE, types);

              eliminationParameter(eqnList.back().B, eqnList.back().A, eqnList.back().l); // eliminate TEC
              obsCount += eqnList.back().l.rows();
              if(GnssType::index(systems, types.at(0)) == NULLINDEX)
                systems.push_back(types.at(0) & GnssType::SYSTEM);
            }
          }

          // setup combined observation equations
          const UInt countStaticPosition = (estimateKinematicPosition ? 0 : 3);
          Vector l(obsCount);
          Matrix A(obsCount, countStaticPosition);                  // pos (if static)
          Matrix B(obsCount, systems.size()+3-countStaticPosition); // clock (per system) + pos (if kinematic)
          UInt idx = 0;
          for(const auto &eqn : eqnList)
          {
            const UInt count = eqn.l.rows();
            copy(eqn.l, l.row(idx, count));
            copy(eqn.A.column(ObservationEquation::idxClockRecv),  B.slice(idx, GnssType::index(systems, eqn.types.at(0)), count, 1)); // clock
            copy(eqn.A.column(ObservationEquation::idxPosRecv, 3), (estimateKinematicPosition ? B.slice(idx, systems.size(), count, 3) : A.slice(idx, 0, count, 3))); // pos
            idx += count;
          }

          listEpoch.push_back(idEpoch);
          listlFull.push_back(l);
          listAFull.push_back(A);
          listBFull.push_back(B);

          if(A.columns())
          {
            eliminationParameter(B, A, l);
            listl.push_back(l);
            listA.push_back(A);
            listObsCount.push_back(obsCount-B.columns());
          }
        } // for(idEpoch)
      if(!listEpoch.size())
        break;

      // estimate static parameters
      // --------------------------
      Vector sigma0(listl.size(), 100.);
      Vector sigma = sigma0;
      const Vector x = robustLeastSquares(huber, huberPower, listA, listl, listObsCount, sigma0, sigma);
      if(!estimateKinematicPosition)
        posDiff(0) = norm(x.row(0,3));

      // reconstruct kinematic parameters
      // --------------------------------
      for(UInt i=0; i<listEpoch.size(); i++)
      {
        if(x.size()) // update with static parameters
          matMult(-1, listAFull.at(i), x, listlFull.at(i));

        Vector sigma0(listlFull.at(i).size(), 1.);
        Vector sigma = sigma0;
        const Vector y = robustLeastSquares(huber, huberPower, listBFull.at(i), listlFull.at(i), sigma0, sigma);

        // set receiver clock error as mean of all system-specific clocks
        const UInt clockCount = y.size()-(estimateKinematicPosition ? 3 : 0);
        updateClockError(listEpoch.at(i), mean(y.row(0, clockCount))/LIGHT_VELOCITY);

        clockDiff(listEpoch.at(i)) = mean(y.row(0, clockCount));
        if(estimateKinematicPosition)
          posDiff(listEpoch.at(i)) = norm(y.row(clockCount, 3));
      }

      if(maxabs(clockDiff)<=maxPosDiff)
        break;
    } // for(iter)

    for(UInt idEpoch=0; idEpoch<times.size(); idEpoch++)
    {
      const UInt idxPos = (estimateKinematicPosition ? idEpoch : 0);
      if((std::fabs(clockDiff(idEpoch)) > maxPosDiff) || (posDiff(idxPos) > maxPosDiff))
      {
// logWarning<<"estimateInitialClockError: "<<name()<<" disabled at "<<times.at(idEpoch).dateTimeStr()<<" due to differences: pos="<<posDiff(idxPos)<<" m, clock="<<clockDiff(idEpoch)<<" m"<<Log::endl;
        disable(idEpoch);
      }
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/

Gnss::Receiver::ObservationEquationList::ObservationEquationList() {}

/***********************************************/

Gnss::Receiver::ObservationEquationList::ObservationEquationList(const Gnss &gnss, const Receiver &receiver, AnalysisType analysisType)
{
  try
  {
    eqn.resize(receiver.idEpochSize());
    for(UInt idEpoch=0; idEpoch<eqn.size(); idEpoch++)
    {
      eqn.at(idEpoch).resize(receiver.idTransmitterSize(idEpoch));
      for(UInt idTrans=0; idTrans<eqn.at(idEpoch).size(); idTrans++)
      {
        Observation *obs = receiver.observation(idTrans, idEpoch);
        std::vector<GnssType> types;
        if(obs && obs->observationList(analysisType, types))
          eqn.at(idEpoch).at(idTrans) = std::unique_ptr<ObservationEquation>(new ObservationEquation(*obs, receiver, *gnss.transmitter.at(idTrans), gnss.ionosphere, idEpoch, FALSE, types));
      }
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Gnss::ObservationEquation *Gnss::Receiver::ObservationEquationList::operator()(UInt idTrans, UInt idEpoch) const
{
  if((idEpoch < eqn.size()) && (idTrans < eqn.at(idEpoch).size()))
    return eqn[idEpoch][idTrans].get();
  return nullptr;
}

/***********************************************/

void Gnss::Receiver::ObservationEquationList::deleteObservationEquation(UInt idTrans, UInt idEpoch)
{
  if((*this)(idTrans, idEpoch))
    eqn[idEpoch][idTrans] = nullptr;
}

/***********************************************/
/***********************************************/

void Gnss::Receiver::disableEpochsWithGrossCodeObservationOutliers(ObservationEquationList &eqnList, Double threshold)
{
  try
  {
    for(UInt idEpoch=0; idEpoch<idEpochSize(); idEpoch++)
      if(useable(idEpoch))
      {
        // delete all observations to a satellite at epoch if they contain a gross code outlier
        UInt outlierCount = 0;
        UInt totalCount   = 0;
        for(UInt idTrans=0; idTrans<idTransmitterSize(idEpoch); idTrans++)
          if(eqnList(idTrans, idEpoch))
          {
            const ObservationEquation &eqn = *eqnList(idTrans, idEpoch);
            for(UInt idType=0; idType<eqn.types.size(); idType++)
              if((eqn.types.at(idType) == GnssType::RANGE) && (std::fabs(eqn.l.at(idType)) >= threshold))
              {
                deleteObservation(idTrans, idEpoch);
                eqnList.deleteObservationEquation(idTrans, idEpoch);
                outlierCount++;
                break;
              }
            totalCount++;
          }

        // disable receiver at epoch if 50% or more of the observed satellites have gross code outliers
        const Double outlierRatio = static_cast<Double>(outlierCount)/totalCount;
        if(outlierRatio >= 0.5)
        {
//          logWarning << "grossCodeOutlier: "<< name() << " disabled at " << times.at(idEpoch).dateTimeStr() << " due to "
//                     << (outlierRatio*100)%"%.1f"s << "% of observed satellites having gross code outliers" << Log::endl;
          for(UInt idTrans=0; idTrans<idTransmitterSize(idEpoch); idTrans++)
            eqnList.deleteObservationEquation(idTrans, idEpoch);
          disable(idEpoch);
        }
      }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void Gnss::Receiver::createTracks(UInt minObsCountPerTrack, const std::vector<GnssType> &extraTypes)
{
  try
  {
    // delete old tracks
    // -----------------
    track.clear();

    for(UInt idTrans=0; idTrans<gnss().transmitter.size(); idTrans++)
    {
      UInt idEpochStart = 0;
      for(;;)
      {
        // find continuous track
        // ---------------------
        // find first epoch of the track
        while((idEpochStart < idEpochSize()) && !observation(idTrans, idEpochStart))
          idEpochStart++;

        Observation *obs = observation(idTrans, idEpochStart);
        if(!obs) // at end?
          break;

        // phase types of start epoch
        std::vector<GnssType> types;
        for(UInt idType= 0; idType <obs->size(); idType++)
          if(obs->at(idType).type == GnssType::PHASE)
            types.push_back(obs->at(idType).type);
        std::sort(types.begin(), types.end());

        // find last epoch of the track
        UInt countEpoch = 1;
        UInt idEpochEnd = idEpochStart;
        for(UInt idEpoch=idEpochStart+1; idEpoch<idEpochSize(); idEpoch++)
        {
          Observation *obs = observation(idTrans, idEpoch);
          if(obs)
          {
            if(times.at(idEpoch) > times.at(idEpochEnd)+seconds2time(1.5*observationSampling))
              break;

            // test types
            std::vector<GnssType> typesNew;
            for(UInt idType=0; idType<obs->size(); idType++)
              if(obs->at(idType).type == GnssType::PHASE)
                typesNew.push_back(obs->at(idType).type);
            if(!GnssType::allEqual(types, typesNew))
              break;

            idEpochEnd = idEpoch;
            countEpoch++;
          }
        } // for(idEpoch)

        // need phases at two frequencies (additional to extraTypes)
        std::set<GnssType> typeFrequencies;
        for(GnssType type : types)
          if(GnssType::index(extraTypes, type) == NULLINDEX)
            typeFrequencies.insert(type & GnssType::FREQUENCY);

        // define track
        if((countEpoch >= minObsCountPerTrack) && (typeFrequencies.size() >= 2))
        {
          track.push_back(std::make_shared<Track>(gnss().receiver.at(idRecv()).get(), gnss().transmitter.at(idTrans).get(), idEpochStart, idEpochEnd, types));
          for(UInt idEpoch=idEpochStart; idEpoch<=idEpochEnd; idEpoch++)
            if(observation(idTrans, idEpoch))
              observation(idTrans, idEpoch)->track = track.back().get();
        }
        else
        {
          for(UInt idEpoch=idEpochStart; idEpoch<=idEpochEnd; idEpoch++)
            deleteObservation(idTrans, idEpoch);
        }

        idEpochStart = idEpochEnd + 1;
      } // for(;;)
    } // for(idTrans)
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void Gnss::Receiver::writeTracks(const FileName &fileName, VariableList varList, const ObservationEquationList &eqnList)
{
  try
  {
    addVariable("prn",            varList);
    addVariable("trackTimeStart", varList);
    addVariable("trackTimeEnd",   varList);

    for(const auto &t : track)
    {
      // determine Melbourne-Wuebbena-like linear combinations
      std::vector<UInt>     idEpochs;
      std::vector<Vector>   estimates;
      std::vector<GnssType> typesPhase;
      Matrix Bias;
      t->combinations(eqnList, idEpochs, estimates, typesPhase, Bias);
      if(typesPhase.size() < 2)
        continue;

      // determine range & TEC based on phase observations only
      std::vector<Double> range(idEpochs.size());
      std::vector<Double> tec(idEpochs.size());
      t->rangeAndTecFromPhase(eqnList, typesPhase, idEpochs, range, tec);

      // determine cycle slip size in terms of TEC
      Vector l = Bias.column(0); // cycle slips can only occur in this linear combination anymore
      Matrix A(typesPhase.size(), 2);
      for(UInt idType=0; idType<typesPhase.size(); idType++)
      {
        A(idType, 0) = 1; // range
        A(idType, 1) = -Ionosphere::Ap/std::pow(typesPhase.at(idType).frequency(), 2); // TEC
      }
      const Double cycles2tecu = std::fabs(leastSquares(A, l)(1,0)); // one cycle slip in terms of TEC

      Matrix out(idEpochs.size(), 2+Bias.columns());
      for(UInt i=0; i<idEpochs.size(); i++)
      {
        out(i,0) = times.at(idEpochs.at(i)).mjd();
        out(i,1) = tec.at(i)/cycles2tecu;
        copy(estimates.at(i).slice(2, Bias.columns()).trans(), out.slice(i, 2, 1, Bias.columns()));
      }
      varList["prn"]->setValue(t->transmitter->name());
      varList["trackTimeStart"]->setValue(times.at(t->idEpochStart).mjd());
      varList["trackTimeEnd"]->setValue(times.at(t->idEpochEnd).mjd());
      writeFileMatrix(fileName(varList), out);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void Gnss::Receiver::deleteTrack(ObservationEquationList &eqnList, UInt idTrack)
{
  try
  {
    for(UInt idEpoch=track.at(idTrack)->idEpochStart; idEpoch<=track.at(idTrack)->idEpochEnd; idEpoch++)
    {
      deleteObservation(track.at(idTrack)->transmitter->idTrans(), idEpoch);
      eqnList.deleteObservationEquation(track.at(idTrack)->transmitter->idTrans(), idEpoch);
    }
    track.erase(track.begin()+idTrack);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Gnss::TrackPtr Gnss::Receiver::splitTrack(ObservationEquationList &eqnList, Gnss::TrackPtr track, UInt idEpochSplit)
{
  try
  {
    // new track
    const UInt idTrans = track->transmitter->idTrans();
    TrackPtr trackNew = std::make_shared<Track>(track->receiver, track->transmitter, idEpochSplit, track->idEpochEnd, track->types);
    this->track.push_back(trackNew);

    // shorten old track
    track->idEpochEnd = idEpochSplit-1;

    // connect observations to new track
    for(UInt idEpoch=trackNew->idEpochStart; idEpoch<=trackNew->idEpochEnd; idEpoch++)
      if(observation(idTrans, idEpoch))
        observation(idTrans, idEpoch)->track = trackNew.get();

    // connect observation equations to new track
    for(UInt idEpoch=trackNew->idEpochStart; idEpoch<=trackNew->idEpochEnd; idEpoch++)
      if(eqnList(idTrans, idEpoch))
        eqnList(idTrans, idEpoch)->track = trackNew.get();

    return trackNew;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void Gnss::Receiver::removeLowElevationTracks(ObservationEquationList &eqnList, Angle minElevation)
{
  try
  {
    for(UInt idTrack=track.size(); idTrack-->0;)
    {
      const UInt idTrans = track.at(idTrack)->transmitter->idTrans();
      Bool removeTrack = TRUE;
      for(UInt idEpoch=track.at(idTrack)->idEpochStart; idEpoch<=track.at(idTrack)->idEpochEnd; idEpoch++)
        if(eqnList(idTrans, idEpoch) && (eqnList(idTrans, idEpoch)->elevationRecvAnt >= minElevation))
        {
          removeTrack = FALSE;
          break;
        }

      if(removeTrack)
        deleteTrack(eqnList, idTrack);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Bool Gnss::Receiver::removeTracksWithInsufficientObservationCount(ObservationEquationList &eqnList, UInt minObsCountPerTrack, Gnss::AnalysisType analysisType)
{
  try
  {
    Bool removed = FALSE;
    Bool epochRemoved = FALSE;
    do
    {
      for(UInt idTrack=track.size(); idTrack-->0;)
        if(track.at(idTrack)->countObservations() < minObsCountPerTrack)
          deleteTrack(eqnList, idTrack);

      epochRemoved = FALSE;
      for(UInt idEpoch=0; idEpoch<idEpochSize(); idEpoch++)
        if(!isEpochEstimable(analysisType, idEpoch))
        {
          for(UInt idTrans=0; idTrans<idTransmitterSize(idEpoch); idTrans++)
            if(observation(idTrans, idEpoch))
            {
              eqnList.deleteObservationEquation(idTrans, idEpoch);
              epochRemoved = TRUE;
              removed = TRUE;
            }
          disable(idEpoch);
        }
    }
    while(epochRemoved);

    return removed;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void Gnss::Receiver::trackOutlierDetection(const ObservationEquationList &eqnList, const std::vector<GnssType> &ignoreTypes)
{
  try
  {
    for(auto &track : this->track)
    {
      const UInt idTrans      = track->transmitter->idTrans();
      const UInt idEpochStart = track->idEpochStart;
      const UInt idEpochEnd   = track->idEpochEnd;

      // available observations for this track
      std::vector<GnssType> types;
      for(UInt idEpoch=idEpochStart; idEpoch<=idEpochEnd; idEpoch++)
        if(eqnList(idTrans, idEpoch))
          for(const GnssType &type : eqnList(idTrans, idEpoch)->types)
            if((GnssType::index(types, type) == NULLINDEX) && (GnssType::index(ignoreTypes, type) == NULLINDEX))
              types.push_back(type);
      std::sort(types.begin(), types.end());

      // determine biases (reduced by range and TEC)
      Matrix Bias = identityMatrix(types.size());
      Matrix B(types.size(), 2);
      for(UInt idType=0; idType<types.size(); idType++)
        if(types.at(idType) == GnssType::RANGE)
        {
          B(idType, 0) = 1.; //range
          B(idType, 1) = Ionosphere::Ap/std::pow(types.at(idType).frequency(), 2); // TEC
        }
      const Vector tau = QR_decomposition(B);
      QMult(B, tau, Bias);
      Bias = Bias.column(B.columns(), Bias.rows()-B.columns());

      // setup observation equations: range, TEC, ambiguities
      // ----------------------------------------------------
      std::vector<Matrix> listl, listA;
      std::vector<Matrix> listlFull, listAFull, listBFull;
      std::vector<UInt>   listEpoch, listObsCount;

      for(UInt idEpoch=idEpochStart; idEpoch<=idEpochEnd; idEpoch++)
        if(eqnList(idTrans, idEpoch))
        {
          const ObservationEquation &eqn = *eqnList(idTrans, idEpoch);

          // observations
          Vector l = eqn.l;

          // distance and TEC
          Matrix B(l.rows(), 2);
          copy(eqn.A.column(ObservationEquation::idxRange), B.column(0));
          copy(eqn.B, B.column(1));

          // signal biases (includes ambiguities)
          Matrix A(l.rows(), Bias.columns());
          for(UInt idType=0; idType<eqn.types.size(); idType++)
          {
            const UInt idx = GnssType::index(types, eqn.types.at(idType));
            if(idx != NULLINDEX)
              matMult(1., eqn.A.column(ObservationEquation::idxUnit+idType), Bias.row(idx), A);
          }

          // decorrelate
          for(UInt i=0; i<l.rows(); i++)
          {
            l(i)     *= 1./eqn.sigma(i);
            A.row(i) *= 1./eqn.sigma(i);
            B.row(i) *= 1./eqn.sigma(i);
          }

          // remove ignored types
          UInt obsCount = l.rows();
          for(UInt idType=0; idType<eqn.types.size(); idType++)
            if(GnssType::index(ignoreTypes, eqn.types.at(idType)) != NULLINDEX)
            {
              l.row(idType).setNull();
              A.row(idType).setNull();
              B.row(idType).setNull();
              obsCount--;
            }

          listEpoch.push_back(idEpoch);
          listlFull.push_back(l);
          listAFull.push_back(A);
          listBFull.push_back(B);

          eliminationParameter(B, A, l);
          listl.push_back(l);
          listA.push_back(A);
          listObsCount.push_back(obsCount-B.columns());
        } // for(idEpoch)

      // estimate solution
      // -----------------
      Vector sigma0(listl.size(), 1.);
      Vector sigma = sigma0;
      Vector x = robustLeastSquares(huber, huberPower, listA, listl, listObsCount, sigma0, sigma);

      // downweight outliers
      // -------------------
      for(UInt i=0; i<listEpoch.size(); i++)
      {
        const ObservationEquation &eqn = *eqnList(idTrans, listEpoch.at(i));
        Observation &obs = *observation(idTrans, listEpoch.at(i));
        for(UInt idType=0; idType<eqn.types.size(); idType++)
          obs.at(eqn.types.at(idType)).sigma *= sigma(i);
      }

      // reduce integer ambiguities
      // --------------------------
      Vector b = Bias * x;
      for(UInt i=0; i<listEpoch.size(); i++)
        for(UInt idType=0; idType<types.size(); idType++)
          if(types.at(idType) == GnssType::PHASE)
          {
            const Double lambda = LIGHT_VELOCITY/types.at(idType).frequency();
            b(idType) = lambda * std::round(b(idType)/lambda);
            eqnList(idTrans, listEpoch.at(i))->l(GnssType::index(eqnList(idTrans, listEpoch.at(i))->types, types.at(idType))) -= b(idType);
            observation(idTrans, listEpoch.at(i))->at(types.at(idType)).observation -= b(idType);
          }
    } // for(track)
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Matrix Gnss::Receiver::robustLeastSquares(Double huber, Double huberPower,
                                          const_MatrixSliceRef A, const_MatrixSliceRef l,
                                          const Vector &sigma0, Vector &sigma)
{
  try
  {
    Matrix x;
    UInt   countOutlierOld = 0;
    Double sigmaOld = 1.;
    for(UInt iter=0; iter<30; iter++)
    {
      // weighting
      Matrix Wl = l;
      Matrix WA = A;
      for(UInt i=0; i<sigma.rows(); i++)
      {
        Wl.row(i) *= 1./sigma(i);
        WA.row(i) *= 1./sigma(i);
      }

      // QR decomposition
      Vector tau = QR_decomposition(WA);
      QTransMult(WA, tau, Wl); // transform observations: l:= Q'l
      x = Wl.row(0, WA.columns());
      triangularSolve(1., WA.row(0,WA.columns()), x);
      Wl.row(0, WA.columns()).setNull(); // residuals: remove WB*x
      QMult(WA, tau, Wl); // back transformation
      generateQ(WA, tau); // for redundancies

      // outlier detection
      UInt   countOutlier = 0;
      Double ePeSum       = 0.;
      Double rSum         = 0.;
      for(UInt i=0; i<Wl.rows(); i++)
      {
        const Double e = sigma(i) * Wl(i,0);
        const Double r = 1 - quadsum(WA.row(i));
        const Double s = sqrt(e*e/r)/sigma0(i);
        if((s > huber) && (r > 1e-4)) // redundancy: it is possible to estimate sigma?
        {
          sigma(i) = std::pow(s/huber,huberPower) * sigma0(i);
          countOutlier++;
        }
        else
        {
          sigma(i) = sigma0(i);
          ePeSum += std::pow(e/sigma0(i), 2);
          rSum   += r;
        }
      }

      const Double sigma = sqrt(ePeSum/rSum);
      if((countOutlier==0) || ((countOutlier==countOutlierOld) && (fabs(sigma-sigmaOld)<0.001)))
        break;
      sigmaOld        = sigma;
      countOutlierOld = countOutlier;
    } // for(iter)

    return x;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Matrix Gnss::Receiver::robustLeastSquares(Double huber, Double huberPower,
                                          const std::vector<Matrix> &A, const std::vector<Matrix> &l, const std::vector<UInt> &obsCount,
                                          const Vector &sigma0, Vector &sigma)
{
  try
  {
    const UInt countEpoch = l.size();
    if(!countEpoch)
      return Matrix();
    UInt countObsTotal = 0;
    for(UInt idEpoch=0; idEpoch<countEpoch; idEpoch++)
      countObsTotal += l.at(idEpoch).rows();

    Matrix x;
    UInt   countOutlierOld = 0;
    Double sigmaOld = 1.;
    const UInt maxIter = 30;
    for(UInt iter=0; iter<maxIter; iter++)
    {
      // Sort into one system
      // --------------------
      Vector Wl(countObsTotal);
      Matrix WA(countObsTotal, A.at(0).columns());
      UInt index = 0;
      for(UInt idEpoch=0; idEpoch<countEpoch; idEpoch++)
      {
        axpy(1/sigma(idEpoch), l.at(idEpoch), Wl.row(index, l.at(idEpoch).rows()));
        axpy(1/sigma(idEpoch), A.at(idEpoch), WA.row(index, A.at(idEpoch).rows()));
        index += l.at(idEpoch).rows();
      }

      // solve system
      // ------------
      const Vector tau = QR_decomposition(WA);
      QTransMult(WA, tau, Wl);        // transform observations: l:= Q'l
      x = Wl.row(0, WA.columns());
      triangularSolve(1., WA.row(0,tau.rows()), x);
      Wl.row(0, tau.rows()).setNull();   // residuals: remove A*x
      QMult(WA, tau, Wl);                // back transformation
      generateQ(WA, tau);                // for redundancy computation

      UInt countOutlier = 0;
      Double ePeSum = 0;
      Double rSum   = 0;
      index = 0;
      for(UInt idEpoch=0; idEpoch<countEpoch; idEpoch++)
      {
        const Double ePe      = quadsum(Wl.row(index, l.at(idEpoch).rows()));
        const Double r        = obsCount.at(idEpoch)-quadsum(WA.row(index, l.at(idEpoch).rows())); // redundancies
        const Double s        = sqrt(ePe/r) * sigma(idEpoch)/sigma0(idEpoch);

        sigma(idEpoch) = sigma0(idEpoch);
        if((s > huber) && (r > 1e-4))
        {
          sigma(idEpoch) *= std::pow(s/huber, huberPower);
          countOutlier   += obsCount.at(idEpoch);
        }
        else
        {
          ePeSum += ePe;
          rSum   += r;
        }

        index += l.at(idEpoch).rows();
      } // for(idEpoch)

      const Double sigma = sqrt(ePeSum/rSum);
      if((countOutlier==0) || ((countOutlier==countOutlierOld) && (fabs(sigma-sigmaOld)<0.001)))
        break;
      sigmaOld = sigma;
      countOutlierOld = countOutlier;

      if(iter>=maxIter-1)
        logWarning<<"Gnss::Receiver::robustLeastSquares: no convergence, sigma="<<sigma<<", outlier="<<countOutlier<<" of "<<countObsTotal<<Log::endl;
    } // for(iter)

    return x;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/

void Gnss::Receiver::cycleSlipsDetection(ObservationEquationList &eqnList, Double lambda, UInt windowSize, Double tecSigmaFactor)
{
  try
  {
    for(UInt idTrack=0; idTrack<track.size(); idTrack++)
      cycleSlipsDetection(eqnList, track.at(idTrack), lambda, windowSize, tecSigmaFactor);

    // sort new created tracks: transmitter by transmitter, secondly by epochStart
    std::sort(track.begin(), track.end(), [](auto &t1, auto &t2)
              {return (t1->transmitter->idTrans() == t2->transmitter->idTrans()) ? (t1->idEpochStart < t2->idEpochStart) : (t1->transmitter->idTrans() < t2->transmitter->idTrans());});
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void Gnss::Receiver::cycleSlipsDetection(ObservationEquationList &eqnList, TrackPtr track, Double lambda, UInt windowSize, Double tecSigmaFactor)
{
  try
  {
    const UInt idTrans = track->transmitter->idTrans();

    // determine Melbourne-Wuebbena-like linear combinations
    // -----------------------------------------------------
    std::vector<UInt>     idEpochs;
    std::vector<Vector>   estimates;
    std::vector<GnssType> typesPhase;
    Matrix Bias;
    track->combinations(eqnList, idEpochs, estimates, typesPhase, Bias);

    // loop over all Melbourne-Wuebbena-like linear combinations
    // starting with the most accurate, skip the last, inaccurate one
    // --------------------------------------------------------------
    for(UInt k=2+Bias.columns(); k-->3;) // k=0: range, k=1: TEC
    {
      std::vector<Double> value(idEpochs.size());
      for(UInt i=0; i<value.size(); i++)
        value.at(i) = estimates.at(i)(k);
      if(k < 4) // only denoise inaccurate linear combinations
        value = totalVariationDenoising(value, lambda);
      // split track if denoised difference exceeds 3/4 cycle
      Bool split = FALSE;
      for(UInt i=value.size(); i-->1;)
        if(std::fabs(value.at(i)-value.at(i-1)) > 0.75)
        {
          splitTrack(eqnList, track, idEpochs.at(i));
          split = TRUE;
        }
      if(split)
      {
        cycleSlipsDetection(eqnList, track, lambda, windowSize, tecSigmaFactor); // reprocess (shorter) track again
        return;
      }
    }

    const GnssType typeL5 = GnssType::PHASE + GnssType::L5 + GnssType::GPS;
    const Bool isGpsL5 = (GnssType::index(track->types, typeL5) != NULLINDEX);

    if(windowSize || isGpsL5)
    {
      // determine range & TEC based on phase observations only
      // ------------------------------------------------------
      std::vector<Double> range(idEpochs.size());
      std::vector<Double> tec(idEpochs.size());
      track->rangeAndTecFromPhase(eqnList, typesPhase, idEpochs, range, tec);

      // split track at detected cycle slips based on smoothness analysis of TEC
      // -----------------------------------------------------------------------
      if(windowSize)
      {
        // determine cycle slip size in terms of TEC
        // -----------------------------------------
        Vector l = Bias.column(0); // cycle slips can only occur in this linear combination anymore
        Matrix A(typesPhase.size(), 2);
        for(UInt idType=0; idType<typesPhase.size(); idType++)
        {
          A(idType, 0) = 1; // range
          A(idType, 1) = -Ionosphere::Ap/std::pow(typesPhase.at(idType).frequency(), 2); // TEC
        }
        const Double cycles2tecu = std::fabs(leastSquares(A, l)(1,0)); // one cycle slip in terms of TEC

        // find cycle slips in TEC based on moving window over autoregressive model residuals
        // ----------------------------------------------------------------------------------
        const UInt order = 3; // AR model order
        std::vector<UInt> slips;
        if(tec.size() >= order+windowSize)
        {
          // high pass filter via AR model
          const Vector data(tec);
          Vector l = data.row(order, data.size()-order);
          Matrix A = Matrix(l.rows(), order);
          for(UInt k = 0; k<order; k++)
            copy(data.row(order-k-1, data.rows()-order), A.column(k));
          leastSquares(A, l); // l contains AR model residuals after function call

          // peak/outlier detection using moving standard deviation over AR model residuals.
          // automatic threshold scaling via standard deviation is used to prevent excessive
          // splitting during periods with high ionospheric variations/scintillations
          for(UInt idxStart = 0; idxStart < l.size()-windowSize; idxStart++)
            if(l(idxStart+windowSize/2) > std::max(0.9*cycles2tecu, tecSigmaFactor*standardDeviation(l.row(idxStart, windowSize))))
              slips.push_back(idEpochs.at(idxStart+windowSize/2+order));
        }

        Bool split = FALSE;
        for(UInt i=slips.size(); i-->0;)
        {
          splitTrack(eqnList, track, slips.at(i));
          split = TRUE;
        }
        if(split && isGpsL5)
        {
          cycleSlipsDetection(eqnList, track, lambda, windowSize, tecSigmaFactor); // reprocess (shorter) track again
          return;
        }
      }

      // repair GPS L5 cycle slips
      // -------------------------
      if(isGpsL5)
      {
        const Double wavelength = wavelengthFactor() * LIGHT_VELOCITY/typeL5.frequency();
        const Double TEC = (-Ionosphere::Ap/std::pow(typeL5.frequency(), 2));
        for(UInt idType=0; idType<track->types.size(); idType++)
          if(track->types.at(idType) == typeL5)
          {
            // reduce l5 by estimated range and tec
            std::vector<Double> l5(idEpochs.size());
            for(UInt i=0; i<idEpochs.size(); i++)
            {
              const ObservationEquation &eqn = *eqnList(idTrans, idEpochs.at(i));
              l5.at(i) = eqn.l(GnssType::index(eqn.types, track->types.at(idType))) - range.at(i) - TEC * tec.at(i);
            }

            // fix jumps
            Double jump = 0;
            for(UInt i=1; i<l5.size(); i++)
            {
              jump += wavelength * std::round((l5.at(i)-l5.at(i-1))/wavelength);
              eqnList(idTrans, idEpochs.at(i))->l(GnssType::index(eqnList(idTrans, idEpochs.at(i))->types, track->types.at(idType))) -= jump;
              observation(idTrans, idEpochs.at(i))->at(track->types.at(idType)).observation -= jump;
            }
          }
      } // if(isGpsL5)
    } // if(windowSize || isGpsL5)
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void Gnss::Receiver::cycleSlipsRepairAtSameFrequency(ObservationEquationList &eqnList)
{
  try
  {
    // get all phase types
    std::vector<GnssType> types;
    for(UInt idEpoch=0; idEpoch<idEpochSize(); idEpoch++)
      for(UInt idTrans=0; idTrans<idTransmitterSize(idEpoch); idTrans++)
        if(eqnList(idTrans, idEpoch))
        {
          const ObservationEquation &eqn = *eqnList(idTrans, idEpoch);
          for(UInt idType=0; idType<eqn.types.size(); idType++)
            if((eqn.types.at(idType) == GnssType::PHASE) && (GnssType::index(types, eqn.types.at(idType)) == NULLINDEX))
              types.push_back(eqn.types.at(idType) & ~GnssType::PRN);
        }

    // find two phase observations with same system and frequency
    for(UInt idType1=0; idType1<types.size(); idType1++)
      for(UInt idType2=idType1+1; idType2<types.size(); idType2++)
        if(types.at(idType1) == (types.at(idType2) & (GnssType::FREQUENCY + GnssType::SYSTEM + GnssType::FREQ_NO)))
        {
          const Double wavelength = wavelengthFactor() * LIGHT_VELOCITY/types.at(idType1).frequency();

          // compute difference
          std::vector<UInt>   idxTrans, idxEpoch;
          std::vector<Double> value;
          for(UInt idEpoch=0; idEpoch<idEpochSize(); idEpoch++)
            for(UInt idTrans=0; idTrans<idTransmitterSize(idEpoch); idTrans++)
              if(eqnList(idTrans, idEpoch))
              {
                const ObservationEquation &eqn = *eqnList(idTrans, idEpoch);
                const UInt idx1 = GnssType::index(eqn.types, types.at(idType1));
                const UInt idx2 = GnssType::index(eqn.types, types.at(idType2));
                if((idx1 == NULLINDEX) || (idx2 == NULLINDEX))
                  continue;
                idxTrans.push_back(idTrans);
                idxEpoch.push_back(idEpoch);
                value.push_back((eqn.l(idx1) - eqn.l(idx2))/wavelength);
              }

          if(!value.size())
            throw(Exception(name()+": "+types.at(idType1).str()+" and "+types.at(idType2).str()+" must be observed together"));

          // fix jumps
          const Double v0 = value.at(0)-std::round(value.at(0)); // consider bias (e.g. quarter cycles)
          for(UInt i=0; i<value.size(); i++)
          {
            const Double v = wavelength * std::round(value.at(i)-v0);
            eqnList(idxTrans.at(i), idxEpoch.at(i))->l(GnssType::index(eqnList(idxTrans.at(i), idxEpoch.at(i))->types, types.at(idType2))) += v;
            observation(idxTrans.at(i), idxEpoch.at(i))->at(types.at(idType2)).observation += v;
          }
        } // for(idType1, idType2)
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void Gnss::Receiver::addRotiPseudoObservations(ObservationEquationList &eqnList, Double windowSize)
{
  try
  {
    for(UInt idTrack=0; idTrack<track.size(); idTrack++)
    {
      const UInt     idTrans = track.at(idTrack)->transmitter->idTrans();
      const GnssType prn     = track.at(idTrack)->transmitter->PRN();

      // get valid epochs
      std::vector<UInt> idEpochs;
      for(UInt idEpoch=track.at(idTrack)->idEpochStart; idEpoch<=track.at(idTrack)->idEpochEnd; idEpoch++)
        if(eqnList(idTrans, idEpoch) && observation(idTrans, idEpoch))
           idEpochs.push_back(idEpoch);

      // estimate TEC
      std::vector<Double> range, STEC;
      track.at(idTrack)->rangeAndTecFromPhase(eqnList, track.at(idTrack)->types, idEpochs, range, STEC);
      if(STEC.size() < 2)
        continue;

      // compute TEC change as single difference
      std::vector<Double> dTEC(idEpochs.size());
      dTEC.at(0) = (STEC.at(1)-STEC.at(0))/(times.at(idEpochs.at(1))-times.at(idEpochs.at(0))).seconds();
      for(UInt i=1; i<STEC.size(); i++)
        dTEC.at(i) = (STEC.at(i)-STEC.at(i-1))/(times.at(idEpochs.at(i))-times.at(idEpochs.at(i-1))).seconds();

      // compute ROTI
      std::vector<Double> ROTI(idEpochs.size(), NAN_EXPR);
      UInt idxStart = 0;
      for(UInt i=0; i<ROTI.size(); i++)
      {
        while((idxStart<dTEC.size()) && ((times.at(idEpochs.at(idxStart))-times.at(idEpochs.at(i))).seconds() < -windowSize/2))
          idxStart++;
        UInt   count   = 0;
        Double sum     = 0;
        Double quadsum = 0;
        while((idxStart+count<dTEC.size()) && ((times.at(idEpochs.at(idxStart+count))-times.at(idEpochs.at(i))).seconds() < windowSize/2))
        {
          sum     += dTEC.at(idxStart+count);
          quadsum += std::pow(dTEC.at(idxStart+count), 2);
          count++;
        }
        if(count)
          ROTI.at(i) = std::sqrt(quadsum/count - std::pow(sum/count, 2));
      }

      // add (pseudo) observation
      for(UInt i=0; i<ROTI.size(); i++)
        observation(idTrans, idEpochs.at(i))->push_back(SingleObservation(GnssType::ROTI + prn, ROTI.at(i)));
    } // for(idTrack)
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

// Total variation denoising algorithm source:
// Laurent Condat. A Direct Algorithm for 1D Total Variation Denoising. IEEE Signal Processing Letters,
// Institute of Electrical and Electronics Engineers, 2013, 20 (11), pp.1054-1057. DOI: 10.1109/LSP.2013.2278339.
std::vector<Double> Gnss::Receiver::totalVariationDenoising(const std::vector<Double> &y, Double lambda)
{
  try
  {
    // initialize total variation denoising algorithm
    std::vector<Double> x(y.size(), 0.);  // output
    UInt N = y.size()-1;
    UInt k  = 0;                    // current sample location
    UInt k0 = 0;                    // beginning of current segment
    UInt km = 0;                    // last position where umax = -lambda
    UInt kp = 0;                    // last position where umin =  lambda
    Double vMin = y.at(0) - lambda; // lower bound for the segment's value
    Double vMax = y.at(0) + lambda; // upper bound for the segment's value
    Double uMin =  lambda;          // u is the dual variable
    Double uMax = -lambda;          // u is the dual variable

    // total variation denoising algorithm
    for(;;)
    {
      if(k == N)
      {
        x.at(N) = vMin + uMin;
        break;
      }

      if(y.at(k+1) + uMin < vMin - lambda)      // negative jump necessary
      {
        for(UInt i = k0; i <= km; i++)
          x.at(i) = vMin;
        k = k0 = km = kp = km+1;
        vMin = y.at(k);
        vMax = y.at(k) + 2*lambda;
        uMin =  lambda;
        uMax = -lambda;
      }
      else if(y.at(k+1) + uMax > vMax + lambda) // positive jump necessary
      {
        for(UInt i = k0; i <= kp; i++)
          x.at(i) = vMax;
        k = k0 = km = kp = kp+1;
        vMin = y.at(k) - 2*lambda;
        vMax = y.at(k);
        uMin =  lambda;
        uMax = -lambda;
      }
      else  // no jump necessary
      {
        k = k+1;
        uMin = uMin + y.at(k) - vMin;
        uMax = uMax + y.at(k) - vMax;
        if(uMin >= lambda)  // update of vMin
        {
          vMin = vMin + (uMin-lambda)/(k-k0+1);
          uMin = lambda;
          km = k;
        }
        if(uMax <= -lambda) // update of vMax
        {
          vMax = vMax + (uMax+lambda)/(k-k0+1);
          uMax = -lambda;
          kp = k;
        }
      }

      if(k < N)
        continue;

      if(uMin < 0.)       // vMin is too high ==> negative jump necessary
      {
        for(UInt i = k0; i <= km; i++)
          x.at(i) = vMin;
        k = k0 = km = km+1;
        vMin = y.at(k);
        uMin = lambda;
        uMax = y.at(k) + lambda - vMax;
        continue;
      }
      else if(uMax > 0.)  // vMax is too low ==> positive jump necessary
      {
        for(UInt i = k0; i <= kp; i++)
          x.at(i) = vMax;
        k = k0 = kp = kp+1;
        vMax = y.at(k);
        uMax = -lambda;
        uMin = y.at(k) - lambda - vMin;
        continue;
      }
      else
      {
        for(UInt i = k0; i <= N; i++)
          x.at(i) = vMin + uMin/(k-k0+1);
        break;
      }
    }

    return x;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/

UInt Gnss::Track::countObservations() const
{
  try
  {
    UInt count = 0;
    for(UInt idEpoch=idEpochStart; idEpoch<=idEpochEnd; idEpoch++)
      if(transmitter->useable(idEpoch) && receiver->observation(transmitter->idTrans(), idEpoch))
        count++;
    return count;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void Gnss::Track::removeAmbiguitiesFromObservations(const std::vector<GnssType> &types, const std::vector<Double> &value)
{
  try
  {
    for(UInt idEpoch=idEpochStart; idEpoch<=idEpochEnd; idEpoch++)
      if(receiver->observation(transmitter->idTrans(), idEpoch))
      {
        Observation &obs = *receiver->observation(transmitter->idTrans(), idEpoch);
        for(UInt idType=0; idType<obs.size(); idType++)
        {
          const UInt idx = GnssType::index(types, obs.at(idType).type);
          if(idx != NULLINDEX)
            obs.at(idType).observation -= value.at(idx);
        }
      }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void Gnss::Track::combinations(const Gnss::Receiver::ObservationEquationList &eqnList, std::vector<UInt> &idEpochs, std::vector<Vector> &estimates, std::vector<GnssType> &typesPhase, Matrix &Bias) const
{
  try
  {
    const UInt idTrans = transmitter->idTrans();

    // available observations for this track
    std::vector<GnssType> types;
    for(UInt idEpoch=idEpochStart; idEpoch<=idEpochEnd; idEpoch++)
      if(eqnList(idTrans, idEpoch) && receiver->observation(idTrans, idEpoch))
        for(const auto &type : eqnList(idTrans, idEpoch)->types)
          if(type != GnssType::L5 + GnssType::GPS && GnssType::index(types, type) == NULLINDEX) // ignore time variable GPS L5 signals
            types.push_back(type);
    std::sort(types.begin(), types.end());

    typesPhase = types;
    typesPhase.erase(std::remove_if(typesPhase.begin(), typesPhase.end(), [](const GnssType &x){return x != GnssType::PHASE;}), typesPhase.end());
    Bias = GnssParametrizationAmbiguities::phaseDecorrelation(typesPhase, receiver->wavelengthFactor());

    idEpochs.clear();
    estimates.clear();
    for(UInt idEpoch=idEpochStart; idEpoch<=idEpochEnd; idEpoch++)
      if(eqnList(idTrans, idEpoch) && receiver->observation(idTrans, idEpoch))
      {
        const ObservationEquation &eqn = *eqnList(idTrans, idEpoch);

// Code biases can generate apparent cycle slips
// if(!GnssType::allEqual(eqn.types, types))
//   logWarning<<"types not equal"<<Log::endl;

        // observations
        Vector l = eqn.l;

        // design matrix
        Matrix A(l.rows(), 2+Bias.columns());
        copy(eqn.A.column(ObservationEquation::idxRange), A.column(0)); // distance
        copy(eqn.B, A.column(1));                                       // TEC
        for(UInt idType=0; idType<eqn.types.size(); idType++)            // signal biases (includes ambiguities)
        {
          const UInt idx = GnssType::index(typesPhase, eqn.types.at(idType));
          if(idx != NULLINDEX)
            matMult(1., eqn.A.column(ObservationEquation::idxUnit+idType), Bias.row(idx), A.column(2, Bias.columns()));
        }

        // decorrelate
        for(UInt i=0; i<l.rows(); i++)
        {
          l(i)     *= 1./eqn.sigma0(i);
          A.row(i) *= 1./eqn.sigma0(i);
        }

        idEpochs.push_back(idEpoch);
        estimates.push_back(leastSquares(A, l));
      } // for(idEpoch)
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void Gnss::Track::rangeAndTecFromPhase(const Gnss::Receiver::ObservationEquationList &eqnList, const std::vector<GnssType> &typesPhase,
                                       const std::vector<UInt> &idEpochs, std::vector<Double> &range, std::vector<Double> &tec) const
{
  try
  {
    range.resize(idEpochs.size());
    tec.resize(idEpochs.size());
    for(UInt i=0; i<idEpochs.size(); i++)
    {
      const ObservationEquation &eqn = *eqnList(transmitter->idTrans(), idEpochs.at(i));

      Vector l(typesPhase.size());
      Matrix A(typesPhase.size(), 2);
      for(UInt idType=0; idType<eqn.types.size(); idType++)
      {
        const UInt idx = GnssType::index(typesPhase, eqn.types.at(idType));
        if(idx == NULLINDEX)
          continue;
        l(idx) = eqn.l(idType);
        A(idx, 0) = 1; // range
        A(idx, 1) = -Ionosphere::Ap/std::pow(eqn.types.at(idType).frequency(), 2); // TEC
      }

      const Vector x = leastSquares(A, l);
      range.at(i) = x(0);
      tec.at(i)   = x(1);
    } // for(idEpoch)
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
