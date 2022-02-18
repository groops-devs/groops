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

#include <random>
#include "base/import.h"
#include "base/string.h"
#include "parser/expressionParser.h"
#include "files/fileMatrix.h"
#include "files/fileInstrument.h"
#include "inputOutput/logging.h"
#include "misc/varianceComponentEstimation.h"
#include "gnss/gnssLambda.h"
#include "gnss/gnssObservation.h"
#include "gnss/gnssTransmitter.h"
#include "gnss/gnssReceiver.h"

/***********************************************/

GnssReceiver::GnssReceiver(Bool isMyRank, Bool isEarthFixed, const std::string &name, const GnssStationInfo &info,
                           GnssAntennaDefinition::NoPatternFoundAction noPatternFoundAction, const Vector &useableEpochs,
                           Bool integerAmbiguities, Double wavelengthFactor)
  : GnssTransceiver(name, info, noPatternFoundAction, useableEpochs), isMyRank_(isMyRank),
    isEarthFixed_(isEarthFixed), integerAmbiguities(integerAmbiguities), wavelengthFactor(wavelengthFactor)
{
}

/***********************************************/

// this seems to improve the performance
void GnssReceiver::copyObservations2ContinuousMemoryBlock()
{
  try
  {
    UInt count = 0;
    for(const auto &obsEpoch : observations_)
      for(const GnssObservation *obs : obsEpoch)
        if(obs)
          count++;
    obsMem.resize(count);
    obsMem.shrink_to_fit();
    // copy
    count = 0;
    for(auto &obsEpoch : observations_)
      for(GnssObservation *&obs : obsEpoch)
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

void GnssReceiver::disable(UInt idEpoch)
{
  try
  {
    GnssTransceiver::disable(idEpoch);
    if((idEpoch != NULLINDEX) && (idEpoch < observations_.size()))
      observations_.at(idEpoch).clear();
    if(!useable())
    {
      isMyRank_ = FALSE;
      obsMem.clear();
      obsMem.shrink_to_fit();
      observations_.clear();
      observations_.shrink_to_fit();
      tracks.clear();
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

GnssObservation *GnssReceiver::observation(UInt idTrans, UInt idEpoch) const
{
  if((idEpoch < observations_.size()) && (idTrans < observations_.at(idEpoch).size()))
    return observations_[idEpoch][idTrans];
  return nullptr;
}

/***********************************************/

void GnssReceiver::deleteObservation(UInt idTrans, UInt idEpoch)
{
  try
  {
    if(!observation(idTrans, idEpoch))
      return;
    observations_[idEpoch][idTrans] = nullptr;
    if(std::all_of(observations_[idEpoch].begin(), observations_[idEpoch].end(), [](auto obs) {return obs == nullptr;}))
      disable(idEpoch);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssReceiver::signalComposition(UInt /*idEpoch*/, const std::vector<GnssType> &types, std::vector<GnssType> &typesTrans, Matrix &A) const
{
  try
  {
    // composed type = factor1 * type1 + factor2 * type2
    static const std::vector<std::tuple<GnssType, GnssType, GnssType, Double, Double>> composites =
      {{GnssType::C1XG,  GnssType::C1SG, GnssType::C1LG, 0.5, 0.5},
       {GnssType::C2XG,  GnssType::C2SG, GnssType::C2LG, 0.5, 0.5},
       {GnssType::C5XG,  GnssType::C5IG, GnssType::C5QG, 0.5, 0.5},

       {GnssType::C4XR,  GnssType::C4AR, GnssType::C4BR, 0.5, 0.5},
       {GnssType::C6XR,  GnssType::C6AR, GnssType::C6BR, 0.5, 0.5},
       {GnssType::C3XR,  GnssType::C3IR, GnssType::C3QR, 0.5, 0.5},

       {GnssType::C1XE,  GnssType::C1BE, GnssType::C1CE, 0.5, 0.5},
       {GnssType::C5XE,  GnssType::C5IE, GnssType::C5QE, 0.5, 0.5},
       {GnssType::C7XE,  GnssType::C7IE, GnssType::C7QE, 0.5, 0.5},
       {GnssType::C8XE,  GnssType::C8IE, GnssType::C8QE, 0.5, 0.5},
       {GnssType::C6XE,  GnssType::C6BE, GnssType::C6CE, 0.5, 0.5},

       {GnssType::C2XC,  GnssType::C2IC, GnssType::C2QC, 0.5, 0.5},
       {GnssType::C1XC,  GnssType::C1DC, GnssType::C1PC, 0.5, 0.5},
       {GnssType::C1ZC,  GnssType::C1SC, GnssType::C1LC, 0.5, 0.5},
       {GnssType::C5XC,  GnssType::C5DC, GnssType::C5PC, 0.5, 0.5},
       {GnssType::C7XC,  GnssType::C7IC, GnssType::C7QC, 0.5, 0.5},
       {GnssType::C7ZC,  GnssType::C7DC, GnssType::C7PC, 0.5, 0.5},
       {GnssType::C8XC,  GnssType::C8DC, GnssType::C8PC, 0.5, 0.5},
       {GnssType::C6XC,  GnssType::C6IC, GnssType::C6QC, 0.5, 0.5},
       {GnssType::C6ZC,  GnssType::C6DC, GnssType::C6PC, 0.5, 0.5},

       {GnssType::C1XJ,  GnssType::C1SJ, GnssType::C1LJ, 0.5, 0.5},
       {GnssType::C2XJ,  GnssType::C2SJ, GnssType::C2LJ, 0.5, 0.5},
       {GnssType::C5XJ,  GnssType::C5IJ, GnssType::C5QJ, 0.5, 0.5},
       {GnssType::C5ZJ,  GnssType::C5DJ, GnssType::C5PJ, 0.5, 0.5},
       {GnssType::C6XJ,  GnssType::C6SJ, GnssType::C6LJ, 0.5, 0.5},
       {GnssType::C6ZJ,  GnssType::C6SJ, GnssType::C6EJ, 0.5, 0.5},

       // unknown attributes
       {GnssType::C2UG,  GnssType::C2SG, GnssType::C2LG, 0.5, 0.5},
       {GnssType::C5UG,  GnssType::C5IG, GnssType::C5QG, 0.5, 0.5},
       {GnssType::C1UE,  GnssType::C1BE, GnssType::C1CE, 0.5, 0.5},
       {GnssType::C5UE,  GnssType::C5IE, GnssType::C5QE, 0.5, 0.5},
       {GnssType::C7UE,  GnssType::C7IE, GnssType::C7QE, 0.5, 0.5},
       {GnssType::C8UE,  GnssType::C8IE, GnssType::C8QE, 0.5, 0.5},
       {GnssType::C6UE,  GnssType::C6BE, GnssType::C6CE, 0.5, 0.5}};

    typesTrans = GnssType::replaceCompositeSignals(types);

    A = Matrix(types.size(), typesTrans.size());
    for(UInt idType=0; idType<types.size(); idType++)
      if((types.at(idType) == GnssType::PHASE) || (types.at(idType) == GnssType::RANGE)) // only phase and code signals are transmitted (what about doppler?)
      {
        GnssType type = types.at(idType);

        UInt idx;
        if(type.isInList(typesTrans, idx))
        {
          A(idType, idx) = 1.; // signal observed directly
          continue;
        }

        if(type == GnssType::C2DG)
        {
          A(idType, GnssType::index(typesTrans, GnssType::C1CG)) = +1.;
          A(idType, GnssType::index(typesTrans, GnssType::C1WG)) = -1.;
          A(idType, GnssType::index(typesTrans, GnssType::C2WG)) = +1.;
          continue;
        }

        const auto composite = std::find_if(composites.begin(), composites.end(), [&](const auto &composite) {return (std::get<0>(composite) == type);});
        if(composite != composites.end())
        {
          A(idType, GnssType::index(typesTrans, std::get<1>(*composite))) = std::get<3>(*composite);
          A(idType, GnssType::index(typesTrans, std::get<2>(*composite))) = std::get<4>(*composite);
          continue;
        }

        throw(Exception("composite signal not implemented: "+type.str()));
      } // for(idType)
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/

GnssReceiver::ObservationEquationList::ObservationEquationList(const GnssReceiver &receiver, const std::vector<GnssTransmitterPtr> &transmitters,
                                                               const std::function<Rotary3d(const Time &time)> &rotationCrf2Trf,
                                                               const std::function<void(GnssObservationEquation &eqn)> &reduceModels,
                                                               GnssObservation::Group group)
{
  try
  {
    eqn.resize(receiver.idEpochSize());
    for(UInt idEpoch=0; idEpoch<eqn.size(); idEpoch++)
    {
      eqn.at(idEpoch).resize(receiver.idTransmitterSize(idEpoch));
      for(UInt idTrans=0; idTrans<eqn.at(idEpoch).size(); idTrans++)
      {
        GnssObservation *obs = receiver.observation(idTrans, idEpoch);
        std::vector<GnssType> types;
        if(obs && obs->observationList(group, types))
        {
          auto e = new GnssObservationEquation(*obs, receiver, *transmitters.at(idTrans), rotationCrf2Trf, reduceModels, idEpoch, FALSE, types);
          eqn.at(idEpoch).at(idTrans) = std::unique_ptr<GnssObservationEquation>(e);
        }
      }
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

GnssObservationEquation *GnssReceiver::ObservationEquationList::operator()(UInt idTrans, UInt idEpoch) const
{
  if((idEpoch < eqn.size()) && (idTrans < eqn.at(idEpoch).size()))
    return eqn[idEpoch][idTrans].get();
  return nullptr;
}

/***********************************************/

void GnssReceiver::ObservationEquationList::deleteObservationEquation(UInt idTrans, UInt idEpoch)
{
  if((*this)(idTrans, idEpoch))
    eqn[idEpoch][idTrans] = nullptr;
}

/***********************************************/
/***********************************************/

void GnssReceiver::readObservations(const FileName &fileName, const std::vector<GnssTransmitterPtr> &transmitters,
                                    std::function<Rotary3d(const Time &time)> rotationCrf2Trf, const Time &timeMargin, Angle elevationCutOff,
                                    const std::vector<GnssType> &useType, const std::vector<GnssType> &ignoreType, GnssObservation::Group group)
{
  try
  {
    GnssReceiverArc arc = InstrumentFile::read(fileName);

    std::vector<Time> observationTimes;
    Vector phaseWindup(transmitters.size());
    std::map<GnssType, UInt> removedTypes;

    UInt idEpoch = 0;
    for(UInt arcEpoch=0; arcEpoch<arc.size(); arcEpoch++)
    {
      // search time slot
      while((idEpoch < times.size()) && (times.at(idEpoch)+timeMargin < arc.at(arcEpoch).time))
        disable(idEpoch++);
      if(idEpoch >= times.size())
        break;
      if((arc.at(arcEpoch).time+timeMargin < times.at(idEpoch)) || !useable(idEpoch))
        continue;
      times.at(idEpoch) = arc.at(arcEpoch).time;
//    clk.at(idEpoch)   = arc.at(arcEpoch).clockError;
      observationTimes.push_back(arc.at(arcEpoch).time);

      const std::vector<GnssType> receiverTypes = definedTypes(times.at(idEpoch));

      // create observation class for each satellite
      UInt idObs  = 0;
      for(UInt k=0; k<arc.at(arcEpoch).satellite.size(); k++)
      {
        // find list of observation types for this satellite
        GnssType satType = arc.at(arcEpoch).satellite.at(k);
        UInt idType = 0;
        while(arc.at(arcEpoch).obsType.at(idType) != satType)
          idType++;

        // search transmitter index for satellite number (PRN)
        const UInt idTrans = std::distance(transmitters.begin(), std::find_if(transmitters.begin(), transmitters.end(),
                                                                              [&](auto t) {return t->PRN() == satType;}));
        std::vector<GnssType> transmitterTypes;
        if(idTrans < transmitters.size())
          transmitterTypes = transmitters.at(idTrans)->definedTypes(times.at(idEpoch));

        // repair GLONASS frequency number
        if((satType == GnssType::GLONASS) && transmitterTypes.size() && (transmitterTypes.front().frequencyNumber() != 9999))
          satType.setFrequencyNumber(transmitterTypes.front().frequencyNumber());

        GnssObservation *obs = new GnssObservation();
        for(; (idType<arc.at(arcEpoch).obsType.size()) && (arc.at(arcEpoch).obsType.at(idType)==satType); idType++, idObs++)
          if((idTrans < transmitters.size()) && arc.at(arcEpoch).observation.at(idObs)  && !std::isnan(arc.at(arcEpoch).observation.at(idObs)))
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

            // check useType and ignoreType
            Bool use = (useType.size()==0) ? TRUE : FALSE;
            if(type.isInList(useType))
              use = TRUE;
            if(type.isInList(ignoreType))
              use = FALSE;

            // check against receiver and transmitter types
            if(use && receiverTypes.size() && !type.isInList(receiverTypes))
            {
              use = FALSE;
              removedTypes[type]++;
            }
            if(use && transmitterTypes.size())
              for(const GnssType &typeTrans : GnssType::replaceCompositeSignals({type}))
                if(!typeTrans.isInList(transmitterTypes))
                {
                  use = FALSE;
                  removedTypes[type]++;
                }

            if(use)
              obs->push_back(GnssSingleObservation(type, arc.at(arcEpoch).observation.at(idObs)));
          }

        std::vector<GnssType> types;
        if((obs->size() == 0) || (idTrans >= transmitters.size()) ||
           !obs->init(*this, *transmitters.at(idTrans), rotationCrf2Trf, idEpoch, elevationCutOff, phaseWindup(idTrans)) ||
           !obs->observationList(group, types))
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

      if((observations_.size() <= idEpoch) || (observations_[idEpoch].size() == 0))
        disable(idEpoch);
      idEpoch++;
    } // for(arcEpoch)

    for(UInt idEpoch=observations_.size(); idEpoch<times.size(); idEpoch++)
      disable(idEpoch);

    if(removedTypes.size())
    {
      std::stringstream ss;
      for(const auto &type : removedTypes)
        ss<<"  "<<type.first.str()<<"="<<type.second;
      logWarning<<name()<<": removed undefined observations"<<ss.str()<<Log::endl;
    }

    observationSampling = medianSampling(observationTimes).seconds();
    copyObservations2ContinuousMemoryBlock();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssReceiver::simulateObservations(const std::vector<GnssType> &types,
                                        NoiseGeneratorPtr noiseClock, NoiseGeneratorPtr noiseObs,
                                        const std::vector<GnssTransmitterPtr> &transmitters,
                                        std::function<Rotary3d(const Time &time)> rotationCrf2Trf,
                                        const std::function<void(GnssObservationEquation &eqn)> &reduceModels,
                                        UInt minObsCountPerTrack, Angle elevationCutOff, Angle elevationTrackMinimum,
                                        const std::vector<GnssType> &useType, const std::vector<GnssType> &ignoreType, GnssObservation::Group group)
{
  try
  {
    // Simulate clock error
    // --------------------
    Vector clock = noiseClock->noise(times.size());
    for(UInt idEpoch=0; idEpoch<times.size(); idEpoch++)
      if(useable(idEpoch))
        updateClockError(idEpoch, clock(idEpoch)/LIGHT_VELOCITY);

    // Simulate zero observations
    // --------------------------
    Vector phaseWindup(transmitters.size());
    std::vector<std::vector<GnssType>> typesTrans(transmitters.size());
    for(UInt idEpoch=0; idEpoch<times.size(); idEpoch++)
    {
      const std::vector<GnssType> receiverTypes = definedTypes(times.at(idEpoch));

      // create observation class for each satellite
      for(UInt idTrans=0; idTrans<transmitters.size(); idTrans++)
      {
        // find list of observation types for this satellite
        GnssType satType = transmitters.at(idTrans)->PRN();
        std::vector<GnssType> transmitterTypes;
        if(idTrans < transmitters.size())
          transmitterTypes = transmitters.at(idTrans)->definedTypes(times.at(idEpoch));

        // repair GLONASS frequency number
        if((satType == GnssType::GLONASS) && transmitterTypes.size() && (transmitterTypes.front().frequencyNumber() != 9999))
          satType.setFrequencyNumber(transmitterTypes.front().frequencyNumber());

        GnssObservation *obs = new GnssObservation();
        for(UInt idType=0; idType<types.size(); idType++)
          if(types.at(idType) == satType)
          {
            GnssType type = types.at(idType) + satType;
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

            // check useType and ignoreType
            Bool use = (useType.size()==0) ? TRUE : FALSE;
            if(type.isInList(useType))
              use = TRUE;
            if(type.isInList(ignoreType))
              use = FALSE;

            // check against receiver and transmitter types
            if(use && receiverTypes.size() && !type.isInList(receiverTypes))
              use = FALSE;
            if(use && transmitterTypes.size())
              for(const GnssType &typeTrans : GnssType::replaceCompositeSignals({type}))
                if(!typeTrans.isInList(transmitterTypes))
                  use = FALSE;
            if(!use)
              continue;

            obs->push_back(GnssSingleObservation(type, 0.0));
            if(!type.isInList(typesTrans.at(idTrans)))
              typesTrans.at(idTrans).push_back(type);
          }

        std::vector<GnssType> types;
        if((obs->size() == 0) || (idTrans >= transmitters.size()) ||
           !obs->init(*this, *transmitters.at(idTrans), rotationCrf2Trf, idEpoch, elevationCutOff, phaseWindup(idTrans)) ||
           !obs->observationList(group, types))
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

      if((observations_.size() <= idEpoch) || (observations_[idEpoch].size() == 0))
        disable(idEpoch);
    } // for(arcEpoch)

    observationSampling = medianSampling(times).seconds();
    copyObservations2ContinuousMemoryBlock();

    // ambiguities
    // -----------
    class Ambiguity : public GnssAmbiguity
    {
    public:
      std::vector<GnssType> types;
      Vector                value; // ambiguities in meter

      explicit Ambiguity(GnssTrack *track, const Vector &value) : GnssAmbiguity(track), types(track->types), value(value) {}

      Vector ambiguities(const std::vector<GnssType> &types) const override
      {
        Vector value(types.size());
        UInt idx;
        for(UInt idType=0; idType<types.size(); idType++)
          if(types.at(idType).isInList(this->types, idx))
            value(idType) = this->value(idx);
        return value;
      }
    };

    // init random phase ambiguities
    std::random_device randomDevice;
    std::mt19937_64 generator; // for ambiguities
    generator.seed(randomDevice());
    auto ambiguityRandom = std::uniform_int_distribution<Int>(-10000, 10000);

    createTracks(transmitters, minObsCountPerTrack, {});
    for(auto &track : tracks)
    {
      Vector value(track->types.size());
      for(UInt i=0; i<value.size(); i++)
        value(i) = wavelengthFactor*track->types.at(i).wavelength() * ambiguityRandom(generator); // cycles to meter
      new Ambiguity(track.get(), value); // track is owner of ambiguity
    }

    // reduced observations
    // --------------------
    ObservationEquationList eqnList(*this, transmitters, rotationCrf2Trf, reduceModels, group);
    removeLowElevationTracks(eqnList, elevationTrackMinimum);

    for(UInt idTrans=0; idTrans<transmitters.size(); idTrans++)
      if(typesTrans.at(idTrans).size())
      {
        const Matrix eps = noiseObs->noise(times.size(), typesTrans.at(idTrans).size()); // obs noise
        UInt idx;
        for(UInt idEpoch=0; idEpoch<times.size(); idEpoch++)
          if(eqnList(idTrans, idEpoch))
          {
            const GnssObservationEquation &eqn = *eqnList(idTrans, idEpoch);
            GnssObservation *obs = observation(idTrans, idEpoch);
            for(UInt idType=0; idType<obs->size(); idType++)
              if(obs->at(idType).type.isInList(eqn.types, idx))
                obs->at(idType).observation = -eqn.l(idx) + eqn.sigma(idx) * eps(idEpoch, GnssType::index(typesTrans.at(idTrans), obs->at(idType).type));
          }
      }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssReceiver::estimateInitialClockErrorFromCodeObservations(const std::vector<GnssTransmitterPtr> &transmitters, std::function<Rotary3d(const Time &time)> rotationCrf2Trf,
                                                                 const std::function<void(GnssObservationEquation &eqn)> &reduceModels,
                                                                 Double huber, Double huberPower, Double maxPosDiff, Bool estimateKinematicPosition)
{
  try
  {
    const UInt countStaticPosition = (estimateKinematicPosition ? 0 : 3);
    Double posDiffStatic = 0;
    Vector posDiff(idEpochSize());
    Vector clockDiff(idEpochSize()); // in meters

    for(UInt iter=0; iter<10; iter++)
    {
      // setup observation equations: position, clock
      // --------------------------------------------
      std::vector<Matrix> listl, listA;
      std::vector<Matrix> listlFull, listAFull, listBFull;
      std::vector<UInt>   listEpoch, listObsCount;
      UInt maxSat = 0;
      for(UInt idEpoch=0; idEpoch<idEpochSize(); idEpoch++)
        if(useable(idEpoch))
        {
          // count observations and setup observation equations for each transmitter
          UInt                                 obsCount = 0;
          std::vector<GnssObservationEquation> eqnList;
          std::vector<GnssType>                types, systems;
          for(UInt idTrans=0; idTrans<idTransmitterSize(idEpoch); idTrans++)
          {
            GnssObservation *obs = observation(idTrans, idEpoch);
            if(obs && obs->observationList(GnssObservation::RANGE, types))
            {
              eqnList.emplace_back(*obs, *this, *transmitters.at(idTrans), rotationCrf2Trf, reduceModels, idEpoch, TRUE, types);
              eqnList.back().eliminateGroupParameters();
              obsCount += eqnList.back().l.rows();
              if(!types.at(0).isInList(systems))
                systems.push_back(types.at(0) & GnssType::SYSTEM);
            }
          }

          // setup combined observation equations
          Vector l(obsCount);
          Matrix A(obsCount, countStaticPosition);                  // pos (if static)
          Matrix B(obsCount, systems.size()+3-countStaticPosition); // clock (per system) + pos (if kinematic)

          if(!obsCount || (eqnList.size() <= systems.size()+3-countStaticPosition)) // if not enough observations -> delete epoch
          {
            disable(idEpoch);
            continue;
          }

          UInt idx = 0;
          for(const auto &eqn : eqnList)
          {
            const UInt count = eqn.l.rows();
            copy(eqn.l, l.row(idx, count));
            copy(eqn.A.column(GnssObservationEquation::idxClockRecv),  B.slice(idx, GnssType::index(systems, eqn.types.at(0)), count, 1)); // clock
            copy(eqn.A.column(GnssObservationEquation::idxPosRecv, 3), (estimateKinematicPosition ? B.slice(idx, systems.size(), count, 3) : A.slice(idx, 0, count, 3))); // pos
            idx += count;
          }

          listEpoch.push_back(idEpoch);
          listlFull.push_back(l);
          listAFull.push_back(A);
          listBFull.push_back(B);
          maxSat = std::max(maxSat, eqnList.size());

          if(countStaticPosition)
          {
            eliminationParameter(B, {A, l});
            listl.push_back(l);
            listA.push_back(A);
            listObsCount.push_back(l.rows());
          }
        } // for(idEpoch)

      if(!listEpoch.size() || (maxSat < 4))
      {
        logWarning<<"estimateInitialClockError: "<<name()<<" disabled because only "<<maxSat<<" satellites tracked"<<Log::endl;
        disable();
        return;
      }

      // estimate static parameters
      // --------------------------
      Vector x;
      if(countStaticPosition)
      {
        Vector sigma;
        x = robustLeastSquares(listA, listl, listObsCount, huber, huberPower, 30, sigma);
        posDiffStatic = norm(x.row(0,3));
      }

      // reconstruct kinematic parameters
      // --------------------------------
      for(UInt i=0; i<listEpoch.size(); i++)
      {
        if(x.size()) // update with static parameters
          matMult(-1, listAFull.at(i), x, listlFull.at(i));

        Vector sigma;
        const Vector y = Vce::robustLeastSquares(listBFull.at(i), listlFull.at(i), 1, huber, huberPower, 10, sigma);

        // set receiver clock error as mean of all system-specific clocks
        const UInt clockCount = y.size()-(estimateKinematicPosition ? 3 : 0);
        clockDiff(listEpoch.at(i)) = mean(y.row(0, clockCount));
        updateClockError(listEpoch.at(i), clockDiff(listEpoch.at(i))/LIGHT_VELOCITY);
        if(estimateKinematicPosition)
          posDiff(listEpoch.at(i)) = norm(y.row(clockCount, 3));
      }

      if(maxabs(clockDiff) <= maxPosDiff)
        break;
    } // for(iter)

    if(posDiffStatic > maxPosDiff)
    {
      logWarning<<"estimateInitialClockError: "<<name()<<" disabled due to position differences: "<<posDiffStatic<<" m"<<Log::endl;
      disable();
      return;
    }

    for(UInt idEpoch=0; idEpoch<idEpochSize(); idEpoch++)
      if((std::fabs(clockDiff(idEpoch)) > maxPosDiff) || (posDiff(idEpoch) > maxPosDiff))
      {
        logWarning<<"estimateInitialClockError: "<<name()<<" disabled at "<<times.at(idEpoch).dateTimeStr()<<" due to differences: pos="<<posDiff(idEpoch)<<" m, clock="<<clockDiff(idEpoch)<<" m"<<Log::endl;
        disable(idEpoch);
      }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssReceiver::disableEpochsWithGrossCodeObservationOutliers(ObservationEquationList &eqnList, Double threshold, Double outlierRatio)
{
  try
  {
    UInt disabledEpochs = 0;
    UInt outlierTotal = 0;
    UInt countTotal   = 0;
    for(UInt idEpoch=0; idEpoch<idEpochSize(); idEpoch++)
      if(useable(idEpoch))
      {
        // delete all observations to a satellite at epoch if they contain a gross code outlier
        UInt outlierCount = 0;
        UInt count   = 0;
        for(UInt idTrans=0; idTrans<idTransmitterSize(idEpoch); idTrans++)
          if(eqnList(idTrans, idEpoch))
          {
            const GnssObservationEquation &eqn = *eqnList(idTrans, idEpoch);
            for(UInt idType=0; idType<eqn.types.size(); idType++)
              if((eqn.types.at(idType) == GnssType::RANGE) && (std::fabs(eqn.l.at(idType)) >= threshold))
              {
                deleteObservation(idTrans, idEpoch);
                eqnList.deleteObservationEquation(idTrans, idEpoch);
                outlierCount++;
                break;
              }
            count++;
          }

        // disable receiver at epoch if 50% or more of the observed satellites have gross code outliers
        if(outlierCount >= outlierRatio * count)
        {
          for(UInt idTrans=0; idTrans<idTransmitterSize(idEpoch); idTrans++)
            eqnList.deleteObservationEquation(idTrans, idEpoch);
          disable(idEpoch);
          disabledEpochs++;
        }

        outlierTotal += outlierCount;
        countTotal   += count;
      }

//     if(disabledEpochs)
//       logWarning<<"grossCodeOutlier: "<<name()<<" deleted "<<outlierTotal<<" of "<<countTotal<<" observations and disabled "<<disabledEpochs<<" epochs"<<Log::endl;
//     else if(outlierTotal > 0.001*countTotal)
//       logWarning<<"grossCodeOutlier: "<<name()<<" deleted "<<outlierTotal<<" of "<<countTotal<<" observations"<<Log::endl;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssReceiver::createTracks(const std::vector<GnssTransmitterPtr> &transmitters, UInt minObsCountPerTrack, const std::vector<GnssType> &extraTypes)
{
  try
  {
    tracks.clear();
    for(UInt idTrans=0; idTrans<transmitters.size(); idTrans++)
    {
      UInt idEpochStart = 0;
      for(;;)
      {
        // find continuous track
        // ---------------------
        // find first epoch of the track
        while((idEpochStart < idEpochSize()) && !observation(idTrans, idEpochStart))
          idEpochStart++;

        GnssObservation *obs = observation(idTrans, idEpochStart);
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
          GnssObservation *obs = observation(idTrans, idEpoch);
          if(obs)
          {
            if((times.at(idEpoch)-times.at(idEpochEnd)).seconds() > 1.5*observationSampling)
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
        std::vector<GnssType> typeFrequencies;
        for(GnssType type : types)
          if(!type.isInList(extraTypes) && !type.isInList(typeFrequencies))
            typeFrequencies.push_back(type & GnssType::FREQUENCY);

        // define track
        if((countEpoch >= minObsCountPerTrack) && (typeFrequencies.size() >= 2))
        {
          tracks.push_back(std::make_shared<GnssTrack>(this, transmitters.at(idTrans).get(), idEpochStart, idEpochEnd, types));
          for(UInt idEpoch=idEpochStart; idEpoch<=idEpochEnd; idEpoch++)
            if(observation(idTrans, idEpoch))
              observation(idTrans, idEpoch)->track = tracks.back().get();
        }
        else
          for(UInt idEpoch=idEpochStart; idEpoch<=idEpochEnd; idEpoch++)
            deleteObservation(idTrans, idEpoch);

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

void GnssReceiver::deleteTrack(ObservationEquationList &eqnList, UInt idTrack)
{
  try
  {
    for(UInt idEpoch=tracks.at(idTrack)->idEpochStart; idEpoch<=tracks.at(idTrack)->idEpochEnd; idEpoch++)
    {
      deleteObservation(tracks.at(idTrack)->transmitter->idTrans(), idEpoch);
      eqnList.deleteObservationEquation(tracks.at(idTrack)->transmitter->idTrans(), idEpoch);
    }
    tracks.erase(tracks.begin()+idTrack);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssReceiver::deleteEmptyTracks()
{
  try
  {
    auto isEmpty = [](GnssTrackPtr t)
    {
      for(UInt idEpoch=t->idEpochStart; idEpoch<=t->idEpochEnd; idEpoch++)
        if(t->transmitter->useable(idEpoch) && t->receiver->observation(t->transmitter->idTrans(), idEpoch))
          return FALSE;
      return TRUE;
    };

    tracks.erase(std::remove_if(tracks.begin(), tracks.end(), isEmpty), tracks.end());
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssReceiver::removeLowElevationTracks(ObservationEquationList &eqnList, Angle minElevation)
{
  try
  {
    for(UInt idTrack=tracks.size(); idTrack-->0;)
    {
      const UInt idTrans = tracks.at(idTrack)->transmitter->idTrans();
      Bool removeTrack = TRUE;
      for(UInt idEpoch=tracks.at(idTrack)->idEpochStart; idEpoch<=tracks.at(idTrack)->idEpochEnd; idEpoch++)
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

GnssTrackPtr GnssReceiver::splitTrack(ObservationEquationList &eqnList, GnssTrackPtr track, UInt idEpochSplit)
{
  try
  {
    // new track
    const UInt idTrans = track->transmitter->idTrans();
    GnssTrackPtr trackNew = std::make_shared<GnssTrack>(track->receiver, track->transmitter, idEpochSplit, track->idEpochEnd, track->types);
    tracks.push_back(trackNew);

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

// determine Melbourne-Wuebbena-like linear combinations
void GnssReceiver::linearCombinations(ObservationEquationList &eqnList, GnssTrackPtr track, const std::vector<GnssType> &extraTypes,
                                      std::vector<GnssType> &typesPhase, std::vector<UInt> &idEpochs, Matrix &combinations, Double &cycles2tecu) const
{
  try
  {
    // available observations for this track
    typesPhase.clear();
    idEpochs.clear();
    for(UInt idEpoch=track->idEpochStart; idEpoch<=track->idEpochEnd; idEpoch++)
      if(eqnList(track->transmitter->idTrans(), idEpoch) && observation(track->transmitter->idTrans(), idEpoch))
      {
        idEpochs.push_back(idEpoch);
        for(GnssType type : eqnList(track->transmitter->idTrans(), idEpoch)->types)
          if((type == GnssType::PHASE) && !type.isInList(typesPhase) && !type.isInList(extraTypes)) // ignore time variable GPS L5 signals
            typesPhase.push_back(type);
      }

    const Matrix Bias = GnssLambda::phaseDecorrelation(typesPhase, wavelengthFactor);
    combinations = Matrix(idEpochs.size(), Bias.columns()-1);
    UInt row = 0;
    for(UInt idEpoch : idEpochs)
    {
      const GnssObservationEquation &eqn = *eqnList(track->transmitter->idTrans(), idEpoch);
      Vector l = eqn.l;
      Matrix A(l.rows(), Bias.columns()+2);
      UInt idx;
      for(UInt idType=0; idType<eqn.types.size(); idType++) // ambiguities
        if(eqn.types.at(idType).isInList(typesPhase, idx) || (eqn.types.at(idType) == GnssType::RANGE))
        {
          l(idType) = eqn.l(idType)/eqn.sigma0(idType);
          A(idType, 0) = 1.; // range
          A(idType, 1) = eqn.types.at(idType).ionosphericFactor(); // TEC
          if(idx != NULLINDEX)
            copy(Bias.row(idx), A.slice(idType, 2, 1, Bias.columns()));
          A.row(idType) *= 1./eqn.sigma0(idType);
        }

      // skip the first, inaccurate one
      copy(leastSquares(A, l).row(2+1, Bias.columns()-1).trans(), combinations.row(row++));
    }

    // determine cycle slip size in terms of TEC
    Vector l = Bias.column(0); // cycle slips can only occur in this linear combination anymore
    Matrix A(typesPhase.size(), 2, 1.); // first column range
    for(UInt idType=0; idType<typesPhase.size(); idType++)
      A(idType, 1) = typesPhase.at(idType).ionosphericFactor(); // TEC
    cycles2tecu = std::fabs(leastSquares(A, l)(1,0)); // one cycle slip in terms of TEC
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

// determine range & TEC based on phase observations only
void GnssReceiver::rangeAndTec(ObservationEquationList &eqnList, UInt idTrans, const std::vector<UInt> &idEpochs,
                               const std::vector<GnssType> &typesPhase, Vector &range, Vector &tec) const
{
  try
  {
    range = Vector(idEpochs.size());
    tec   = Vector(idEpochs.size());
    for(UInt i=0; i<idEpochs.size(); i++)
    {
      const GnssObservationEquation &eqn = *eqnList(idTrans, idEpochs.at(i));
      Vector l(typesPhase.size());
      for(UInt k=0; k<l.rows(); k++)
        l(k) = eqn.l(GnssType::index(eqn.types, typesPhase.at(k)));
      Matrix A(typesPhase.size(), 2, 1.);
      for(UInt k=0; k<l.rows(); k++)
        A(k, 1) = typesPhase.at(k).ionosphericFactor();
      const Vector x = leastSquares(A, l);
      range(i) = x(0);
      tec(i)   = x(1);
    } // for(idEpoch)
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssReceiver::writeTracks(const FileName &fileName, ObservationEquationList &eqnList, const std::vector<GnssType> &extraTypes) const
{
  try
  {
    if(fileName.empty())
      return;

    for(const auto &track : tracks)
      if(track->countObservations())
      {
        std::vector<GnssType> typesPhase;
        std::vector<UInt>     idEpochs;
        Matrix                combinations;
        Double                cycles2tecu;
        Vector                range, tec;
        linearCombinations(eqnList, track, extraTypes, typesPhase, idEpochs, combinations, cycles2tecu);
        rangeAndTec(eqnList, track->transmitter->idTrans(), idEpochs, typesPhase, range, tec);

        Matrix A(idEpochs.size(), 2+combinations.columns());
        axpy(1./cycles2tecu, tec, A.column(1));
        copy(combinations, A.column(2, combinations.columns()));

        std::vector<Time> timesTrack;
        for(UInt idEpoch : idEpochs)
          timesTrack.push_back(times.at(idEpoch));
        for(UInt i=1; i<A.columns(); i++)
          A.column(i) -= median(A.column(i));

        std::string typeStr;
        for(GnssType type : typesPhase)
          typeStr += type.str().substr(0, 3);
        typeStr = String::replaceAll(typeStr, "?", "");
        VariableList varList;
        addVariable("station",        name(),                     varList);
        addVariable("prn",            track->transmitter->name(), varList);
        addVariable("trackTimeStart", timesTrack.front().mjd(),   varList);
        addVariable("trackTimeEnd",   timesTrack.back().mjd(),    varList);
        addVariable("types",          typeStr,                    varList);
        InstrumentFile::write(fileName(varList), Arc(timesTrack, A));
      }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssReceiver::cycleSlipsDetection(ObservationEquationList &eqnList, UInt minObsCountPerTrack, Double lambda, UInt windowSize, Double tecSigmaFactor, const std::vector<GnssType> &extraTypes)
{
  try
  {
    for(UInt idTrack=0; idTrack<tracks.size(); idTrack++)
    {
      if(tracks.at(idTrack)->countObservations() >= std::max(minObsCountPerTrack, windowSize))
        cycleSlipsDetection(eqnList, tracks.at(idTrack), lambda, windowSize, tecSigmaFactor, extraTypes);
      if(tracks.at(idTrack)->countObservations() < minObsCountPerTrack)
        deleteTrack(eqnList, idTrack--);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssReceiver::cycleSlipsDetection(ObservationEquationList &eqnList, GnssTrackPtr track, Double lambda, UInt windowSize, Double tecSigmaFactor, const std::vector<GnssType> &extraTypes)
{
  try
  {
    // determine Melbourne-Wuebbena-like linear combinations
    // -----------------------------------------------------
    std::vector<GnssType> typesPhase;
    std::vector<UInt>     idEpochs;
    Matrix                combinations;
    Double                cycles2tecu;
    linearCombinations(eqnList, track, extraTypes, typesPhase, idEpochs, combinations, cycles2tecu);

    Vector slips(idEpochs.size());
    for(UInt k=0; k<combinations.columns(); k++)
    {
      Vector smoothed = combinations.column(k);
      if(k == 0)
        smoothed = totalVariationDenoising(smoothed, lambda);
      for(UInt i=1; i<idEpochs.size(); i++) // cycle slip if denoised difference exceeds 3/4 cycle
        if(std::fabs(smoothed(i)-smoothed(i-1)) > 0.75)
          slips(i) = TRUE;
    }

    for(UInt i=slips.rows(); i-->0;)
      if(slips(i))
      {
        splitTrack(eqnList, track, idEpochs.at(i));
        idEpochs.resize(i); // shorten data to new track length
      }

    // find cycle slips in TEC based on moving window over autoregressive model residuals
    // ----------------------------------------------------------------------------------
    Vector range, tec;
    rangeAndTec(eqnList, track->transmitter->idTrans(), idEpochs, typesPhase, range, tec);
    if(windowSize)
    {
      const UInt order = 3; // AR model order
      std::vector<UInt> slips;
      if(tec.size() >= order+windowSize)
      {
        // high pass filter via AR model
        Vector l = tec.row(order, tec.size()-order);
        Matrix A = Matrix(l.rows(), order);
        for(UInt k=0; k<order; k++)
          copy(tec.row(order-k-1, tec.rows()-order), A.column(k));
        leastSquares(A, l); // l contains AR model residuals after function call

        // peak/outlier detection using moving standard deviation over AR model residuals.
        // automatic threshold scaling via standard deviation is used to prevent excessive
        // splitting during periods with high ionospheric variations/scintillations
        for(UInt idxStart=0; idxStart<l.size()-windowSize; idxStart++)
          if(l(idxStart+windowSize/2) > std::max(0.9*cycles2tecu, tecSigmaFactor*standardDeviation(l.row(idxStart, windowSize))))
            slips.push_back(idxStart+windowSize/2+order);
      }

      for(UInt i=slips.size(); i-->0;)
      {
        splitTrack(eqnList, track, idEpochs.at(slips.at(i)));
        idEpochs.resize(slips.at(i)); // shorten data to new track length
      }
    }

    // repair GPS L5 cycle slips
    // -------------------------
    for(GnssType type : extraTypes)
    {
      const UInt   idTrans    = track->transmitter->idTrans();
      const Double wavelength = wavelengthFactor * type.wavelength();
      const Double TEC        = type.ionosphericFactor();
      for(UInt idType=0; idType<track->types.size(); idType++)
        if(track->types.at(idType) == type)
        {
          // reduce l by estimated range and tec
          Vector l(idEpochs.size());
          for(UInt i=0; i<idEpochs.size(); i++)
            l(i) = eqnList(idTrans, idEpochs.at(i))->l(GnssType::index(eqnList(idTrans, idEpochs.at(i))->types, track->types.at(idType))) - range(i) - TEC * tec(i);

          // fix jumps
          Double jump = 0;
          for(UInt i=1; i<l.rows(); i++)
          {
            jump += wavelength * std::round((l(i)-l(i-1))/wavelength);
            eqnList(idTrans, idEpochs.at(i))->l(GnssType::index(eqnList(idTrans, idEpochs.at(i))->types, track->types.at(idType))) -= jump;
            observation(idTrans, idEpochs.at(i))->at(track->types.at(idType)).observation -= jump;
          }
        }
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssReceiver::cycleSlipsRepairAtSameFrequency(ObservationEquationList &eqnList)
{
  try
  {
    // get all phase types
    std::vector<GnssType> types;
    for(UInt idEpoch=0; idEpoch<idEpochSize(); idEpoch++)
      for(UInt idTrans=0; idTrans<idTransmitterSize(idEpoch); idTrans++)
        if(eqnList(idTrans, idEpoch))
        {
          const GnssObservationEquation &eqn = *eqnList(idTrans, idEpoch);
          for(UInt idType=0; idType<eqn.types.size(); idType++)
            if((eqn.types.at(idType) == GnssType::PHASE) && !eqn.types.at(idType).isInList(types))
              types.push_back(eqn.types.at(idType) & ~(GnssType::PRN+GnssType::FREQ_NO));
        }
    std::sort(types.begin(), types.end());

    // find two phase observations with same system and frequency
    for(UInt idType=1; idType<types.size(); idType++)
      if(types.at(idType) == (types.at(idType-1) & (GnssType::FREQUENCY + GnssType::SYSTEM)))
      {
        // compute difference
        std::vector<UInt>   idxTrans, idxEpoch;
        std::vector<Double> values;
        for(UInt idEpoch=0; idEpoch<idEpochSize(); idEpoch++)
          for(UInt idTrans=0; idTrans<idTransmitterSize(idEpoch); idTrans++)
            if(eqnList(idTrans, idEpoch))
            {
              const GnssObservationEquation &eqn = *eqnList(idTrans, idEpoch);
              UInt idx, idx1;
              if(!types.at(idType).isInList(eqn.types, idx) || !types.at(idType-1).isInList(eqn.types, idx1))
                continue;
              idxTrans.push_back(idTrans);
              idxEpoch.push_back(idEpoch);
              values.push_back((eqn.l(idx)-eqn.l(idx1))/eqn.types.at(idx).wavelength()); // diff in cycles
            }

        if(!values.size())
          throw(Exception(name()+": "+types.at(idType).str()+" and "+types.at(idType-1).str()+" must be observed together"));

        // consider bias (e.g. quarter cycles)
        Vector v0s(values.size());
        for(UInt i=0; i<values.size(); i++)
          v0s(i) = values.at(i)-std::round(values.at(i));
        const Double v0 = median(v0s);

        if(std::fabs(std::fabs(v0)-0.5) < 0.05)
          logWarning<<name()<<": a phase bias difference of near a half cycle ("<<v0<<") between "<<types.at(idType).str()<<" and "<<types.at(idType-1).str()<<" might causes problems"<<Log::endl;

        // fix jumps
        for(UInt i=0; i<values.size(); i++)
        {
          GnssObservationEquation &eqn = *eqnList(idxTrans.at(i), idxEpoch.at(i));
          const UInt   idx = GnssType::index(eqn.types, types.at(idType));
          const Double v   = eqn.types.at(idx).wavelength() * std::round(values.at(i)-v0);
          eqn.l(idx) -= v;
          observation(idxTrans.at(i), idxEpoch.at(i))->at(eqn.types.at(idx)).observation -= v;
        }
      } // for(idType)
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssReceiver::trackOutlierDetection(const ObservationEquationList &eqnList, const std::vector<GnssType> &ignoreTypes, Double huber, Double huberPower)
{
  try
  {
    for(auto track : tracks)
    {
      const UInt idTrans      = track->transmitter->idTrans();
      const UInt idEpochStart = track->idEpochStart;
      const UInt idEpochEnd   = track->idEpochEnd;

      // available observations for this track
      std::vector<GnssType> types;
      for(UInt idEpoch=idEpochStart; idEpoch<=idEpochEnd; idEpoch++)
        if(eqnList(idTrans, idEpoch))
          for(const GnssType &type : eqnList(idTrans, idEpoch)->types)
            if(!type.isInList(types) && !type.isInList(ignoreTypes))
              types.push_back(type);

      // determine biases (reduced by range and TEC)
      Matrix Bias = identityMatrix(types.size());
      Matrix B(types.size(), 2);
      for(UInt idType=0; idType<types.size(); idType++)
        if(types.at(idType) == GnssType::RANGE)
        {
          B(idType, 0) = 1.; //range
          B(idType, 1) = types.at(idType).ionosphericFactor(); // TEC
        }
      const Vector tau = QR_decomposition(B);
      QMult(B, tau, Bias);
      Bias = Bias.column(B.columns(), Bias.rows()-B.columns());

      // setup observation equations: range, TEC, ambiguities
      // ----------------------------------------------------
      std::vector<Matrix> listl, listA;
      std::vector<UInt>   listEpoch, listObsCount;
      for(UInt idEpoch=idEpochStart; idEpoch<=idEpochEnd; idEpoch++)
        if(eqnList(idTrans, idEpoch))
        {
          const GnssObservationEquation &eqn = *eqnList(idTrans, idEpoch);

          // observations
          Vector l = eqn.l;

          // distance and TEC
          Matrix B(l.rows(), 2);
          copy(eqn.A.column(GnssObservationEquation::idxRange), B.column(0));
          copy(eqn.A.column(GnssObservationEquation::idxSTEC),  B.column(1));

          // signal biases (includes ambiguities)
          UInt idx;
          Matrix A(l.rows(), Bias.columns());
          for(UInt idType=0; idType<eqn.types.size(); idType++)
            if(eqn.types.at(idType).isInList(types, idx))
              matMult(1., eqn.A.column(GnssObservationEquation::idxUnit+idType), Bias.row(idx), A);

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
            if(eqn.types.at(idType).isInList(ignoreTypes))
            {
              l.row(idType).setNull();
              A.row(idType).setNull();
              B.row(idType).setNull();
              obsCount--;
            }

          eliminationParameter(B, A, l);
          listEpoch.push_back(idEpoch);
          listl.push_back(l);
          listA.push_back(A);
          listObsCount.push_back(obsCount-B.columns());
        } // for(idEpoch)

      // estimate solution
      // -----------------
      Vector sigma;
      Vector x = robustLeastSquares(listA, listl, listObsCount, huber, huberPower, 30, sigma);

      // downweight outliers
      // -------------------
      for(UInt i=0; i<listEpoch.size(); i++)
      {
        const GnssObservationEquation &eqn = *eqnList(idTrans, listEpoch.at(i));
        GnssObservation &obs = *observation(idTrans, listEpoch.at(i));
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
            const Double lambda = types.at(idType).wavelength();
            b(idType) = lambda * std::round(b(idType)/lambda);
            eqnList(idTrans, listEpoch.at(i))->l(GnssType::index(eqnList(idTrans, listEpoch.at(i))->types, types.at(idType))) -= b(idType);
            observation(idTrans, listEpoch.at(i))->at(types.at(idType)).observation -= b(idType);
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

Matrix GnssReceiver::robustLeastSquares(const std::vector<Matrix> &A, const std::vector<Matrix> &l, const std::vector<UInt> &obsCount,
                                        Double huber, Double huberPower, UInt maxIter, Vector &sigma)
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
    UInt   countOutlier = 0;
    Double sigma0 = 0;
    sigma = Vector(countEpoch, 1.);
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

      // QR decomposition
      Vector tau = QR_decomposition(WA);
      QTransMult(WA, tau, Wl); // transform observations: l:= Q'l
      x = Wl.row(0, WA.columns());
      triangularSolve(1., WA.row(0, WA.columns()), x);
      Wl.row(0, WA.columns()).setNull(); // residuals: remove WB*x
      QMult(WA, tau, Wl); // back transformation
      generateQ(WA, tau); // for redundancies

      if(sigma0 == 0.)
        sigma0 = std::sqrt(quadsum(Wl)/(Wl.size()-x.size()));

      // outlier detection
      UInt   countOutlierNew = 0;
      Double ePeSum = 0.;
      Double rSum   = 0.;
      index = 0;
      for(UInt i=0; i<sigma.rows(); i++)
      {
        const UInt   count = l.at(i).rows();
        const Double ePe   = quadsum(Wl.row(index, count))/Wl.columns();
        const Double r     = obsCount.at(i) - quadsum(WA.row(index, count));
        const Double s     = std::sqrt(ePe/r)*sigma(i)/sigma0;
        ePeSum += ePe;
        rSum   += r;
        sigma(i) = 1.;
        if((s > huber) && (r > 1e-4)) // redundancy: it is possible to estimate sigma?
        {
          sigma(i) = std::pow(s/huber, huberPower);
          countOutlierNew++;
        }
        index += count;
      }

      const Double sigma0New = Vce::standardDeviation(ePeSum, rSum, huber, huberPower);
      if((countOutlierNew == 0) || ((countOutlier == countOutlierNew) && (std::fabs(sigma0New-sigma0)/sigma0 < 0.001)))
        break;
      sigma0       = sigma0New;
      countOutlier = countOutlierNew;

//       if(iter >= maxIter-1)
//         logWarning<<"GnssReceiver::robustLeastSquares: no convergence, sigma="<<sigma0<<", outlier="<<countOutlier<<" of "<<countEpoch<<Log::endl;
    } // for(iter)

    return x;
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
Matrix GnssReceiver::totalVariationDenoising(const_MatrixSliceRef y, Double lambda)
{
  try
  {
    Matrix x(y.rows(), y.columns());
    for(UInt col=0; col<y.columns(); col++)
    {
      // initialize total variation denoising algorithm
      UInt N = y.rows()-1;
      UInt k  = 0;                      // current sample location
      UInt k0 = 0;                      // beginning of current segment
      UInt km = 0;                      // last position where umax = -lambda
      UInt kp = 0;                      // last position where umin =  lambda
      Double vMin = y(0, col) - lambda; // lower bound for the segment's value
      Double vMax = y(0, col) + lambda; // upper bound for the segment's value
      Double uMin =  lambda;            // u is the dual variable
      Double uMax = -lambda;            // u is the dual variable

      // total variation denoising algorithm
      for(;;)
      {
        if(k == N)
        {
          x(N, col) = vMin + uMin;
          break;
        }

        if(y(k+1, col) + uMin < vMin - lambda)      // negative jump necessary
        {
          for(UInt i=k0; i<=km; i++)
            x(i, col) = vMin;
          k = k0 = km = kp = km+1;
          vMin = y(k, col);
          vMax = y(k, col) + 2*lambda;
          uMin =  lambda;
          uMax = -lambda;
        }
        else if(y(k+1, col) + uMax > vMax + lambda) // positive jump necessary
        {
          for(UInt i=k0; i<=kp; i++)
            x(i, col) = vMax;
          k = k0 = km = kp = kp+1;
          vMin = y(k, col) - 2*lambda;
          vMax = y(k, col);
          uMin =  lambda;
          uMax = -lambda;
        }
        else  // no jump necessary
        {
          k = k+1;
          uMin = uMin + y(k, col) - vMin;
          uMax = uMax + y(k, col) - vMax;
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
          for(UInt i=k0; i<=km; i++)
            x(i, col) = vMin;
          k = k0 = km = km+1;
          vMin = y(k, col);
          uMin = lambda;
          uMax = y(k, col) + lambda - vMax;
          continue;
        }
        else if(uMax > 0.)  // vMax is too low ==> positive jump necessary
        {
          for(UInt i=k0; i<=kp; i++)
            x(i, col) = vMax;
          k = k0 = kp = kp+1;
          vMax = y(k, col);
          uMax = -lambda;
          uMin = y(k, col) - lambda - vMin;
          continue;
        }
        else
        {
          for(UInt i=k0; i<=N; i++)
            x(i, col) = vMin + uMin/(k-k0+1);
          break;
        }
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

GnssTrack::GnssTrack(GnssReceiver *_receiver, GnssTransmitter *_transmitter, UInt _idEpochStart, UInt _idEpochEnd, const std::vector<GnssType> &_types) :
    receiver(_receiver), transmitter(_transmitter), idEpochStart(_idEpochStart), idEpochEnd(_idEpochEnd), types(_types), ambiguity(nullptr)
{
}

/***********************************************/

GnssTrack::~GnssTrack()
{
  delete ambiguity;
}

/***********************************************/

UInt GnssTrack::countObservations() const
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

void GnssTrack::removeAmbiguitiesFromObservations(const std::vector<GnssType> &types, const std::vector<Double> &value)
{
  try
  {
    for(UInt idEpoch=idEpochStart; idEpoch<=idEpochEnd; idEpoch++)
      if(receiver->observation(transmitter->idTrans(), idEpoch))
      {
        GnssObservation &obs = *receiver->observation(transmitter->idTrans(), idEpoch);
        UInt idx;
        for(UInt idType=0; idType<obs.size(); idType++)
          if(obs.at(idType).type.isInList(types, idx))
            obs.at(idType).observation -= value.at(idx);
      }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
