/***********************************************/
/**
* @file gnssSimulateReceiver.cpp
*
* @brief GNSS receiver simulation.
*
* @author Sebastian Strasser
* @author Torsten Mayer-Gürr
* @date 2014-09-25
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program simulates observations from receivers to GNSS satellites.
These simulated observations can then be used in \program{GnssProcessing}, for example to conduct closed-loop simulations.

One or more GNSS constellations must be defined via \configClass{transmitter}{gnssParametrizationTransmitterType}.
Receivers such as ground station networks or Low Earth Orbit (LEO) satellites can be defined via \configClass{receiver}{gnssParametrizationReceiverType}.
In case of \configClass{receiver:stationNetwork}{gnssParametrizationReceiverType:stationNetwork} \configFile{outputfileGnssReceiver}{instrument} and
\configFile{outputfileClock}{instrument} are template file names with the \reference{variable}{general.parser} \verb|{station}| being replaced by the station name.

A list of simulated observation types can be defined via \configClass{observationType}{gnssType}. Noise can be added to both observations and clock errors
via \configClass{noiseObervation}{noiseGeneratorType} and \configClass{noiseClockReceiver}{noiseGeneratorType}, respectively. Observation noise is
interpreted as a factor that is multiplied to the accuracy derived from the accuracy pattern of the respective observation type
(see \configFile{inputfileAccuracyDefinition}{gnssAntennaDefinition} in \configClass{receiver}{gnssParametrizationReceiverType}).

If the program is run on multiple processes in \reference{parallel}{general.parallelization}, \config{parallelIntervals}=\verb|yes| distributes the
defined \configClass{intervals}{timeSeriesType} over all processes. Otherwise the intervals are processed consecutively while the
\configClass{receiver}{gnssParametrizationReceiverType}s (stations or LEO satellites) are distributed over the processes.
)";

/***********************************************/

#include <random>
#include "programs/program.h"
#include "base/planets.h"
#include "parser/dataVariables.h"
#include "files/fileStringTable.h"
#include "files/fileParameterName.h"
#include "files/fileMatrix.h"
#include "classes/timeSeries/timeSeries.h"
#include "classes/noiseGenerator/noiseGenerator.h"
#include "gnss/gnss.h"
#include "gnss/gnssObservation.h"
#include "gnss/gnssReceiver.h"
#include "gnss/gnssParametrizationTransmitter.h"
#include "gnss/gnssParametrizationReceiver.h"
#include "gnss/gnssParametrizationIonosphere.h"
#include "gnss/gnssParametrizationEarthRotation.h"

/***** CLASS ***********************************/

/** @brief GNSS receiver simulation.
* @ingroup programsGroup */
class GnssSimulateReceiver
{
  FileName                           fileNameReceiver, fileNameClock;
  std::vector<Time>                  timeSeries, timesInterval;
  VariableList                       fileNameVariableList;
  Double                             marginSeconds;
  Gnss                               gnss;
  GnssParametrizationIonospherePtr   gnssIonosphere;
  std::vector<GnssType>              typesObs;
  NoiseGeneratorPtr                  noiseClock, noiseObs;
  std::mt19937_64                    generator; // for ambiguities
  std::uniform_int_distribution<Int> ambiguityRandom;

  void computeInterval(UInt idInterval, Parallel::CommunicatorPtr comm);
  Bool computeObservation(const Gnss::Receiver &receiver, const Gnss::Transmitter &transmitter, GnssParametrizationIonospherePtr ionosphere,
                          UInt idEpoch, const std::vector<GnssType> &types, Double &phaseWindupOld,
                          Vector &obs, Vector &sigma);

public:
  void run(Config &config);
};

GROOPS_REGISTER_PROGRAM(GnssSimulateReceiver, PARALLEL, "GNSS receiver simulation", Gnss, Simulation)

/***********************************************/

void GnssSimulateReceiver::run(Config &config)
{
  try
  {
    TimeSeriesPtr                                  timeSeriesPtr, timesIntervalPtr;
    std::vector<GnssParametrizationTransmitterPtr> gnssTransmitters;
    std::vector<GnssParametrizationReceiverPtr>    gnssReceivers;
    GnssParametrizationEarthRotationPtr            gnssEarthRotation;
    Bool                                           parallelIntervals;

    readConfig(config, "outputfileGnssReceiver",  fileNameReceiver,  Config::MUSTSET,  "gnssReceiver_{loopTime:%D}.{station}.dat.gz", "simulated observations");
    readConfig(config, "outputfileClock",         fileNameClock,     Config::OPTIONAL, "clock_{loopTime:%D}.{station}.dat", "simulated receiver clock errors");
    readConfig(config, "intervals",               timesIntervalPtr,  Config::OPTIONAL, "",    "each interval is processed independently");
    readConfig(config, "timeSeries",              timeSeriesPtr,     Config::MUSTSET,  "",    "defines epochs within intervals");
    readConfig(config, "timeMargin",              marginSeconds,     Config::DEFAULT,  "0.1", "[seconds] margin to consider two times identical");
    readConfig(config, "transmitter",             gnssTransmitters,  Config::MUSTSET,  "",    "constellation of GNSS satellites");
    readConfig(config, "receiver",                gnssReceivers,     Config::MUSTSET,  "",    "ground station network or LEO satellite");
    readConfig(config, "ionosphere",              gnssIonosphere,    Config::MUSTSET,  "",    "ionosphere settings");
    readConfig(config, "earthRotation",           gnssEarthRotation, Config::MUSTSET,  "",    "Earth rotation");
    readConfig(config, "observationType",         typesObs,          Config::MUSTSET,  "",    "simulated observation types");
    readConfig(config, "noiseObservation",        noiseObs,          Config::DEFAULT,  "",    "[-] noise is multplied with type accuracy pattern of receiver");
    readConfig(config, "noiseClockReceiver",      noiseClock,        Config::DEFAULT,  "",    "[m] noise added to the simulated receiver clock");
    readConfig(config, "parallelIntervals",       parallelIntervals, Config::DEFAULT,  "1",   "parallelize intervals instead of receivers per interval");
    if(isCreateSchema(config)) return;

    // ============================

    // init random phase ambiguities
    std::random_device randomDevice;
    generator.seed(randomDevice()+1234*Parallel::myRank());
    ambiguityRandom = std::uniform_int_distribution<Int>(-10000, 10000);

    // ============================

    // init time intervals
    // -------------------
    timeSeries = timeSeriesPtr->times();
    if(timesIntervalPtr)
      timesInterval = timesIntervalPtr->times();
    if(timesInterval.size()==0)
    {
      timesInterval.push_back(timeSeries.at(0));
      timesInterval.push_back(timeSeries.back()+seconds2time(1));
    }

    addTimeVariables(fileNameVariableList);
    addVariable("station", fileNameVariableList);

    std::sort(typesObs.begin(), typesObs.end());

    // ============================

    // init
    // ----
    std::vector<Gnss::TransmitterPtr>     transmitters;
    std::vector<Gnss::ReceiverPtr>        receivers;
    std::vector<Gnss::ParametrizationPtr> parametrizations;

    for(auto &gnssTransmitter : gnssTransmitters)
    {
      auto t = gnssTransmitter->transmitters();
      auto p = gnssTransmitter->parametrizations();
      transmitters.insert(transmitters.end(), t.begin(), t.end());
      parametrizations.insert(parametrizations.end(), p.begin(), p.end());
    }

    for(auto &gnssReceiver : gnssReceivers)
    {
      auto r = gnssReceiver->receivers();
      auto p = gnssReceiver->parametrizations();
      receivers.insert(receivers.end(), r.begin(), r.end());
      parametrizations.insert(parametrizations.end(), p.begin(), p.end());
    }

    logInfo<<"Init GNSS"<<Log::endl;
    gnss.init(receivers, transmitters, parametrizations, gnssIonosphere, gnssEarthRotation, GnssParametrizationGravityFieldPtr(), GnssParametrizationAmbiguitiesPtr());

    logInfo<<"  transmitter: "<<transmitters.size()<<Log::endl;
    logInfo<<"  receiver:    "<<receivers.size()<<Log::endl;

    // ============================

    logStatus<<"compute intervals"<<Log::endl;
    if(parallelIntervals)
    {
      Parallel::forEach(timesInterval.size()-1, [this](UInt idInterval) {computeInterval(idInterval, Parallel::selfCommunicator());});
    }
    else
    {
      for(UInt idInterval=0; idInterval<timesInterval.size()-1; idInterval++)
        computeInterval(idInterval, nullptr); // use all processes for each interval
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssSimulateReceiver::computeInterval(UInt idInterval, Parallel::CommunicatorPtr comm)
{
  try
  {
    logStatus<<"================================================================"<<Log::endl;
    logStatus<<idInterval+1<<". "<<timesInterval.at(idInterval).dateTimeStr()<<" - "<<timesInterval.at(idInterval+1).dateTimeStr()<<Log::endl;
    evaluateTimeVariables(idInterval, timesInterval.at(idInterval), timesInterval.at(idInterval+1), fileNameVariableList);

    // init time series
    // ----------------
    std::vector<Time> times;
    for(UInt i=0; i<timeSeries.size(); i++)
      if(timeSeries.at(i).isInInterval(timesInterval.at(idInterval), timesInterval.at(idInterval+1)))
        times.push_back(timeSeries.at(i));

    if(!times.size())
    {
      logWarning<<"empty interval"<<Log::endl;
      return;
    }

    // init the GNSS system
    // --------------------
    Gnss::AnalysisType analysisType = Gnss::ANALYSIS_RAW;
    gnss.initInterval(analysisType, times, seconds2time(marginSeconds), comm);


    // ==========================================

    // compute station wise solutions
    // ------------------------------
    logStatus<<"simulate observations"<<Log::endl;
    for(UInt idRecv=0; idRecv<gnss.receiver.size(); idRecv++)
      if(gnss.receiver.at(idRecv)->useable())
      {
        logStatus<<"################ station: "<<idRecv<<" ("<<gnss.receiver.at(idRecv)->name()<<") #############"<<Log::endl;
        fileNameVariableList["station"]->setValue(gnss.receiver.at(idRecv)->name());

        // Simulate clock error
        // --------------------
        Vector clock = noiseClock->noise(gnss.times.size());
        for(UInt idEpoch=0; idEpoch<gnss.times.size(); idEpoch++)
          gnss.receiver.at(idRecv)->updateClockError(idEpoch, clock(idEpoch)/LIGHT_VELOCITY);

        if(!fileNameClock.empty())
        {
          logStatus<<"write clock errors to <"<<fileNameClock(fileNameVariableList)<<">"<<Log::endl;
          MiscValueArc arc;
          for(UInt idEpoch=0; idEpoch<gnss.times.size(); idEpoch++)
          {
            MiscValueEpoch epoch;
            epoch.time  = gnss.times.at(idEpoch);
            epoch.value = gnss.receiver.at(idRecv)->clockError(idEpoch);
            arc.push_back(epoch);
          }
          InstrumentFile::write(fileNameClock(fileNameVariableList), arc);
        }

        // Prepare observations
        GnssReceiverArc receiverArc;
        for(UInt idEpoch=0; idEpoch<gnss.times.size(); idEpoch++)
        {
          GnssReceiverEpoch epoch;
          epoch.time       = gnss.receiver.at(idRecv)->timeReceiver(idEpoch);
          epoch.obsType    = typesObs;
          receiverArc.push_back(epoch);
        }

        // simulate tracks
        // ---------------
        for(UInt idTrans=0; idTrans<gnss.transmitter.size(); idTrans++)
          if(gnss.transmitter.at(idTrans)->useable())
          {
            // find types of this system
            std::vector<GnssType> typesTrans;
            for(UInt idType=0; idType<typesObs.size(); idType++)
              if(typesObs.at(idType) == (gnss.transmitter.at(idTrans)->PRN() & GnssType::SYSTEM))
                typesTrans.push_back(typesObs.at(idType));

            std::vector<UInt>   epochList;
            std::vector<Vector> obsList, sigmaList;
            Vector ambiguities(typesTrans.size());
            Double phaseWindup = 0;
            Bool   phaseBreak  = TRUE;
            for(UInt idEpoch=0; idEpoch<gnss.times.size(); idEpoch++)
              if(gnss.transmitter.at(idTrans)->useable())
              {
                Vector obs, sigma;
                if(!computeObservation(*gnss.receiver.at(idRecv), *gnss.transmitter.at(idTrans),
                                       gnssIonosphere, idEpoch, typesTrans, phaseWindup, obs, sigma))
                {
                  phaseBreak = TRUE;
                  continue;
                }

                // new integer ambiguities after phase break
                if(phaseBreak)
                  for(UInt i=0; i<ambiguities.rows(); i++)
                    ambiguities(i) = ambiguityRandom(generator);
                phaseBreak = FALSE;

                // add integer ambiguities to phase observations
                for(UInt idType=0; idType<typesTrans.size(); idType++)
                  if(typesTrans.at(idType) == GnssType::PHASE)
                    obs(idType) += ambiguities(idType) * LIGHT_VELOCITY/typesTrans.at(idType).frequency();

                epochList.push_back(idEpoch);
                obsList.push_back(obs);
                sigmaList.push_back(sigma);
              } // for(idEpoch)

            // add noise
            Matrix eps = noiseObs->noise(obsList.size(), typesTrans.size());
            for(UInt i=0; i<obsList.size(); i++)
              for(UInt idType=0; idType<typesTrans.size(); idType++)
                obsList.at(i)(idType) += sigmaList.at(i)(idType) * eps(i, idType);

            // append to receiver epochs
            for(UInt i=0; i<epochList.size(); i++)
            {
              const UInt idEpoch = epochList.at(i);
              receiverArc.at(idEpoch).satellite.push_back( gnss.transmitter.at(idTrans)->PRN() );
              for(UInt idType=0; idType<typesTrans.size(); idType++)
                receiverArc.at(idEpoch).observation.push_back( obsList.at(i)(idType) );
            }
          } // for(idTrans)

      // remove empty epochs
      GnssReceiverArc receiverArc2;
      for(UInt i=0; i<receiverArc.size(); i++)
        if(receiverArc.at(i).satellite.size())
          receiverArc2.push_back(receiverArc.at(i));

      // write simulated receiver
      logStatus<<"write gnss data to <"<<fileNameReceiver(fileNameVariableList)<<">"<<Log::endl;
      logInfo<<"  number of epochs: "<<receiverArc2.size()<<Log::endl;
      InstrumentFile::write(fileNameReceiver(fileNameVariableList), receiverArc2);
    } // for(idRecv)
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Bool GnssSimulateReceiver::computeObservation(const Gnss::Receiver &receiver, const Gnss::Transmitter &transmitter, GnssParametrizationIonospherePtr /*ionosphere*/,
                                              UInt idEpoch, const std::vector<GnssType> &types, Double &phaseWindupOld,
                                              Vector &obs, Vector &sigma)
{
  try
  {
    // position, time of transmitter & receiver
    // ----------------------------------------
    Time     timeTrans, timeRecv = receiver.timeCorrected(idEpoch);
    Vector3d posTrans,  posRecv  = receiver.position(idEpoch);
    transmitter.transmitTime(idEpoch, timeRecv, posRecv, timeTrans, posTrans);

    // orientation of antennas
    // -----------------------
    const Vector3d k                  = normalize(posRecv - posTrans);                              // line of sight from transmitter to receiver
    const Vector3d kRecvLocal         = receiver.celestial2localFrame(idEpoch).transform(-k);       // line of sight in local frame (north, east, up or vehicle)
    const Angle    azimutRecvLocal    = kRecvLocal.lambda();
    const Angle    elevationRecvLocal = kRecvLocal.phi();
    const Vector3d kRecvAnt           = receiver.local2antennaFrame(idEpoch).transform(kRecvLocal); // line of sight in left-handed antenna system (useally north, east, up)
    const Angle    azimutRecvAnt      = kRecvAnt.lambda();
    const Angle    elevationRecvAnt   = kRecvAnt.phi();

    if(elevationRecvAnt < receiver.elevationCutOff())
      return FALSE;

    const Vector3d kTrans             = transmitter.celestial2antennaFrame(idEpoch, timeTrans).transform(k);
    const Angle    azimutTrans        = kTrans.lambda();
    const Angle    elevationTrans     = kTrans.phi();

    // Observed range
    // --------------
    Double r12;
    const Double r1 = posTrans.r();
    const Double r2 = posRecv.r();
    r12  = (posRecv - posTrans).r();
    r12 += 2*DEFAULT_GM/pow(LIGHT_VELOCITY,2)*log((r1+r2+r12)/(r1+r2-r12));            // curved space-time
    r12 += 2*inner(posTrans, transmitter.velocity(idEpoch, timeTrans))/LIGHT_VELOCITY; // relativistic clock correction
    r12 += receiver.troposphere(idEpoch, azimutRecvLocal, elevationRecvLocal);
    r12 -= LIGHT_VELOCITY * transmitter.clockError(idEpoch, timeTrans);
    r12 += LIGHT_VELOCITY * receiver.clockError(idEpoch);

    obs = Vector(types.size());
    for(UInt i=0; i<types.size(); i++)
      if((types.at(i) == GnssType::RANGE) || (types.at(i) == GnssType::PHASE))
        obs(i) += r12;

    // Composed signals (e.g. C2DG)
    std::vector<GnssType> typesTransmitted;
    Matrix T;
    receiver.signalComposition(idEpoch, types, typesTransmitted, T);

    // antenna correction
    // ------------------
    obs += receiver.antennaVariations(idEpoch, azimutRecvAnt,  elevationRecvAnt,  types);
    obs += T * transmitter.antennaVariations(idEpoch, azimutTrans, elevationTrans, typesTransmitted);

    // ionospheric effects
    // -------------------
    // TODO

    // phase wind-up
    // Carrier phase wind-up in GPS reflectometry, Georg Beyerle, Springer Verlag 2008
    // -------------------------------------------------------------------------------
    const Transform3d crf2arfRecv  = receiver.local2antennaFrame(idEpoch) * receiver.celestial2localFrame(idEpoch);
    const Transform3d crf2arfTrans = transmitter.celestial2antennaFrame(idEpoch, timeTrans);
    const Vector3d Tx = crf2arfRecv.transform(crossProduct(crossProduct(k, crf2arfTrans.inverseTransform(Vector3d(1,0,0))), k));
    const Vector3d Ty = crf2arfRecv.transform(crossProduct(crossProduct(k, crf2arfTrans.inverseTransform(Vector3d(0,1,0))), k));
    Double phaseWindup = atan2(Tx.y()+Ty.x(), Ty.y()-Tx.x()); // both left-handed systems
    while((phaseWindupOld-phaseWindup) > PI)
      phaseWindup += 2*PI;
    while((phaseWindupOld-phaseWindup) < -PI)
      phaseWindup -= 2*PI;
    phaseWindupOld = phaseWindup;
    for(UInt i=0; i<types.size(); i++)
      if(types.at(i) == GnssType::PHASE)
        obs(i) += phaseWindup/(2*PI) * (LIGHT_VELOCITY/types.at(i).frequency());

    // add observation noise
    // ---------------------
    sigma = receiver.accuracy(idEpoch, azimutRecvAnt, elevationRecvAnt, types);

    return TRUE;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
