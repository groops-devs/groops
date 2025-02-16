/***********************************************/
/**
* @file gnss.cpp
*
* @brief global navigation satellite system.
*
* @author Torsten Mayer-Guerr
* @date 2010-08-03
*
*/
/***********************************************/

#include "base/import.h"
#include "base/planets.h"
#include "config/config.h"
#include "inputOutput/logging.h"
#include "classes/earthRotation/earthRotation.h"
#include "classes/platformSelector/platformSelector.h"
#include "gnss.h"
#include "gnss/gnssObservation.h"
#include "gnss/gnssDesignMatrix.h"
#include "gnss/gnssTransmitter.h"
#include "gnss/gnssReceiver.h"
#include "gnss/gnssParametrization/gnssParametrization.h"
#include "gnss/gnssTransmitterGenerator/gnssTransmitterGenerator.h"
#include "gnss/gnssReceiverGenerator/gnssReceiverGenerator.h"

/***********************************************/

void Gnss::init(std::vector<GnssType> simulationTypes, const std::vector<Time> &times, const Time &timeMargin,
                GnssTransmitterGeneratorPtr transmitterGenerator, GnssReceiverGeneratorPtr receiverGenerator,
                EarthRotationPtr earthRotation, GnssParametrizationPtr parametrization, Parallel::CommunicatorPtr comm)
{
  try
  {
    this->times = times;

    // init earth rotation
    // -------------------
    eop = Matrix(times.size(), 8); // Matrix eop columns: xp, yp, sp, deltaUT, LOD, X, Y, S
    for(UInt i=0; i<times.size(); i++)
      earthRotation->earthOrientationParameter(times.at(i), eop(i,0), eop(i,1), eop(i,2), eop(i,3), eop(i,4), eop(i,5), eop(i,6), eop(i,7));
    // UT1-UTC => UT1-GPS (avoid leap seconds jumps for interpolation)
    for(UInt i=0; i<times.size(); i++)
      eop(i,3) -= (times.at(i)-timeGPS2UTC(times.at(i))).seconds();

    funcRotationCrf2Trf = std::bind(&Gnss::rotationCrf2Trf, this, std::placeholders::_1);

    // init transmitters
    // -----------------
    transmitters = transmitterGenerator->transmitters(times);
    for(UInt idTrans=0; idTrans<transmitters.size(); idTrans++)
      transmitters.at(idTrans)->id_ = idTrans;

    // init receivers
    // --------------
    receivers = receiverGenerator->receivers(simulationTypes, times, timeMargin, transmitters, earthRotation, comm);
    for(UInt idRecv=0; idRecv<receivers.size(); idRecv++)
      receivers.at(idRecv)->id_ = idRecv;
    synchronizeTransceivers(comm);

    // init parametrization
    // --------------------
    this->parametrization = parametrization;
    if(parametrization)
    {
      parametrization->init(this, comm);
      funcReduceModels = std::bind(&GnssParametrization::observationCorrections, parametrization, std::placeholders::_1);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Rotary3d Gnss::rotationCrf2Trf(const Time &time) const
{
  try
  {
    const UInt idx = std::min((times.size()-1), static_cast<UInt>(std::distance(times.begin(),
                     std::upper_bound(times.begin(), times.end(), time, [](const Time &t, const Time &s) {return (t-s).seconds() < 0.5;}))));
    const Double xp      = eop(idx, 0);
    const Double yp      = eop(idx, 1);
    const Double sp      = eop(idx, 2);
    const Double deltaUT = eop(idx, 3) + (time-timeGPS2UTC(time)).seconds();
    const Double X       = eop(idx, 5);
    const Double Y       = eop(idx, 6);
    const Double S       = eop(idx, 7);

    const Double ERA = Planets::ERA(timeGPS2UTC(time) + seconds2time(deltaUT));
    const Double r2  = X*X + Y*Y;
    const Double E   = (r2!=0.) ? std::atan2(Y, X) : 0.;
    const Double D   = std::atan(std::sqrt(r2/(1-r2)));

    return  rotaryX(Angle(-yp)) * rotaryY(Angle(-xp)) *
            rotaryZ(Angle(sp+ERA-S-E)) *
            rotaryY(Angle(D)) * rotaryZ(Angle(E));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void Gnss::synchronizeTransceivers(Parallel::CommunicatorPtr comm)
{
  try
  {
    // distribute process id of receivers
    // ----------------------------------
    Vector recvProcess(receivers.size());
    for(UInt idRecv=0; idRecv<receivers.size(); idRecv++)
      if(receivers.at(idRecv)->isMyRank())
        recvProcess(idRecv) = Parallel::myRank(comm)+1;
    Parallel::reduceSum(recvProcess, 0, comm);
    Parallel::broadCast(recvProcess, 0, comm);

    // synchronize transceivers
    // ------------------------
    for(UInt idRecv=0; idRecv<receivers.size(); idRecv++)
      if(recvProcess(idRecv))
        Parallel::broadCast(static_cast<GnssTransceiver&>(*receivers.at(idRecv)), static_cast<UInt>(recvProcess(idRecv)-1), comm);
      else if(receivers.at(idRecv)->useable())
        receivers.at(idRecv)->disable("");

    // collect observation types
    // -------------------------
    typesRecvTrans.clear();
    typesRecvTrans.resize(receivers.size(), std::vector<std::vector<GnssType>>(transmitters.size()));
    for(auto recv : receivers)
    {
      if(recv->isMyRank())
        for(UInt idTrans=0; idTrans<transmitters.size(); idTrans++)
        {
          for(UInt idEpoch=0; idEpoch<recv->idEpochSize(); idEpoch++)
          {
            auto obs = recv->observation(idTrans, idEpoch);
            if(obs)
              for(UInt idType=0; idType<obs->size(); idType++)
                if(!obs->at(idType).type.isInList(typesRecvTrans.at(recv->idRecv()).at(idTrans)))
                  typesRecvTrans.at(recv->idRecv()).at(idTrans).push_back(obs->at(idType).type);
          }
          std::sort(typesRecvTrans.at(recv->idRecv()).at(idTrans).begin(), typesRecvTrans.at(recv->idRecv()).at(idTrans).end());
        }
      if(recv->useable())
        Parallel::broadCast(typesRecvTrans.at(recv->idRecv()), static_cast<UInt>(recvProcess(recv->idRecv())-1), comm);
    }

    // adjust signal biases to available observation types
    // ---------------------------------------------------
    for(auto trans : transmitters)
    {
      std::vector<GnssType> types;
      for(auto &typesTrans : typesRecvTrans)
        for(GnssType type : typesTrans.at(trans->idTrans()))
          if((type == GnssType::PHASE) || (type == GnssType::RANGE))
            if(!type.isInList(types))
              types.push_back(type + trans->PRN());
      types = GnssType::replaceCompositeSignals(types);

      if(types.size())
      {
        trans->signalBias.biases = trans->signalBias.compute(types); // apriori signal bias
        trans->signalBias.types  = types;
      }
    }

    for(auto recv : receivers)
    {
      std::vector<GnssType> types;
      for(auto &typesTrans : typesRecvTrans.at(recv->idRecv()))
        for(GnssType type : typesTrans)
          if((type == GnssType::PHASE) || (type == GnssType::RANGE))
            if(!type.isInList(types))
              types.push_back(type & ~GnssType::PRN);
      std::sort(types.begin(), types.end());

      if(types.size())
      {
        recv->signalBias.biases = recv->signalBias.compute(types); // apriori signal bias
        recv->signalBias.types  = types;
      }
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void Gnss::initParameter(GnssNormalEquationInfo &normalEquationInfo)
{
  try
  {
    logStatus<<"setup parameters"<<Log::endl;
    if(!parametrization)
      throw(Exception("no parametrization given"));

    // disable unuseable transmitters/receivers/epochs
    // -----------------------------------------------
    // check number of required observations
    std::vector<UInt> transCount(transmitters.size(), 0), transCountEpoch(transmitters.size(), 0);
    std::vector<UInt> recvCount(receivers.size(), 0),     recvCountEpoch(receivers.size(), 0);
    parametrization->requirements(normalEquationInfo, transCount, transCountEpoch, recvCount, recvCountEpoch);
    UInt disabledEpochsTrans = 0;
    UInt disabledEpochsRecv = 0;
    for(;;)
    {
      UInt mustSync = 0;

      // disable transmitters
      for(auto trans : transmitters)
        if(trans->useable() && (transCount.at(trans->idTrans()) || transCountEpoch.at(trans->idTrans())))
        {
          Vector countEpoch(times.size());
          for(const auto &recv : receivers)
            if(normalEquationInfo.estimateReceiver.at(recv->idRecv()) && recv->isMyRank())
              for(UInt idEpoch : normalEquationInfo.idEpochs)
                if(trans->useable(idEpoch) && recv->useable(idEpoch) && recv->observation(trans->idTrans(), idEpoch))
                  countEpoch(idEpoch)++;
          Parallel::reduceSum(countEpoch, 0, normalEquationInfo.comm);
          Parallel::broadCast(countEpoch, 0, normalEquationInfo.comm);

          for(UInt idEpoch : normalEquationInfo.idEpochs)
            if(trans->useable(idEpoch) && (countEpoch(idEpoch) < transCountEpoch.at(trans->idTrans())))
            {
              // logWarningOnce<<trans->name()<<" disabled epoch "<<times.at(idEpoch).dateTimeStr()<<Log::endl;
              disabledEpochsTrans++;
              trans->disable(idEpoch, "failed parametrization requirements");
              mustSync = TRUE;
              for(const auto &recv : receivers)
                if(recv->isMyRank() && recv->observation(trans->idTrans(), idEpoch))
                  recv->deleteObservation(trans->idTrans(), idEpoch);
            }

          UInt epochCount = 0;
          for(UInt idEpoch : normalEquationInfo.idEpochs)
            if(countEpoch(idEpoch) > transCountEpoch.at(trans->idTrans()))
              epochCount++;
          if(epochCount < transCount.at(trans->idTrans()))
          {
            logWarningOnce<<trans->name()<<" disabled: not enough estimable epochs ("<<epochCount<<")"<<Log::endl;
            trans->disable("not enough estimable epochs ("+epochCount%"%i)"s);
            mustSync = TRUE;
            for(UInt idEpoch=0; idEpoch<times.size(); idEpoch++)
              for(const auto &recv : receivers)
                if(recv->isMyRank() && recv->observation(trans->idTrans(), idEpoch))
                  recv->deleteObservation(trans->idTrans(), idEpoch);
          }
        }

      // disable receivers
      for(auto recv : receivers)
        if(normalEquationInfo.estimateReceiver.at(recv->idRecv()) && recv->isMyRank() && (recvCount.at(recv->idRecv()) || recvCountEpoch.at(recv->idRecv())))
        {
          UInt epochCount = 0;
          for(UInt idEpoch : normalEquationInfo.idEpochs)
            if(recv->useable(idEpoch))
            {
              UInt count = 0;
              for(auto trans : transmitters)
                if(trans->useable(idEpoch) && recv->observation(trans->idTrans(), idEpoch))
                  count++;
              if(count < recvCountEpoch.at(recv->idRecv()))
              {
                // logWarning<<recv->name()<<" disabled epoch "<<times.at(idEpoch).dateTimeStr()<<Log::endl;
                disabledEpochsRecv++;
                recv->disable(idEpoch, "failed parametrization requirements");
                mustSync = TRUE;
              }
              else if(count > recvCountEpoch.at(recv->idRecv()))
                epochCount++;
            }

          if(epochCount < recvCount.at(recv->idRecv()))
          {
            logWarning<<recv->name()<<" disabled: not enough estimable epochs ("<<epochCount<<")"<<Log::endl;
            recv->disable("not enough estimable epochs ("+epochCount%"%i)"s);
            mustSync = TRUE;
          }
        }

      // is something disabled?
      Parallel::reduceSum(mustSync, 0, normalEquationInfo.comm);
      Parallel::broadCast(mustSync, 0, normalEquationInfo.comm);
      if(!mustSync)
        break;

      // synchronize transceivers
      synchronizeTransceivers(normalEquationInfo.comm);
      Parallel::reduceSum(disabledEpochsRecv, 0, normalEquationInfo.comm);
      if(!Parallel::isMaster(normalEquationInfo.comm))
        disabledEpochsRecv = 0;
    } // for(;;)

    if(disabledEpochsTrans) logWarningOnce<<disabledEpochsTrans<<" disabled transmitter epochs"<<Log::endl;
    if(disabledEpochsRecv)  logWarningOnce<<disabledEpochsRecv <<" disabled receiver epochs"<<Log::endl;

    for(auto recv : receivers)
      if(!recv->useable())
        normalEquationInfo.estimateReceiver.at(recv->idRecv()) = FALSE;

    // distribute process id of receivers
    // ----------------------------------
    Vector recvProcess(receivers.size());
    for(UInt idRecv=0; idRecv<receivers.size(); idRecv++)
      if(receivers.at(idRecv)->isMyRank())
        recvProcess(idRecv) = Parallel::myRank(normalEquationInfo.comm)+1;
    Parallel::reduceSum(recvProcess, 0, normalEquationInfo.comm);
    Parallel::broadCast(recvProcess, 0, normalEquationInfo.comm);

    // init parameters
    // ---------------
    normalEquationInfo.initNewParameterNames();
    parametrization->initParameter(normalEquationInfo);
    normalEquationInfo.calculateIndex(recvProcess);
    logInfo<<"+ ======="<<Log::endl;
    logInfo<<normalEquationInfo.parameterCount()%"%9i parameters in "s<<normalEquationInfo.blockCount()<<" normal equation matrix blocks"<<Log::endl;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Vector Gnss::aprioriParameter(const GnssNormalEquationInfo &normalEquationInfo) const
{
  try
  {
    Vector x0 = parametrization->aprioriParameter(normalEquationInfo);
    Parallel::reduceSum(x0, 0, normalEquationInfo.comm);
    return x0;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Bool Gnss::basicObservationEquations(const GnssNormalEquationInfo &/*normalEquationInfo*/, UInt idRecv, UInt idTrans, UInt idEpoch, GnssObservationEquation &eqn) const
{
  try
  {
    std::vector<GnssType> type;
    if(!receivers.at(idRecv)->observation(idTrans, idEpoch) ||
       !receivers.at(idRecv)->observation(idTrans, idEpoch)->observationList(GnssObservation::RANGE | GnssObservation::PHASE, type))
      return FALSE;
    eqn.compute(*receivers.at(idRecv)->observation(idTrans, idEpoch), *receivers.at(idRecv), *transmitters.at(idTrans),
                funcRotationCrf2Trf, funcReduceModels, idEpoch, TRUE, type);
    return TRUE;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void Gnss::designMatrix(const GnssNormalEquationInfo &normalEquationInfo, const GnssObservationEquation &eqn, GnssDesignMatrix &A) const
{
  try
  {
    if(eqn.l.rows())
      parametrization->designMatrix(normalEquationInfo, eqn, A);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void Gnss::constraintsEpoch(const GnssNormalEquationInfo &normalEquationInfo, UInt idEpoch, MatrixDistributed &normals, std::vector<Matrix> &n, Double &lPl, UInt &obsCount) const
{
  try
  {
    parametrization->constraintsEpoch(normalEquationInfo, idEpoch, normals, n, lPl, obsCount);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void Gnss::constraints(const GnssNormalEquationInfo &normalEquationInfo, MatrixDistributed &normals, std::vector<Matrix> &n, Double &lPl, UInt &obsCount) const
{
  try
  {
    parametrization->constraints(normalEquationInfo, normals, n, lPl, obsCount);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Double Gnss::ambiguityResolve(const GnssNormalEquationInfo &normalEquationInfo, MatrixDistributed &normals, std::vector<Matrix> &n, Double &lPl, UInt &obsCount,
                              const std::vector<Byte> &selectedTransmitters, const std::vector<Byte> &selectedReceivers,
                              const std::function<Vector(const_MatrixSliceRef xFloat, MatrixSliceRef W, const_MatrixSliceRef d, Vector &xInt, Double &sigma)> &searchInteger)
{
  try
  {
    return parametrization->ambiguityResolve(normalEquationInfo, normals, n, lPl, obsCount,
                                             selectedTransmitters, selectedReceivers, searchInteger);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Double Gnss::updateParameter(const GnssNormalEquationInfo &normalEquationInfo, const_MatrixSliceRef x, const_MatrixSliceRef Wz)
{
  try
  {
    return parametrization->updateParameter(normalEquationInfo, x, Wz);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void Gnss::updateCovariance(const GnssNormalEquationInfo &normalEquationInfo, const MatrixDistributed &covariance)
{
  try
  {
    parametrization->updateCovariance(normalEquationInfo, covariance);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void Gnss::writeResults(const GnssNormalEquationInfo &normalEquationInfo, const std::string &suffix)
{
  try
  {
    parametrization->writeResults(normalEquationInfo, suffix);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/

std::vector<GnssType> Gnss::types(const GnssType mask) const
{
  try
  {
    std::vector<GnssType> types;
    for(UInt idRecv=0; idRecv<receivers.size(); idRecv++)
      for(UInt idTrans=0; idTrans<transmitters.size(); idTrans++)
        for(GnssType type : typesRecvTrans.at(idRecv).at(idTrans))
          if(!type.isInList(types))
             types.push_back(type & mask);
    std::sort(types.begin(), types.end());
    return types;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/

std::vector<Byte> Gnss::selectTransmitters(PlatformSelectorPtr selector)
{
  try
  {
    std::vector<const Platform*> platforms(transmitters.size(), nullptr);
    for(UInt idTrans=0; idTrans<transmitters.size(); idTrans++)
      if(transmitters.at(idTrans)->useable())
        platforms.at(idTrans) = &transmitters.at(idTrans)->platform;
    return selector->select(times.front(), times.back(), platforms);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

std::vector<Byte> Gnss::selectReceivers(PlatformSelectorPtr selector)
{
  try
  {
    std::vector<const Platform*> platforms(receivers.size(), nullptr);
    for(UInt idRecv=0; idRecv<receivers.size(); idRecv++)
      if(receivers.at(idRecv)->useable())
        platforms.at(idRecv) = &receivers.at(idRecv)->platform;
    return selector->select(times.front(), times.back(), platforms);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/

Bool Gnss::InfoParameterChange::update(Double change)
{
  try
  {
    count++;
    rms += change*change;
    if(std::fabs(change) > std::fabs(maxChange))
    {
      maxChange = change;
      return TRUE;
    }
    return FALSE;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void Gnss::InfoParameterChange::synchronizeAndPrint(Parallel::CommunicatorPtr comm, Double convertToMeter, Double &maxChangeTotal)
{
  try
  {
    Vector change(Parallel::size(comm));
    change(Parallel::myRank(comm)) = maxChange;
    Parallel::reduceSum(change, 0, comm);
    Parallel::broadCast(change, 0, comm);
    UInt idProcess = 0;
    maxChange = 0;
    for(UInt i=0; i<change.rows(); i++)
      if(std::fabs(change(i)) > std::fabs(maxChange))
      {
        maxChange = change(i);
        idProcess = i;
      }
    maxChangeTotal = std::max(maxChangeTotal, std::fabs(convertToMeter*maxChange));

    Parallel::reduceSum(count, 0, comm);
    Parallel::reduceSum(rms,   0, comm);

    if((idProcess != 0) && (idProcess == Parallel::myRank(comm)))
      Parallel::send(info, 0, comm);
    else if(Parallel::isMaster(comm))
    {
      if(idProcess != 0)
        Parallel::receive(info, idProcess, comm);
      rms = std::sqrt(rms/count);
      if(!info.empty())
      {
        std::string space(5-std::min(unit.size(), std::size_t(4)), ' ');
        logInfo<<"  rms ="<<rms%"%7.1f "s<<unit<<","<<space<<"max ="<<maxChange%"%8.1f "s<<unit<<","<<space<<info<<Log::endl;
      }
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
