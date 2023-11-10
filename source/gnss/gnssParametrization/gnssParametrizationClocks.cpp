/***********************************************/
/**
* @file gnssParametrizationClocks.h
*
* @brief Clocks errors.
* @see GnssParametrization
*
* @author Torsten Mayer-Guerr
* @date 2021-01-23
*
*/
/***********************************************/

#include "base/import.h"
#include "config/config.h"
#include "files/fileInstrument.h"
#include "classes/platformSelector/platformSelector.h"
#include "gnss/gnssParametrization/gnssParametrization.h"
#include "gnss/gnssParametrization/gnssParametrizationClocks.h"

/***********************************************/

GnssParametrizationClocks::GnssParametrizationClocks(Config &config)
{
  try
  {
    readConfig(config, "name",                       name,                       Config::OPTIONAL, "parameter.clocks", "used for parameter selection");
    readConfig(config, "selectTransmitters",         selectTransmitters,         Config::DEFAULT,  R"(["all"])", "");
    readConfig(config, "selectReceivers",            selectReceivers,            Config::DEFAULT,  R"(["all"])", "");
    readConfig(config, "outputfileClockTransmitter", fileNameTransmitter,        Config::OPTIONAL, "", "variable {prn} available");
    readConfig(config, "outputfileClockReceiver",    fileNameReceiver,           Config::OPTIONAL, "", "variable {station} available");
    readConfig(config, "nameConstraint",             nameConstraint,             Config::OPTIONAL, "constraintEpoch.clocks", "used for parameter selection");
    readConfig(config, "selectTransmittersZeroMean", selectTransmittersZeroMean, Config::DEFAULT,  R"(["all"])", "");
    readConfig(config, "selectReceiversZeroMean",    selectReceiversZeroMean,    Config::DEFAULT,  "", "");
    readConfig(config, "sigmaZeroMeanConstraint",    sigmaZeroMean,              Config::DEFAULT,  "0.0001", "(0 = unconstrained) sigma [m] for zero-mean constraint over all selected clocks");
    if(isCreateSchema(config)) return;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssParametrizationClocks::init(Gnss *gnss, Parallel::CommunicatorPtr comm)
{
  try
  {
    this->gnss                   = gnss;
    selectedTransmitters         = gnss->selectTransmitters(selectTransmitters);
    selectedReceivers            = gnss->selectReceivers(selectReceivers);
    selectedTransmittersZeroMean = gnss->selectTransmitters(selectTransmittersZeroMean);
    selectedReceiversZeroMean    = gnss->selectReceivers(selectReceiversZeroMean);

    x0Trans.resize(gnss->transmitters.size());
    if(Parallel::isMaster(comm))
      for(auto trans : gnss->transmitters)
      {
        const UInt idTrans = trans->idTrans();
        if(trans->useable() && selectedTransmittersZeroMean.at(idTrans))
        {
          x0Trans.at(idTrans).resize(gnss->times.size(), 0);
          for(UInt idEpoch=0; idEpoch<gnss->times.size(); idEpoch++)
            if(trans->useable(idEpoch))
              x0Trans.at(idTrans).at(idEpoch) = LIGHT_VELOCITY * trans->clockError(idEpoch);
        }
      }

    x0Recv.resize(gnss->receivers.size());
    for(auto recv : gnss->receivers)
    {
      const UInt idRecv = recv->idRecv();
      if(recv->isMyRank() && selectedReceiversZeroMean.at(idRecv))
      {
        x0Recv.at(idRecv).resize(gnss->times.size(), 0);
        for(UInt idEpoch=0; idEpoch<gnss->times.size(); idEpoch++)
          if(recv->useable(idEpoch))
            x0Recv.at(idRecv).at(idEpoch) = LIGHT_VELOCITY * recv->clockError(idEpoch);
      }
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssParametrizationClocks::requirements(GnssNormalEquationInfo &normalEquationInfo,
                                             std::vector<UInt> &/*transCount*/, std::vector<UInt> &transCountEpoch,
                                             std::vector<UInt> &/*recvCount*/,  std::vector<UInt> &recvCountEpoch)
{
  try
  {
    if(!isEnabled(normalEquationInfo, name))
      return;

    // transmitter epochs
    if(!normalEquationInfo.isEachReceiverSeparately)
      for(auto trans : gnss->transmitters)
        if(trans->useable() && selectedTransmitters.at(trans->idTrans()))
          transCountEpoch.at(trans->idTrans())++;

    // receiver clocks
    for(auto recv : gnss->receivers)
      if(recv->useable() && normalEquationInfo.estimateReceiver.at(recv->idRecv()) && selectedReceivers.at(recv->idRecv()))
        recvCountEpoch.at(recv->idRecv())++;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssParametrizationClocks::initParameter(GnssNormalEquationInfo &normalEquationInfo)
{
  try
  {
    indexTrans.clear();
    indexTrans.resize(gnss->transmitters.size());
    indexRecv.clear();
    indexRecv.resize(gnss->receivers.size());
    applyConstraint = FALSE;
    if(!isEnabled(normalEquationInfo, name))
      return;

    UInt countParaTrans = 0;
    if(!normalEquationInfo.isEachReceiverSeparately)
      for(auto trans : gnss->transmitters)
      {
        const UInt idTrans = trans->idTrans();
        if(trans->useable() && selectedTransmitters.at(idTrans))
        {
          indexTrans.at(idTrans).resize(gnss->times.size());
          for(UInt idEpoch : normalEquationInfo.idEpochs)
            if(trans->useable(idEpoch))
            {
              indexTrans.at(idTrans).at(idEpoch) = normalEquationInfo.parameterNamesEpochTransmitter(idEpoch, idTrans, {ParameterName(trans->name(), "clock", "", gnss->times.at(idEpoch))});
              countParaTrans++;
            }
        }
      }
    if(countParaTrans)
      logInfo<<countParaTrans%"%9i transmitter clock error parameters"s<<Log::endl;

    UInt countParaRecv = 0;
    for(auto recv : gnss->receivers)
    {
      const UInt idRecv = recv->idRecv();
      if(recv->useable() && normalEquationInfo.estimateReceiver.at(idRecv) && selectedReceivers.at(idRecv))
      {
        indexRecv.at(idRecv).resize(gnss->times.size());
        for(UInt idEpoch : normalEquationInfo.idEpochs)
          if(recv->useable(idEpoch))
          {
            indexRecv.at(idRecv).at(idEpoch) = normalEquationInfo.parameterNamesEpochReceiver(idEpoch, idRecv, {ParameterName(recv->name(), "clock", "", gnss->times.at(idEpoch))});
            countParaRecv++;
          }
      }
    }
    if(countParaRecv)
      logInfo<<countParaRecv%"%9i receiver clock error parameters"s<<Log::endl;

    applyConstraint = isEnabled(normalEquationInfo, nameConstraint) && sigmaZeroMean
                    && (countParaTrans+countParaRecv) &&  !normalEquationInfo.isEachReceiverSeparately;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssParametrizationClocks::aprioriParameter(const GnssNormalEquationInfo &normalEquationInfo, MatrixSliceRef x0) const
{
  try
  {
    if(Parallel::isMaster(normalEquationInfo.comm))
      for(auto trans : gnss->transmitters)
      {
        const UInt idTrans = trans->idTrans();
        for(UInt idEpoch : normalEquationInfo.idEpochs)
          if(indexTrans.at(idTrans).size() && indexTrans.at(idTrans).at(idEpoch))
            x0(normalEquationInfo.index(indexTrans.at(idTrans).at(idEpoch)), 0) = LIGHT_VELOCITY * trans->clockError(idEpoch);
      }

    for(auto recv : gnss->receivers)
      if(recv->isMyRank())
      {
        const UInt idRecv = recv->idRecv();
        for(UInt idEpoch : normalEquationInfo.idEpochs)
          if(indexRecv.at(idRecv).size() && indexRecv.at(idRecv).at(idEpoch))
            x0(normalEquationInfo.index(indexRecv.at(idRecv).at(idEpoch)), 0) = LIGHT_VELOCITY * recv->clockError(idEpoch);
      }
    }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssParametrizationClocks::designMatrix(const GnssNormalEquationInfo &/*normalEquationInfo*/, const GnssObservationEquation &eqn, GnssDesignMatrix &A) const
{
  try
  {
    // transmitter clock
    if(indexTrans.at(eqn.transmitter->idTrans()).size() && indexTrans.at(eqn.transmitter->idTrans()).at(eqn.idEpoch))
      copy(eqn.A.column(GnssObservationEquation::idxClockTrans,1), A.column(indexTrans.at(eqn.transmitter->idTrans()).at(eqn.idEpoch)));

    // receiver clock
    if(indexRecv.at(eqn.receiver->idRecv()).size() && indexRecv.at(eqn.receiver->idRecv()).at(eqn.idEpoch))
      copy(eqn.A.column(GnssObservationEquation::idxClockRecv,1), A.column(indexRecv.at(eqn.receiver->idRecv()).at(eqn.idEpoch)));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssParametrizationClocks::constraintsEpoch(const GnssNormalEquationInfo &normalEquationInfo, UInt idEpoch, MatrixDistributed &normals, std::vector<Matrix> &n, Double &lPl, UInt &obsCount) const
{
  try
  {
    if(!applyConstraint)
      return;

    // zero-mean constraint of clocks
    // ------------------------------
    UInt countTrans = 0, countRecv = 0;
    for(UInt idTrans=0; idTrans<indexTrans.size(); idTrans++)
      if(selectedTransmittersZeroMean.at(idTrans) && indexTrans.at(idTrans).size() && indexTrans.at(idTrans).at(idEpoch))
        countTrans++;
    for(UInt idRecv=0; idRecv<indexRecv.size(); idRecv++)
      if(selectedReceiversZeroMean.at(idRecv) && indexRecv.at(idRecv).size() && indexRecv.at(idRecv).at(idEpoch))
        countRecv++;
    UInt count = countTrans + countRecv;
    if(!count)
      return;

    // collect apriori clock errors
    Vector lRecvEpoch;
    if(countRecv)
    {
      lRecvEpoch = Vector(gnss->receivers.size());
      for(UInt idRecv=0; idRecv<indexRecv.size(); idRecv++)
        if(selectedReceiversZeroMean.at(idRecv) && indexRecv.at(idRecv).size() && indexRecv.at(idRecv).at(idEpoch) && gnss->receivers.at(idRecv)->isMyRank())
          lRecvEpoch(idRecv) = x0Recv.at(idRecv).at(idEpoch) - LIGHT_VELOCITY*gnss->receivers.at(idRecv)->clockError(idEpoch);
      Parallel::reduceSum(lRecvEpoch, 0, normalEquationInfo.comm);
    }

    if(Parallel::isMaster(normalEquationInfo.comm))
    {
      GnssDesignMatrix A(normalEquationInfo, Vector(1));
      for(UInt idTrans=0; idTrans<indexTrans.size(); idTrans++)
        if(selectedTransmittersZeroMean.at(idTrans) && indexTrans.at(idTrans).size() && indexTrans.at(idTrans).at(idEpoch))
        {
          A.l(0) += (x0Trans.at(idTrans).at(idEpoch) - LIGHT_VELOCITY*gnss->transmitters.at(idTrans)->clockError(idEpoch))/count/sigmaZeroMean; // remove apriori value -> regularization towards 0
          A.column(indexTrans.at(idTrans).at(idEpoch))(0,0) = 1./count/sigmaZeroMean;
        }
      for(UInt idRecv=0; idRecv<indexRecv.size(); idRecv++)
        if(selectedReceiversZeroMean.at(idRecv) && indexRecv.at(idRecv).size() && indexRecv.at(idRecv).at(idEpoch))
        {
          A.l(0) += lRecvEpoch(idRecv)/count/sigmaZeroMean; // remove apriori value -> regularization towards 0
          A.column(indexRecv.at(idRecv).at(idEpoch))(0,0) = 1./count/sigmaZeroMean;
        }
      A.accumulateNormals(normals, n, lPl, obsCount);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Double GnssParametrizationClocks::updateParameter(const GnssNormalEquationInfo &normalEquationInfo, const_MatrixSliceRef x, const_MatrixSliceRef /*Wz*/)
{
  try
  {
    Double maxChange = 0;
    Gnss::InfoParameterChange infoTrans("mm");
    for(auto trans : gnss->transmitters)
    {
      const UInt idTrans = trans->idTrans();
      for(UInt idEpoch : normalEquationInfo.idEpochs)
        if(indexTrans.at(idTrans).size() && indexTrans.at(idTrans).at(idEpoch))
        {
          const Double dClock = x(normalEquationInfo.index(indexTrans.at(idTrans).at(idEpoch)), 0);
          trans->updateClockError(idEpoch, dClock/LIGHT_VELOCITY);
          if(infoTrans.update(1e3*dClock))
            infoTrans.info = "clock transmitter ("+trans->name()+", "+gnss->times.at(idEpoch).dateTimeStr()+")";
        }
    }
    infoTrans.synchronizeAndPrint(normalEquationInfo.comm, 1e-3, maxChange);

    Gnss::InfoParameterChange infoRecv("mm");
    for(auto recv : gnss->receivers)
      if(recv->isMyRank())
      {
        const UInt idRecv = recv->idRecv();
        for(UInt idEpoch : normalEquationInfo.idEpochs)
          if(indexRecv.at(idRecv).size() && indexRecv.at(idRecv).at(idEpoch))
          {
            const Double dClock = x(normalEquationInfo.index(indexRecv.at(idRecv).at(idEpoch)), 0);
            recv->updateClockError(idEpoch, dClock/LIGHT_VELOCITY);
            if(infoRecv.update(1e3*dClock))
              infoRecv.info = "clock receiver ("+recv->name()+", "+gnss->times.at(idEpoch).dateTimeStr()+")";
          }
      }
    infoRecv.synchronizeAndPrint(normalEquationInfo.comm, 1e-3, maxChange);

    return maxChange;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssParametrizationClocks::writeResults(const GnssNormalEquationInfo &normalEquationInfo, const std::string &suffix) const
{
  try
  {
    if(!isEnabled(normalEquationInfo, name))
      return;

    if(!fileNameTransmitter.empty() && Parallel::isMaster(normalEquationInfo.comm))
    {
      VariableList fileNameVariableList;
      fileNameVariableList.setVariable("prn", "***");
      logStatus<<"write transmitter clocks to files <"<<fileNameTransmitter(fileNameVariableList).appendBaseName(suffix)<<">"<<Log::endl;
      for(auto trans : gnss->transmitters)
      {
        const UInt idTrans = trans->idTrans();
        MiscValueArc arc;
        for(UInt idEpoch : normalEquationInfo.idEpochs)
          if(indexTrans.at(idTrans).size() && indexTrans.at(idTrans).at(idEpoch))
          {
            MiscValueEpoch epoch;
            epoch.time  = gnss->times.at(idEpoch);
            epoch.value = trans->clockError(idEpoch);
            arc.push_back(epoch);
          }
        fileNameVariableList.setVariable("prn", trans->name());
        InstrumentFile::write(fileNameTransmitter(fileNameVariableList).appendBaseName(suffix), arc);
      }
    }

    if(!fileNameReceiver.empty())
    {
      VariableList fileNameVariableList;
      fileNameVariableList.setVariable("station", "****");
      logStatus<<"write receiver clocks to files <"<<fileNameReceiver(fileNameVariableList).appendBaseName(suffix)<<">"<<Log::endl;
      for(auto recv : gnss->receivers)
        if(recv->isMyRank())
        {
          const UInt idRecv = recv->idRecv();
          MiscValueArc arc;
          for(UInt idEpoch : normalEquationInfo.idEpochs)
            if(indexRecv.at(idRecv).size() && indexRecv.at(idRecv).at(idEpoch))
            {
              MiscValueEpoch epoch;
              epoch.time  = recv->times.at(idEpoch);
              epoch.value = recv->clockError(idEpoch);
              arc.push_back(epoch);
            }
          fileNameVariableList.setVariable("station", recv->name());
          InstrumentFile::write(fileNameReceiver(fileNameVariableList).appendBaseName(suffix), arc);
        }
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
