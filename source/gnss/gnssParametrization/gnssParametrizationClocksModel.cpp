/***********************************************/
/**
* @file gnssParametrizationClocksModel.h
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
#include "misc/varianceComponentEstimation.h"
#include "gnss/gnssTransceiverSelector/gnssTransceiverSelector.h"
#include "gnss/gnssParametrization/gnssParametrization.h"
#include "gnss/gnssParametrization/gnssParametrizationClocksModel.h"

/***********************************************/

GnssParametrizationClocksModel::GnssParametrizationClocksModel(Config &config)
{
  try
  {
    readConfig(config, "name",                       name,                       Config::OPTIONAL, "parameter.clocks", "used for parameter selection");
    readConfig(config, "selectTransmitters",         selectTransmitters,         Config::DEFAULT,  R"(["all"])", "");
    readConfig(config, "selectReceivers",            selectReceivers,            Config::DEFAULT,  R"(["all"])", "");
    readConfig(config, "outputfileClockTransmitter", fileNameTransmitter,        Config::OPTIONAL, "",    "variable {prn} available");
    readConfig(config, "outputfileClockReceiver",    fileNameReceiver,           Config::OPTIONAL, "",    "variable {station} available");
    readConfig(config, "huber",                      huber,                      Config::DEFAULT,  "2.5", "clock jumps > huber*sigma0 are downweighted");
    readConfig(config, "huberPower",                 huberPower,                 Config::DEFAULT,  "1.5", "clock jumps > huber: sigma=(e/huber)^huberPower*sigma0");
    readConfig(config, "nameConstraint",             nameConstraint,             Config::OPTIONAL, "constraintEpoch.clocksModel", "used for parameter selection");
    readConfig(config, "selectTransmittersZeroMean", selectTransmittersZeroMean, Config::DEFAULT,  R"(["all"])", "use these transmitters for zero-mean constraint");
    readConfig(config, "selectReceiversZeroMean",    selectReceiversZeroMean,    Config::DEFAULT,  "", "use these receivers for zero-mean constraint");
    readConfig(config, "sigmaZeroMeanConstraint",    sigmaZeroMean,              Config::DEFAULT,  "0.0001", "(0 = unconstrained) sigma [m] for zero-mean constraint over all selected clocks");
    if(isCreateSchema(config)) return;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssParametrizationClocksModel::init(Gnss *gnss, Parallel::CommunicatorPtr comm)
{
  try
  {
    this->gnss                   = gnss;
    selectedTransmitters         = selectTransmitters->select(gnss->transmitters);
    selectedReceivers            = selectReceivers->select(gnss->receivers);
    selectedTransmittersZeroMean = selectTransmittersZeroMean->select(gnss->transmitters);
    selectedReceiversZeroMean    = selectReceiversZeroMean->select(gnss->receivers);

    // determine apriori transmitter drift
    // -----------------------------------
    isMyRank     = Vector(gnss->transmitters.size());
    clock0Trans  = Vector(gnss->transmitters.size());
    driftTrans   = Vector(gnss->transmitters.size());
    sigma0Trans  = Vector(gnss->transmitters.size(), 0.02);
    sigmaEpochTrans.resize(gnss->transmitters.size());
    isFirstTrans.resize(gnss->transmitters.size(), TRUE);
    for(auto trans : gnss->transmitters)
    {
      const UInt idTrans = trans->idTrans();
      isMyRank(idTrans) = trans->useable() && ((idTrans % Parallel::size(comm)) == Parallel::myRank(comm));
      if(isMyRank(idTrans))
      {
        sigmaEpochTrans.at(idTrans).resize(gnss->times.size(), 1.);

        Vector l(gnss->times.size());
        Matrix A(gnss->times.size(), 2);
        for(UInt idEpoch=0; idEpoch<gnss->times.size(); idEpoch++)
          if(trans->useable(idEpoch))
          {
            l(idEpoch)    = LIGHT_VELOCITY * trans->clockError(idEpoch);
            A(idEpoch, 0) = 1;
            A(idEpoch, 1) = (gnss->times.at(idEpoch)-gnss->times.front()).seconds();
          }
        Vector x = leastSquares(Matrix(A), l);
        clock0Trans(idTrans) = x(0);
        driftTrans(idTrans)  = x(1);
      }
    }

    // determine apriori receiver drift
    // --------------------------------
    clock0Recv = Vector(gnss->receivers.size());
    driftRecv  = Vector(gnss->receivers.size());
    sigma0Recv = Vector(gnss->receivers.size(), 0.05*LIGHT_VELOCITY); // 50 ms
    sigmaEpochRecv.resize(gnss->receivers.size());
    isFirstRecv.resize(gnss->receivers.size(), TRUE);
    for(auto recv : gnss->receivers)
    {
      const UInt idRecv = recv->idRecv();
      if(recv->isMyRank())
      {
        sigmaEpochRecv.at(idRecv).resize(gnss->times.size(), 1.);

        Vector l(gnss->times.size());
        Matrix A(gnss->times.size(), 2);
        for(UInt idEpoch=0; idEpoch<gnss->times.size(); idEpoch++)
          if(recv->useable(idEpoch))
          {
            l(idEpoch)    = LIGHT_VELOCITY * recv->clockError(idEpoch);
            A(idEpoch, 0) = 1;
            A(idEpoch, 1) = (gnss->times.at(idEpoch)-gnss->times.front()).seconds();
          }
        Vector x = leastSquares(Matrix(A), l);
        clock0Recv(idRecv) = x(0);
        driftRecv(idRecv)  = x(1);
      }
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssParametrizationClocksModel::requirements(GnssNormalEquationInfo &normalEquationInfo,
                                                  std::vector<UInt> &transCount, std::vector<UInt> &/*transCountEpoch*/,
                                                  std::vector<UInt> &recvCount,  std::vector<UInt> &/*recvCountEpoch*/)
{
  try
  {
    if(!isEnabled(normalEquationInfo, name))
      return;

    // transmitter epochs
    if(!normalEquationInfo.isEachReceiverSeparately)
      for(auto trans : gnss->transmitters)
        if(trans->useable() && selectedTransmitters.at(trans->idTrans()))
          transCount.at(trans->idTrans()) += 2; // to determine the drift

    // receiver clocks
    for(auto recv : gnss->receivers)
      if(recv->useable() && normalEquationInfo.estimateReceiver.at(recv->idRecv()) && selectedReceivers.at(recv->idRecv()))
        recvCount.at(recv->idRecv()) += 2; // to determine the drift
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssParametrizationClocksModel::initParameter(GnssNormalEquationInfo &normalEquationInfo)
{
  try
  {
    indexTrans.clear();      indexTrans.resize(gnss->transmitters.size());
    indexDriftTrans.clear(); indexDriftTrans.resize(gnss->transmitters.size());
    indexRecv.clear();       indexRecv.resize(gnss->receivers.size());
    indexDriftRecv.clear();  indexDriftRecv.resize(gnss->receivers.size());
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
          {
            indexTrans.at(idTrans).at(idEpoch) = normalEquationInfo.parameterNamesEpochTransmitter(idEpoch, idTrans, {ParameterName(trans->name(), "clock", "", gnss->times.at(idEpoch))});
            countParaTrans++;
          }
          indexDriftTrans.at(idTrans) = normalEquationInfo.parameterNamesTransmitter(idTrans, {ParameterName(trans->name(), "clockDrift")});
          countParaTrans++;
        }
      }
    if(countParaTrans)
      logInfo<<countParaTrans%"%9i modeled transmitter clock error and drift parameters"s<<Log::endl;

    UInt countParaRecv = 0;
    for(auto recv : gnss->receivers)
    {
      const UInt idRecv = recv->idRecv();
      if(recv->useable() && normalEquationInfo.estimateReceiver.at(idRecv) && selectedReceivers.at(idRecv))
      {
        indexRecv.at(idRecv).resize(gnss->times.size());
        for(UInt idEpoch : normalEquationInfo.idEpochs)
        {
          indexRecv.at(idRecv).at(idEpoch) = normalEquationInfo.parameterNamesEpochReceiver(idEpoch, idRecv, {ParameterName(recv->name(), "clock", "", gnss->times.at(idEpoch))});
          countParaRecv++;
        }
        indexDriftRecv.at(idRecv) = normalEquationInfo.parameterNamesReceiver(idRecv, {ParameterName(recv->name(), "clockDrift")});
        countParaRecv++;
      }
    }
    if(countParaRecv)
      logInfo<<countParaRecv%"%9i modeled receiver clock error and drift parameters"s<<Log::endl;

    applyConstraint = isEnabled(normalEquationInfo, nameConstraint) && sigmaZeroMean
                    && (countParaTrans+countParaRecv);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssParametrizationClocksModel::aprioriParameter(const GnssNormalEquationInfo &normalEquationInfo, MatrixSliceRef x0) const
{
  try
  {
    for(auto trans : gnss->transmitters)
      if(isMyRank(trans->idTrans()))
      {
        const UInt idTrans = trans->idTrans();
        for(UInt idEpoch : normalEquationInfo.idEpochs)
          if(indexTrans.at(idTrans).size() && indexTrans.at(idTrans).at(idEpoch))
            x0(normalEquationInfo.index(indexTrans.at(idTrans).at(idEpoch)), 0) = LIGHT_VELOCITY * trans->clockError(idEpoch);
        if(indexDriftTrans.at(idTrans))
          x0(normalEquationInfo.index(indexDriftTrans.at(idTrans)), 0) = driftTrans(idTrans);
      }

    for(auto recv : gnss->receivers)
      if(recv->isMyRank())
      {
        const UInt idRecv = recv->idRecv();
        for(UInt idEpoch : normalEquationInfo.idEpochs)
          if(indexRecv.at(idRecv).size() && indexRecv.at(idRecv).at(idEpoch))
            x0(normalEquationInfo.index(indexRecv.at(idRecv).at(idEpoch)), 0) = LIGHT_VELOCITY * recv->clockError(idEpoch);
        if(indexDriftRecv.at(idRecv))
          x0(normalEquationInfo.index(indexDriftRecv.at(idRecv)), 0) = driftRecv(idRecv);
      }
    }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssParametrizationClocksModel::designMatrix(const GnssNormalEquationInfo &/*normalEquationInfo*/, const GnssObservationEquation &eqn, GnssDesignMatrix &A) const
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

void GnssParametrizationClocksModel::constraintsEpoch(const GnssNormalEquationInfo &normalEquationInfo, UInt idEpoch, MatrixDistributed &normals, std::vector<Matrix> &n, Double &lPl, UInt &obsCount) const
{
  try
  {
    if(!applyConstraint)
      return;

    // zero mean of first and last epoch
    // ---------------------------------
    if(((idEpoch == normalEquationInfo.idEpochs.front()) || (idEpoch == normalEquationInfo.idEpochs.back())) && !normalEquationInfo.isEachReceiverSeparately)
    {
      UInt countTrans = 0, countRecv = 0;
      for(UInt idTrans=0; idTrans<indexTrans.size(); idTrans++)
        if(selectedTransmittersZeroMean.at(idTrans) && indexTrans.at(idTrans).size() && indexTrans.at(idTrans).at(idEpoch))
          countTrans++;
      for(UInt idRecv=0; idRecv<indexRecv.size(); idRecv++)
        if(selectedReceiversZeroMean.at(idRecv) && indexRecv.at(idRecv).size() && indexRecv.at(idRecv).at(idEpoch))
          countRecv++;

      if(countTrans + countRecv)
      {
        // collect apriori clock errors
        Vector lTransEpoch;
        if(countTrans)
        {
          lTransEpoch = Vector(gnss->transmitters.size());
          for(UInt idTrans=0; idTrans<indexTrans.size(); idTrans++)
            if(selectedTransmittersZeroMean.at(idTrans) && indexTrans.at(idTrans).size() && indexTrans.at(idTrans).at(idEpoch) && isMyRank(idTrans))
              lTransEpoch(idTrans) = clock0Trans(idTrans) + driftTrans(idTrans) * (gnss->times.at(idEpoch)-gnss->times.front()).seconds()
                                   - LIGHT_VELOCITY*gnss->transmitters.at(idTrans)->clockError(idEpoch);
          Parallel::reduceSum(lTransEpoch, 0, normalEquationInfo.comm);
        }

        // collect apriori clock errors
        Vector lRecvEpoch;
        if(countRecv)
        {
          lRecvEpoch = Vector(gnss->receivers.size());
          for(UInt idRecv=0; idRecv<indexRecv.size(); idRecv++)
            if(selectedReceiversZeroMean.at(idRecv) && indexRecv.at(idRecv).size() && indexRecv.at(idRecv).at(idEpoch) && gnss->receivers.at(idRecv)->isMyRank())
              lRecvEpoch(idRecv) = clock0Recv(idRecv) + driftRecv(idRecv) * (gnss->times.at(idEpoch)-gnss->times.front()).seconds()
                                 - LIGHT_VELOCITY*gnss->receivers.at(idRecv)->clockError(idEpoch);
          Parallel::reduceSum(lRecvEpoch, 0, normalEquationInfo.comm);
        }

        if(Parallel::isMaster(normalEquationInfo.comm))
        {
          const Double p = 1./(countTrans + countRecv)/sigmaZeroMean; // weight;
          GnssDesignMatrix A(normalEquationInfo, Vector(1));
          for(UInt idTrans=0; idTrans<indexTrans.size(); idTrans++)
            if(selectedTransmittersZeroMean.at(idTrans) && indexTrans.at(idTrans).size() && indexTrans.at(idTrans).at(idEpoch))
            {
              A.l(0) += lTransEpoch(idTrans) * p; // remove apriori value -> regularization towards 0
              A.column(indexTrans.at(idTrans).at(idEpoch))(0,0) = p;
            }
          for(UInt idRecv=0; idRecv<indexRecv.size(); idRecv++)
            if(selectedReceiversZeroMean.at(idRecv) && indexRecv.at(idRecv).size() && indexRecv.at(idRecv).at(idEpoch))
            {
              A.l(0) += lRecvEpoch(idRecv) * p; // remove apriori value -> regularization towards 0
              A.column(indexRecv.at(idRecv).at(idEpoch))(0,0) = p;
            }
          A.accumulateNormals(normals, n, lPl, obsCount);
        }
      } // if(countTrans + countRecv)
    } // if(first or last epoch)

    // regularize unused receiver clocks
    // ---------------------------------
//     for(UInt idRecv=0; idRecv<indexRecv.size(); idRecv++)
//       if(!gnss->receivers.at(idRecv)->useable(idEpoch))
//         if(indexRecv.at(idRecv).size() && indexRecv.at(idRecv).at(idEpoch) && gnss->receivers.at(idRecv)->isMyRank())
//         {
//           GnssDesignMatrix A(normalEquationInfo, Vector(1));
//           A.column(indexRecv.at(idRecv).at(idEpoch))(0,0) = 1./10.; // sigma = 10m
//           A.accumulateNormals(normals, n, lPl, obsCount);
//        }

    // find next epoch
    auto iter = std::upper_bound(normalEquationInfo.idEpochs.begin(), normalEquationInfo.idEpochs.end(), idEpoch);
    if(iter == normalEquationInfo.idEpochs.end())
      return;
    UInt idEpochNext = *iter;

    // transmitter clock model
    // -----------------------
    for(UInt idTrans=0; idTrans<indexTrans.size(); idTrans++)
      if(indexTrans.at(idTrans).size() && indexTrans.at(idTrans).at(idEpoch) && isMyRank(idTrans))
      {
        const Double p  = 1./sigma0Trans(idTrans)/sigmaEpochTrans.at(idTrans).at(idEpoch);
        const Double dt = (gnss->times.at(idEpochNext)-gnss->times.at(idEpoch)).seconds();
        GnssDesignMatrix A(normalEquationInfo, Vector(1));
        A.l(0) = (LIGHT_VELOCITY * gnss->transmitters.at(idTrans)->clockError(idEpochNext)
                - LIGHT_VELOCITY * gnss->transmitters.at(idTrans)->clockError(idEpoch)
                - driftTrans(idTrans)*dt) * p;
        A.column(indexTrans.at(idTrans).at(idEpochNext))(0,0) = -p;
        A.column(indexTrans.at(idTrans).at(idEpoch))(0,0)     = +p;
        A.column(indexDriftTrans.at(idTrans))(0,0)            = +p*dt;
        A.accumulateNormals(normals, n, lPl, obsCount);
      }

    // receiver clock model
    // --------------------
    for(UInt idRecv=0; idRecv<indexRecv.size(); idRecv++)
      if(indexRecv.at(idRecv).size() && indexRecv.at(idRecv).at(idEpoch) && gnss->receivers.at(idRecv)->isMyRank())
      {
        const Double p  = 1./sigma0Recv(idRecv)/sigmaEpochRecv.at(idRecv).at(idEpoch);
        const Double dt = (gnss->times.at(idEpochNext)-gnss->times.at(idEpoch)).seconds();
        GnssDesignMatrix A(normalEquationInfo, Vector(1));
        A.l(0) = (LIGHT_VELOCITY * gnss->receivers.at(idRecv)->clockError(idEpochNext)
                - LIGHT_VELOCITY * gnss->receivers.at(idRecv)->clockError(idEpoch)
                - driftRecv(idRecv)*dt) * p;
        A.column(indexRecv.at(idRecv).at(idEpochNext))(0,0) = -p;
        A.column(indexRecv.at(idRecv).at(idEpoch))(0,0)     = +p;
        A.column(indexDriftRecv.at(idRecv))(0,0)            = +p*dt;
        A.accumulateNormals(normals, n, lPl, obsCount);
      }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Double GnssParametrizationClocksModel::updateParameter(const GnssNormalEquationInfo &normalEquationInfo, const_MatrixSliceRef x, const_MatrixSliceRef Wz)
{
  try
  {
    Gnss::InfoParameterChange infoTrans("mm");
    for(auto trans : gnss->transmitters)
      if(indexTrans.at(trans->idTrans()).size())
        for(UInt idEpoch : normalEquationInfo.idEpochs)
        {
          const Double dClock = x(normalEquationInfo.index(indexTrans.at(trans->idTrans()).at(idEpoch)), 0);
          trans->updateClockError(idEpoch, dClock/LIGHT_VELOCITY);
          if(infoTrans.update(1e3*dClock))
            infoTrans.info = "clock transmitter ("+trans->name()+", "+gnss->times.at(idEpoch).dateTimeStr()+")";
        }

    Gnss::InfoParameterChange infoRecv("mm");
    for(auto recv : gnss->receivers)
      if(indexRecv.at(recv->idRecv()).size() && recv->isMyRank())
        for(UInt idEpoch : normalEquationInfo.idEpochs)
        {
          const Double dClock = x(normalEquationInfo.index(indexRecv.at(recv->idRecv()).at(idEpoch)), 0);
          recv->updateClockError(idEpoch, dClock/LIGHT_VELOCITY);
          if(recv->useable(idEpoch) && infoRecv.update(1e3*dClock))
            infoRecv.info = "clock receiver ("+recv->name()+", "+gnss->times.at(idEpoch).dateTimeStr()+")";
        }

    // determine variance of models
    // ----------------------------
    // transmitters
    Gnss::InfoParameterChange infoDriftTrans("cm/d");
    for(UInt idTrans=0; idTrans<indexTrans.size(); idTrans++)
      if(indexTrans.at(idTrans).size() && isMyRank(idTrans))
      {
        // update drift parameter
        const Double dDrift = x(normalEquationInfo.index(indexDriftTrans.at(idTrans)), 0);
        driftTrans(idTrans) += dDrift;
        if(infoDriftTrans.update(1e2*86400.*dDrift))
          infoDriftTrans.info = "clock drift transmitter ("+gnss->transmitters.at(idTrans)->name()+")";

        // compute model residuals
        Vector e(normalEquationInfo.idEpochs.size()-1);
        Vector r(normalEquationInfo.idEpochs.size()-1);
        Double ePe = 0;
        for(UInt i=0; i<normalEquationInfo.idEpochs.size()-1; i++)
        {
          const UInt idEpoch     = normalEquationInfo.idEpochs.at(i);
          const UInt idEpochNext = normalEquationInfo.idEpochs.at(i+1);
          const Double dt = (gnss->times.at(idEpochNext)-gnss->times.at(idEpoch)).seconds();
          e(i) = std::fabs(LIGHT_VELOCITY * gnss->transmitters.at(idTrans)->clockError(idEpochNext)
                         - LIGHT_VELOCITY * gnss->transmitters.at(idTrans)->clockError(idEpoch)
                        - driftTrans(idTrans)*dt);

          // redundancy
          const Double p  = 1./sigma0Trans(idTrans)/sigmaEpochTrans.at(idTrans).at(idEpoch);
          Matrix AWz(1, Wz.columns());
          axpy(+p,   Wz.row(normalEquationInfo.index(indexTrans.at(idTrans).at(idEpochNext))), AWz);
          axpy(-p,   Wz.row(normalEquationInfo.index(indexTrans.at(idTrans).at(idEpoch))),     AWz);
          axpy(p*dt, Wz.row(normalEquationInfo.index(indexDriftTrans.at(idTrans))),            AWz);
          ePe  += std::pow(e(i)*p, 2);
          r(i)  = std::max(1-quadsum(AWz), 0.001);
          e(i) /= std::sqrt(r(i));
        }

        // determine variances
        sigma0Trans(idTrans) *= Vce::standardDeviation(ePe, sum(r), huber, huberPower);
        if(isFirstTrans.at(idTrans))
          sigma0Trans(idTrans) = 1.4826 * median(e);
        for(UInt i=0; i<normalEquationInfo.idEpochs.size()-1; i++)
        {
          const UInt idEpoch = normalEquationInfo.idEpochs.at(i);
          sigmaEpochTrans.at(idTrans).at(idEpoch) = 1.;
          if((e(i) > sigma0Trans(idTrans)*huber) && (r(i) > 0.1))
            sigmaEpochTrans.at(idTrans).at(idEpoch) *= std::pow(e(i)/sigma0Trans(idTrans)/huber, huberPower);
        }
        isFirstTrans.at(idTrans) = FALSE;
      }

    // receivers
    Gnss::InfoParameterChange infoDriftRecv("cm/d");
    for(UInt idRecv=0; idRecv<indexRecv.size(); idRecv++)
      if(indexRecv.at(idRecv).size() && gnss->receivers.at(idRecv)->isMyRank())
      {
        // update drift parameter
        Double dDrift = x(normalEquationInfo.index(indexDriftRecv.at(idRecv)), 0);
        if(isFirstRecv.at(idRecv))
        {
          Vector dDrifts(normalEquationInfo.idEpochs.size()-1);
          for(UInt i=0; i<normalEquationInfo.idEpochs.size()-1; i++)
          {
            const UInt idEpoch     = normalEquationInfo.idEpochs.at(i);
            const UInt idEpochNext = normalEquationInfo.idEpochs.at(i+1);
            dDrifts(i) = LIGHT_VELOCITY/(gnss->times.at(idEpochNext)-gnss->times.at(idEpoch)).seconds() * (gnss->receivers.at(idRecv)->clockError(idEpochNext)-gnss->receivers.at(idRecv)->clockError(idEpoch));
          }
          dDrift = median(dDrifts);
        }
        driftRecv(idRecv) += dDrift;
        if((sigma0Recv(idRecv) < 1) && infoDriftRecv.update(1e2*86400.*dDrift))
          infoDriftRecv.info = "clock drift receiver ("+gnss->receivers.at(idRecv)->name()+")";

        // compute model residuals
        Vector e(normalEquationInfo.idEpochs.size()-1);
        Vector r(normalEquationInfo.idEpochs.size()-1);
        Double ePe = 0;
        for(UInt i=0; i<normalEquationInfo.idEpochs.size()-1; i++)
        {
          const UInt idEpoch     = normalEquationInfo.idEpochs.at(i);
          const UInt idEpochNext = normalEquationInfo.idEpochs.at(i+1);
          const Double dt = (gnss->times.at(idEpochNext)-gnss->times.at(idEpoch)).seconds();
          e(i) = std::fabs(LIGHT_VELOCITY * gnss->receivers.at(idRecv)->clockError(idEpochNext)
                         - LIGHT_VELOCITY * gnss->receivers.at(idRecv)->clockError(idEpoch)
                         - driftRecv(idRecv)*dt);

          // redundancy
          const Double p  = 1./sigma0Recv(idRecv)/sigmaEpochRecv.at(idRecv).at(idEpoch);
          Matrix AWz(1, Wz.columns());
          axpy(+p,   Wz.row(normalEquationInfo.index(indexRecv.at(idRecv).at(idEpochNext))), AWz);
          axpy(-p,   Wz.row(normalEquationInfo.index(indexRecv.at(idRecv).at(idEpoch))),     AWz);
          axpy(p*dt, Wz.row(normalEquationInfo.index(indexDriftRecv.at(idRecv))),            AWz);
          ePe  += std::pow(e(i)*p, 2);
          r(i)  = std::max(1-quadsum(AWz), 0.001);
          e(i) /= std::sqrt(r(i));
        }

        // determine variances
        sigma0Recv(idRecv) *= Vce::standardDeviation(ePe, sum(r), huber, huberPower);
        if(isFirstRecv.at(idRecv))
          sigma0Recv(idRecv) = 1.4826 * median(e);
        for(UInt i=0; i<normalEquationInfo.idEpochs.size()-1; i++)
        {
          const UInt idEpoch = normalEquationInfo.idEpochs.at(i);
          sigmaEpochRecv.at(idRecv).at(idEpoch) = 1.;
          if((e(i) > sigma0Recv(idRecv)*huber) && (r(i) > 0.1))
            sigmaEpochRecv.at(idRecv).at(idEpoch) *= std::pow(e(i)/sigma0Recv(idRecv)/huber, huberPower);
        }
        isFirstRecv.at(idRecv) = FALSE;
      }

    Double maxChange = 0;
    infoTrans.synchronizeAndPrint(normalEquationInfo.comm, 1e-3, maxChange);
    infoDriftTrans.synchronizeAndPrint(normalEquationInfo.comm, 1e-2/86400., maxChange);
    infoRecv.synchronizeAndPrint (normalEquationInfo.comm, 1e-3, maxChange);
    infoDriftRecv.synchronizeAndPrint (normalEquationInfo.comm, 1e-2/86400., maxChange);

    return maxChange;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssParametrizationClocksModel::writeResults(const GnssNormalEquationInfo &normalEquationInfo, const std::string &suffix) const
{
  try
  {
    if(!isEnabled(normalEquationInfo, name))
      return;

    if(!fileNameTransmitter.empty())
    {
      VariableList fileNameVariableList;
      addVariable("prn", "***", fileNameVariableList);
      logStatus<<"write transmitter clocks to files <"<<fileNameTransmitter(fileNameVariableList).appendBaseName(suffix)<<">"<<Log::endl;
      for(auto trans : gnss->transmitters)
        if(isMyRank(trans->idTrans()))
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
          fileNameVariableList["prn"]->setValue(trans->name());
          InstrumentFile::write(fileNameTransmitter(fileNameVariableList).appendBaseName(suffix), arc);
        }
    }

    if(!fileNameReceiver.empty())
    {
      VariableList fileNameVariableList;
      addVariable("station", "****", fileNameVariableList);
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
          fileNameVariableList["station"]->setValue(recv->name());
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
