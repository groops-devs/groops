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
#include "classes/platformSelector/platformSelector.h"
#include "misc/varianceComponentEstimation.h"
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
    readConfig(config, "huber",                      huber,                      Config::DEFAULT,  "4",   "clock jumps > huber*sigma0 are downweighted");
    readConfig(config, "huberPower",                 huberPower,                 Config::DEFAULT,  "1.5", "clock jumps > huber: sigma=(e/huber)^huberPower*sigma0");
    readConfig(config, "nameConstraint",             nameConstraint,             Config::OPTIONAL, "constraintEpoch.clocksModel", "used for parameter selection");
    readConfig(config, "selectTransmittersZeroMean", selectTransmittersZeroMean, Config::DEFAULT,  R"(["all"])", "use these transmitters for zero-mean constraint");
    readConfig(config, "selectReceiversZeroMean",    selectReceiversZeroMean,    Config::DEFAULT,  "", "use these receivers for zero-mean constraint");
    readConfig(config, "sigmaZeroMeanConstraint",    sigmaZeroMean,              Config::DEFAULT,  "10", "(0 = unconstrained) sigma [m] for zero-mean constraint over all selected clocks");
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
    selectedTransmitters         = gnss->selectTransmitters(selectTransmitters);
    selectedReceivers            = gnss->selectReceivers(selectReceivers);
    selectedTransmittersZeroMean = gnss->selectTransmitters(selectTransmittersZeroMean);
    selectedReceiversZeroMean    = gnss->selectReceivers(selectReceiversZeroMean);

    // determine apriori transmitter drift
    // -----------------------------------
    isMyRank     = Vector(gnss->transmitters.size());
    clock0Trans  = Vector(gnss->transmitters.size());
    drift0Trans  = Vector(gnss->transmitters.size());
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
        drift0Trans(idTrans) = x(1);
      }
    }
    driftTrans = drift0Trans;

    // determine apriori receiver drift
    // --------------------------------
    clock0Recv = Vector(gnss->receivers.size());
    drift0Recv  = Vector(gnss->receivers.size());
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
        drift0Recv(idRecv) = x(1);
      }
    }
    driftRecv = drift0Recv;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssParametrizationClocksModel::requirements(GnssNormalEquationInfo &normalEquationInfo,
                                                  std::vector<UInt> &transCount, std::vector<UInt> &transCountEpoch,
                                                  std::vector<UInt> &recvCount,  std::vector<UInt> &recvCountEpoch)
{
  try
  {
    if(!isEnabled(normalEquationInfo, name))
      return;

    // transmitter epochs
    if(!normalEquationInfo.isEachReceiverSeparately)
      for(auto trans : gnss->transmitters)
        if(trans->useable() && selectedTransmitters.at(trans->idTrans()))
        {
          transCountEpoch.at(trans->idTrans())++;
          transCount.at(trans->idTrans())++; // to determine the drift
        }

    // receiver clocks
    for(auto recv : gnss->receivers)
      if(recv->useable() && normalEquationInfo.estimateReceiver.at(recv->idRecv()) && selectedReceivers.at(recv->idRecv()))
      {
        recvCountEpoch.at(recv->idRecv())++;
        recvCount.at(recv->idRecv())++; // to determine the drift
      }
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
            if(trans->useable(idEpoch))
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
          if(recv->useable(idEpoch))
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
    // find next epoch
    auto iter = std::upper_bound(normalEquationInfo.idEpochs.begin(), normalEquationInfo.idEpochs.end(), idEpoch);
    if(iter != normalEquationInfo.idEpochs.end())
    {
      const UInt idEpochNext = *iter;

      // clock model
      // -----------
      auto clockModel = [&](const std::function<Double(UInt)> &clockError, Double drift,
                            const std::vector<GnssParameterIndex> &indexEpoch, const GnssParameterIndex &indexDrift,
                            Double sigma0, const std::vector<Double> &sigmaEpoch)
      {
        const Double p  = 1./sigma0/sigmaEpoch.at(idEpoch);
        const Double dt = (gnss->times.at(idEpochNext)-gnss->times.at(idEpoch)).seconds();
        Vector l(1);
        l(0) = p * (LIGHT_VELOCITY * (clockError(idEpochNext)-clockError(idEpoch)) - drift*dt);
        GnssDesignMatrix A(normalEquationInfo, 1);
        A.column(indexEpoch.at(idEpochNext))(0,0) = -p;
        A.column(indexEpoch.at(idEpoch))(0,0)     = +p;
        A.column(indexDrift)(0,0)                 = +p*dt;
        GnssDesignMatrix::accumulateNormals(A, l, normals, n, lPl, obsCount);
      };

      // transmitters
      for(UInt idTrans=0; idTrans<indexTrans.size(); idTrans++)
        if(isMyRank(idTrans) && indexTrans.at(idTrans).size() && indexTrans.at(idTrans).at(idEpoch) && indexTrans.at(idTrans).at(idEpochNext))
          clockModel(std::bind(&GnssTransmitter::clockError, gnss->transmitters.at(idTrans), std::placeholders::_1),
                     driftTrans(idTrans), indexTrans.at(idTrans), indexDriftTrans.at(idTrans),sigma0Trans(idTrans), sigmaEpochTrans.at(idTrans));

      // receivers
      for(UInt idRecv=0; idRecv<indexRecv.size(); idRecv++)
        if(gnss->receivers.at(idRecv)->isMyRank() && indexRecv.at(idRecv).size() && indexRecv.at(idRecv).at(idEpoch) && indexRecv.at(idRecv).at(idEpochNext))
          clockModel(std::bind(&GnssReceiver::clockError, gnss->receivers.at(idRecv), std::placeholders::_1),
                     driftRecv(idRecv), indexRecv.at(idRecv), indexDriftRecv.at(idRecv), sigma0Recv(idRecv), sigmaEpochRecv.at(idRecv));
    }

    // zero-mean constraint of clocks
    // ------------------------------
    if(!applyConstraint || normalEquationInfo.isEachReceiverSeparately)
      return;

    UInt count = 0;
    for(UInt idTrans=0; idTrans<indexTrans.size(); idTrans++)
      if(selectedTransmittersZeroMean.at(idTrans) && indexTrans.at(idTrans).size() && indexTrans.at(idTrans).at(idEpoch))
        count++;
    for(UInt idRecv=0; idRecv<indexRecv.size(); idRecv++)
      if(selectedReceiversZeroMean.at(idRecv) && indexRecv.at(idRecv).size() && indexRecv.at(idRecv).at(idEpoch))
        count++;
    if(!count)
      return;

    // collect apriori clock errors
    Double meanApriori = 0;
    for(UInt idTrans=0; idTrans<indexTrans.size(); idTrans++)
      if(selectedTransmittersZeroMean.at(idTrans) && indexTrans.at(idTrans).size() && indexTrans.at(idTrans).at(idEpoch) && isMyRank(idTrans))
        meanApriori += LIGHT_VELOCITY*gnss->transmitters.at(idTrans)->clockError(idEpoch)
                    -  (clock0Trans(idTrans) + drift0Trans(idTrans) * (gnss->times.at(idEpoch)-gnss->times.front()).seconds());
    for(UInt idRecv=0; idRecv<indexRecv.size(); idRecv++)
      if(selectedReceiversZeroMean.at(idRecv) && indexRecv.at(idRecv).size() && indexRecv.at(idRecv).at(idEpoch) && gnss->receivers.at(idRecv)->isMyRank())
        meanApriori += LIGHT_VELOCITY*gnss->receivers.at(idRecv)->clockError(idEpoch)
                    -  (clock0Recv(idRecv) + drift0Recv(idRecv) * (gnss->times.at(idEpoch)-gnss->times.front()).seconds());
    Parallel::reduceSum(meanApriori, 0, normalEquationInfo.comm);

    if(!Parallel::isMaster(normalEquationInfo.comm))
      return;

    const Double pSum = 1./count/sigmaZeroMean; // weight;
    Vector l(1);
    l(0) = -meanApriori * pSum;  // remove apriori value -> regularization towards 0
    GnssDesignMatrix A(normalEquationInfo, 1);
    for(UInt idTrans=0; idTrans<indexTrans.size(); idTrans++)
      if(selectedTransmittersZeroMean.at(idTrans) && indexTrans.at(idTrans).size() && indexTrans.at(idTrans).at(idEpoch))
        A.column(indexTrans.at(idTrans).at(idEpoch))(0,0) = pSum;
    for(UInt idRecv=0; idRecv<indexRecv.size(); idRecv++)
      if(selectedReceiversZeroMean.at(idRecv) && indexRecv.at(idRecv).size() && indexRecv.at(idRecv).at(idEpoch))
        A.column(indexRecv.at(idRecv).at(idEpoch))(0,0) = pSum;
    GnssDesignMatrix::accumulateNormals(A, l, normals, n, lPl, obsCount);
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
          if(indexTrans.at(trans->idTrans()).at(idEpoch))
          {
            const Double dClock = x(normalEquationInfo.index(indexTrans.at(trans->idTrans()).at(idEpoch)), 0);
            trans->updateClockError(idEpoch, dClock/LIGHT_VELOCITY);
            if(trans->useable(idEpoch) && infoTrans.update(1e3*dClock))
              infoTrans.info = "clock transmitter ("+trans->name()+", "+gnss->times.at(idEpoch).dateTimeStr()+")";
          }

    Gnss::InfoParameterChange infoRecv("mm");
    for(auto recv : gnss->receivers)
      if(indexRecv.at(recv->idRecv()).size() && recv->isMyRank())
        for(UInt idEpoch : normalEquationInfo.idEpochs)
          if(indexRecv.at(recv->idRecv()).at(idEpoch))
          {
            const Double dClock = x(normalEquationInfo.index(indexRecv.at(recv->idRecv()).at(idEpoch)), 0);
            recv->updateClockError(idEpoch, dClock/LIGHT_VELOCITY);
            if(recv->useable(idEpoch) && infoRecv.update(1e3*dClock))
              infoRecv.info = "clock receiver ("+recv->name()+", "+gnss->times.at(idEpoch).dateTimeStr()+")";
          }

    // transmitter drifts
    Gnss::InfoParameterChange infoDriftTrans("cm/d");
    for(UInt idTrans=0; idTrans<indexTrans.size(); idTrans++)
      if(indexDriftTrans.at(idTrans) && isMyRank(idTrans))
      {
        // update drift parameter
        const Double dDrift = x(normalEquationInfo.index(indexDriftTrans.at(idTrans)), 0);
        driftTrans(idTrans) += dDrift;
        if(infoDriftTrans.update(1e2*86400.*dDrift))
          infoDriftTrans.info = "clock drift transmitter ("+gnss->transmitters.at(idTrans)->name()+")";
      }

    // receivers drifts
    Gnss::InfoParameterChange infoDriftRecv("cm/d");
    for(UInt idRecv=0; idRecv<indexRecv.size(); idRecv++)
      if(indexDriftRecv.at(idRecv) && gnss->receivers.at(idRecv)->isMyRank())
      {
        const auto recv = gnss->receivers.at(idRecv);
        // update drift parameter
        Double dDrift = x(normalEquationInfo.index(indexDriftRecv.at(idRecv)), 0);
        if(isFirstRecv.at(idRecv))
        {
          UInt count = 0;
          Vector dDrifts(normalEquationInfo.idEpochs.size()-1);
          for(UInt i=0; i<normalEquationInfo.idEpochs.size()-1; i++)
          {
            const UInt idEpoch     = normalEquationInfo.idEpochs.at(i);
            const UInt idEpochNext = normalEquationInfo.idEpochs.at(i+1);
            if(indexRecv.at(recv->idRecv()).at(idEpoch) && indexRecv.at(recv->idRecv()).at(idEpochNext))
              dDrifts(count++) = LIGHT_VELOCITY/(gnss->times.at(idEpochNext)-gnss->times.at(idEpoch)).seconds()
                               * (recv->clockError(idEpochNext)-recv->clockError(idEpoch));
          }
          dDrift = median(dDrifts.row(0,count));
        }
        driftRecv(idRecv) += dDrift;
        if((sigma0Recv(idRecv) < 1) && infoDriftRecv.update(1e2*86400.*dDrift))
          infoDriftRecv.info = "clock drift receiver ("+recv->name()+")";
      }


    // =================================

    // determine variance of models
    // ----------------------------
    auto calcVariance = [&](const std::function<Double(UInt)> &clockError, Double drift,
                            const std::vector<GnssParameterIndex> &indexEpoch, const GnssParameterIndex &indexDrift,
                            Bool isFirstTime, Double &sigma0, std::vector<Double> &sigmaEpoch)
    {
      // compute model residuals
      Vector e(normalEquationInfo.idEpochs.size()-1);
      Vector r(normalEquationInfo.idEpochs.size()-1);
      Double ePe = 0;
      for(UInt i=0; i<normalEquationInfo.idEpochs.size()-1; i++)
      {
        const UInt idEpoch     = normalEquationInfo.idEpochs.at(i);
        const UInt idEpochNext = normalEquationInfo.idEpochs.at(i+1);
        if(!indexEpoch.at(idEpoch) || !indexEpoch.at(idEpochNext))
          continue;

        // residuals
        const Double dt = (gnss->times.at(idEpochNext)-gnss->times.at(idEpoch)).seconds();
        e(i) = std::fabs(LIGHT_VELOCITY * (clockError(idEpochNext) - clockError(idEpoch)) - drift*dt);

        // redundancy
        const Double p  = 1./sigma0/sigmaEpoch.at(idEpoch);
        Matrix AWz(1, Wz.columns());
        axpy(-p,    Wz.row(normalEquationInfo.index(indexEpoch.at(idEpochNext))), AWz);
        axpy(+p,    Wz.row(normalEquationInfo.index(indexEpoch.at(idEpoch))),     AWz);
        axpy(+p*dt, Wz.row(normalEquationInfo.index(indexDrift)),                 AWz);
        ePe  += std::pow(e(i)*p, 2);
        r(i)  = std::max(1-quadsum(AWz), 0.2);
        e(i) /= std::sqrt(r(i));
      }

      // determine variances
      sigma0 *= Vce::standardDeviation(ePe, sum(r), huber, huberPower);
      if(isFirstTime)
        sigma0 = 1.4826 * median(e);

      for(UInt i=0; i<normalEquationInfo.idEpochs.size()-1; i++)
      {
        const UInt idEpoch = normalEquationInfo.idEpochs.at(i);
        sigmaEpoch.at(idEpoch) = 1.;
        if(e(i) > sigma0*huber)
          sigmaEpoch.at(idEpoch) *= std::pow(e(i)/sigma0/huber, huberPower);
      }
    };

    // transmitters
    for(UInt idTrans=0; idTrans<indexTrans.size(); idTrans++)
      if(indexTrans.at(idTrans).size() && isMyRank(idTrans))
      {
        calcVariance(std::bind(&GnssTransmitter::clockError, gnss->transmitters.at(idTrans), std::placeholders::_1),
                     driftTrans(idTrans), indexTrans.at(idTrans), indexDriftTrans.at(idTrans),
                     isFirstTrans.at(idTrans), sigma0Trans(idTrans), sigmaEpochTrans.at(idTrans));
        isFirstTrans.at(idTrans) = FALSE;
      }

    // receivers
    for(UInt idRecv=0; idRecv<indexRecv.size(); idRecv++)
      if(indexRecv.at(idRecv).size() && gnss->receivers.at(idRecv)->isMyRank())
      {
        calcVariance(std::bind(&GnssReceiver::clockError, gnss->receivers.at(idRecv), std::placeholders::_1),
                     driftRecv(idRecv), indexRecv.at(idRecv), indexDriftRecv.at(idRecv),
                     isFirstRecv.at(idRecv), sigma0Recv(idRecv), sigmaEpochRecv.at(idRecv));
        isFirstRecv.at(idRecv) = FALSE;
      }

    // =================================

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
      fileNameVariableList.setVariable("prn", "***");
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
