/***********************************************/
/**
* @file gnssParametrizationCodeBiases.cpp
*
* @brief Code biases.
* @see GnssParametrization
*
* @author Torsten Mayer-Guerr
* @date 2021-01-23
*
*/
/***********************************************/

#include "base/import.h"
#include "config/config.h"
#include "classes/platformSelector/platformSelector.h"
#include "gnss/gnss.h"
#include "gnss/gnssParametrization/gnssParametrizationCodeBiases.h"

/***********************************************/

GnssParametrizationCodeBiases::GnssParametrizationCodeBiases(Config &config)
{
  try
  {
    readConfig(config, "name",                    name,               Config::OPTIONAL, "parameter.codeBiases", "used for parameter selection");
    readConfig(config, "selectTransmitters",      selectTransmitters, Config::DEFAULT,  R"(["all"])", "");
    readConfig(config, "selectReceivers",         selectReceivers,    Config::DEFAULT,  R"(["all"])", "");
    readConfig(config, "linearGlonassBias",       isLinearBias,       Config::DEFAULT,  "0", "bias depends linear on frequency channel number");
    readConfig(config, "typesClockDatum",         typesClockDatum,    Config::OPTIONAL, R"(["C1WG", "C2WG", "C1CE", "C5QE", "C1PR", "C2PR"])", "first two matching types define the ionosphere free transmitter clock (e.g. C1WG, C2WG)");
    readConfig(config, "nameConstraint",          nameConstraint,     Config::OPTIONAL, "constraint.codeBiases", "used for parameter selection");
    readConfig(config, "sigmaZeroMeanConstraint", sigmaZeroMean,      Config::DEFAULT,  "0.0001", "(0 = unconstrained) sigma [m] for null space constraint");
    if(isCreateSchema(config)) return;

    // mimic behaviour of old GROOPS version
    if(!typesClockDatum.size())
      typesClockDatum = {GnssType::C1WG, GnssType::C2WG};
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

GnssParametrizationCodeBiases::~GnssParametrizationCodeBiases()
{
  for(auto para : paraTrans)
    delete para;
  for(auto para : paraRecv)
    delete para;
}

/***********************************************/

void GnssParametrizationCodeBiases::init(Gnss *gnss, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    this->gnss = gnss;

    auto selectedTransmitters = gnss->selectTransmitters(selectTransmitters);
    paraTrans.resize(gnss->transmitters.size(), nullptr);
    for(UInt idTrans=0; idTrans<gnss->transmitters.size(); idTrans++)
      if(selectedTransmitters.at(idTrans) && gnss->transmitters.at(idTrans)->useable())
      {
        auto para = new Parameter();
        paraTrans.at(idTrans) = para;
        para->trans = gnss->transmitters.at(idTrans);
      }

    auto selectedReceivers = gnss->selectReceivers(selectReceivers);
    paraRecv.resize(gnss->receivers.size(), nullptr);
    for(UInt idRecv=0; idRecv<gnss->receivers.size(); idRecv++)
      if(selectedReceivers.at(idRecv) && gnss->receivers.at(idRecv)->useable())
      {
        auto para = new Parameter();
        paraRecv.at(idRecv) = para;
        para->recv = gnss->receivers.at(idRecv);
      }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssParametrizationCodeBiases::initParameter(GnssNormalEquationInfo &normalEquationInfo)
{
  try
  {
    for(auto para : paraTrans)
      if(para)
        para->index = GnssParameterIndex();
    for(auto para : paraRecv)
      if(para)
        para->index = GnssParameterIndex();
    applyConstraint = FALSE;
    if(!isEnabled(normalEquationInfo, name))
      return;


    auto nullSpace = [](Matrix &N, Bool remove)
    {
      Vector eigen = eigenValueDecomposition(N, TRUE);
      eigen *= (eigen(eigen.rows()-1) > 1e-4) ? 1./eigen(eigen.rows()-1) : 0.;
      UInt countZeros = 0;
      while((countZeros < eigen.rows()) && (eigen(countZeros) < 1e-8))
        countZeros++;
      return (remove) ? N.column(countZeros, N.columns()-countZeros) : N.column(0, countZeros);
    };

    // transmitter parameters
    // ----------------------
    UInt countParaTrans = 0;
    if(!normalEquationInfo.isEachReceiverSeparately)
      for(auto para : paraTrans)
        if(para && para->trans->useable())
        {
          std::vector<GnssType> biasTypes;
          for(GnssType type : para->trans->signalBias.types)
            if((type == GnssType::RANGE) && !type.isInList(biasTypes))
              biasTypes.push_back(type);
          if(biasTypes.size() <= 2)
            continue;

          // transformation matrix
          Matrix T(para->trans->signalBias.types.size(), biasTypes.size());
          for(UInt i=0; i<para->trans->signalBias.types.size(); i++)
            for(UInt k=0; k<biasTypes.size(); k++)
              if(para->trans->signalBias.types.at(i) == biasTypes.at(k))
                T(i, k) = 1.;

          // eliminate two biases which defines the clock
          UInt idx;
          GnssType firstFreq = GnssType::FREQUENCY;
          for(UInt i=0; i<typesClockDatum.size(); i++)
            if(typesClockDatum.at(i).isInList(biasTypes, idx) && (typesClockDatum.at(i) != firstFreq))
            {
              T(GnssType::index(para->trans->signalBias.types, typesClockDatum.at(i)), idx) = 0;
              if(firstFreq != GnssType::FREQUENCY)
                break; // found a second type
              firstFreq = typesClockDatum.at(i) & GnssType::FREQUENCY;
            }

          // special case: Galileo clock/TEC is defined via C1CE/C5QE only
          const UInt idxC1CE = GnssType::index(biasTypes, GnssType::C1CE);
          const UInt idxC5QE = GnssType::index(biasTypes, GnssType::C5QE);
          if((idxC1CE != NULLINDEX) && (idxC5QE != NULLINDEX))
          {
            T(GnssType::index(para->trans->signalBias.types, GnssType::C1CE), idxC1CE) = 0;
            T(GnssType::index(para->trans->signalBias.types, GnssType::C5QE), idxC5QE) = 0;
          }

          // determine nullspace
          Matrix N(1+T.columns(), Matrix::SYMMETRIC);
          for(UInt idRecv=0; idRecv<gnss->receivers.size(); idRecv++)
          {
            std::vector<GnssType> types;
            for(GnssType type : gnss->typesRecvTrans.at(idRecv).at(para->trans->idTrans()))
              if((type == GnssType::RANGE) && !type.isInList(types))
                types.push_back(type);
            if(!types.size())
              continue;

            // Composed signals (e.g. C2DG)
            std::vector<GnssType> typesTrans;
            Matrix Compose;
            gnss->receivers.at(idRecv)->signalComposition(NULLINDEX, types, typesTrans, Compose);

            Matrix A(types.size(), 1+T.columns());
            for(UInt i=0; i<types.size(); i++)
              if(types.at(i) == GnssType::RANGE)
                A(i, 0) = -1.;  // clock
            UInt idx;
            for(UInt i=0; i<typesTrans.size(); i++)
              if(typesTrans.at(i).isInList(para->trans->signalBias.types, idx))
                matMult(1., Compose.column(i), T.row(idx), A.column(1, T.columns()));
            Vector STEC(types.size());
            for(UInt i=0; i<types.size(); i++)
              STEC(i) = types.at(i).ionosphericFactor();
            eliminationParameter(STEC, {A});
            rankKUpdate(1., A, N);
          }

          // eliminate clock
          rankKUpdate(-1./N(0,0), N.slice(0, 1, 1, N.rows()-1), N.slice(1, 1, N.rows()-1, N.rows()-1));
          N = N.slice(1, 1, N.rows()-1, N.rows()-1); // without clock

          para->Bias = T * nullSpace(N, TRUE);
          if(!para->Bias.size())
            continue;

          // determine parameter names
          std::vector<ParameterName> parameterNames;
          for(UInt i=0; i<para->Bias.columns(); i++)
          {
            std::string typeStr;
            for(UInt idType=0; idType<para->Bias.rows(); idType++)
              if(std::fabs(para->Bias(idType, i)) > 1e-4)
                typeStr += ((para->Bias(idType, i) > 0) ? "+" : "") + para->Bias(idType, i)%"%.2f"s + para->trans->signalBias.types.at(idType).str();
            parameterNames.push_back(ParameterName(para->trans->name(), "codeBias"+(i+1)%"%02i("s+typeStr+")"));
          }
          para->index = normalEquationInfo.parameterNamesTransmitter(para->trans->idTrans(), parameterNames);
          countParaTrans += parameterNames.size();
        }
    if(countParaTrans)
      logInfo<<countParaTrans%"%9i transmitter code bias parameters"s<<Log::endl;

    // receiver parameters
    // -------------------
    UInt countParaRecv = 0;
    for(auto para : paraRecv)
      if(para && para->recv->useable() && normalEquationInfo.estimateReceiver.at(para->recv->idRecv()))
      {
        std::vector<GnssType> biasTypes;
        std::vector<GnssType> trendTypes;
        for(GnssType type : para->recv->signalBias.types)
          if(type == GnssType::RANGE)
          {
            if(!type.isInList(biasTypes))
              biasTypes.push_back(type & (isLinearBias ? ~GnssType::FREQ_NO : GnssType::ALL));
            else if(!type.isInList(trendTypes))
              trendTypes.push_back(type & ~GnssType::FREQ_NO);
          }
        if(biasTypes.size() <= 2)
          continue;

        // transformation matrix
        Matrix T(para->recv->signalBias.types.size(), biasTypes.size()+trendTypes.size());
        UInt idx;
        for(UInt i=0; i<para->recv->signalBias.types.size(); i++)
          if(para->recv->signalBias.types.at(i).isInList(biasTypes, idx))
            T(i, idx) = 1.;
        for(UInt i=0; i<para->recv->signalBias.types.size(); i++)
          if(para->recv->signalBias.types.at(i).isInList(trendTypes, idx))
            T(i, biasTypes.size()+idx) = para->recv->signalBias.types.at(i).frequencyNumber();

        // determine nullspace
        Matrix N(1+T.columns(), Matrix::SYMMETRIC);
        for(UInt idTrans=0; idTrans<gnss->transmitters.size(); idTrans++)
        {
          std::vector<GnssType> types;
          for(GnssType type : gnss->typesRecvTrans.at(para->recv->idRecv()).at(idTrans))
            if((type == GnssType::RANGE) && !type.isInList(types))
              types.push_back(type);
          if(!types.size())
            continue;

          Matrix A(types.size(), 1+T.columns());
          for(UInt i=0; i<types.size(); i++)
            if(types.at(i) == GnssType::RANGE)
              A(i, 0) = 1.;  // clock
          for(UInt i=0; i<types.size(); i++)
            copy(T.row(GnssType::index(para->recv->signalBias.types, types.at(i))), A.slice(i, 1, 1, T.columns()));
          Vector STEC(types.size());
          for(UInt i=0; i<types.size(); i++)
            STEC(i) = types.at(i).ionosphericFactor();
          eliminationParameter(STEC, {A});
          rankKUpdate(1., A, N);
        }

        // eliminate clock
        rankKUpdate(-1./N(0,0), N.slice(0, 1, 1, N.rows()-1), N.slice(1, 1, N.rows()-1, N.rows()-1));
        N = N.slice(1, 1, N.rows()-1, N.rows()-1); // without clock

        para->Bias = T * nullSpace(N, TRUE);
        if(!para->Bias.size())
          continue;

        // determine parameter names
        std::vector<ParameterName> parameterNames;
        for(UInt i=0; i<para->Bias.columns(); i++)
        {
          std::string typeStr;
          for(UInt idType=0; idType<para->Bias.rows(); idType++)
            if(std::fabs(para->Bias(idType, i)) > 1e-4)
              typeStr += ((para->Bias(idType, i) > 0) ? "+" : "") + para->Bias(idType, i)%"%.2f"s + para->recv->signalBias.types.at(idType).str();
          parameterNames.push_back(ParameterName(para->recv->name(), "codeBias"+(i+1)%"%02i("s+typeStr+")"));
        }
        para->index = normalEquationInfo.parameterNamesReceiver(para->recv->idRecv(), parameterNames);
        countParaRecv += parameterNames.size();
      }
    if(countParaRecv)
      logInfo<<countParaRecv%"%9i receiver code bias parameters"s<<Log::endl;

    applyConstraint = isEnabled(normalEquationInfo, nameConstraint) && sigmaZeroMean
                    && (countParaTrans+countParaRecv) &&  !normalEquationInfo.isEachReceiverSeparately;


    // calculate constraint equations (zero mean of signal biases)
    // -----------------------------------------------------------
    if(applyConstraint)
    {
      // parameter indices
      UInt parameterCount = 0;
      // indices of transmitter biases
      idxBiasTrans.clear();
      idxBiasTrans.resize(gnss->transmitters.size(), NULLINDEX);
      for(auto para : paraTrans)
        if(para && para->index)
        {
          idxBiasTrans.at(para->trans->idTrans()) = parameterCount;
          parameterCount += para->Bias.columns();
        }
      // indices of receiver biases
      idxBiasRecv.clear();
      idxBiasRecv.resize(gnss->receivers.size(), NULLINDEX);
      std::vector<std::vector<std::vector<GnssType>>> typesRecvTrans(gnss->receivers.size()); // for each receiver and transmitter: used types (receiver types)
      for(auto para : paraRecv)
        if(para && para->index)
        {
          const UInt idRecvOld = std::distance(typesRecvTrans.begin(), std::find(typesRecvTrans.begin(), typesRecvTrans.end(), gnss->typesRecvTrans.at(para->recv->idRecv())));
          if(idRecvOld >= typesRecvTrans.size())
          {
            typesRecvTrans.at(para->recv->idRecv()) = gnss->typesRecvTrans.at(para->recv->idRecv());
            idxBiasRecv.at(para->recv->idRecv()) = parameterCount;
            parameterCount += para->Bias.columns();
          }
          else
            idxBiasRecv.at(para->recv->idRecv()) = idxBiasRecv.at(idRecvOld);
        }

      if(parameterCount)
      {
        UInt idxClocks = parameterCount;
        // indicies of transmitter clocks
        std::vector<UInt> idxClockTrans(gnss->transmitters.size(), NULLINDEX);
        for(auto trans : gnss->transmitters)
          if(trans->useable())
            idxClockTrans.at(trans->idTrans()) = parameterCount++;
        // indices of receiver clocks
        std::vector<UInt> idxClockRecv(gnss->receivers.size(), NULLINDEX);
        typesRecvTrans.clear();
        typesRecvTrans.resize(gnss->receivers.size());
        for(auto recv : gnss->receivers)
          if(normalEquationInfo.estimateReceiver.at(recv->idRecv()))
          {
            const UInt idRecvOld = std::distance(typesRecvTrans.begin(), std::find(typesRecvTrans.begin(), typesRecvTrans.end(), gnss->typesRecvTrans.at(recv->idRecv())));
            if(idRecvOld >= typesRecvTrans.size())
            {
              typesRecvTrans.at(recv->idRecv()) = gnss->typesRecvTrans.at(recv->idRecv());
              idxClockRecv.at(recv->idRecv())   = parameterCount++;
            }
            else
              idxClockRecv.at(recv->idRecv()) = idxClockRecv.at(idRecvOld);
          }

        // normals of pseudo observations
        Matrix N(parameterCount, Matrix::SYMMETRIC);
        for(auto recv : gnss->receivers)
          if(recv->isMyRank() && normalEquationInfo.estimateReceiver.at(recv->idRecv()))
            for(auto trans : gnss->transmitters)
            {
              // observation types
              std::vector<GnssType> types;
              for(GnssType type : gnss->typesRecvTrans.at(recv->idRecv()).at(trans->idTrans()))
                if((type == GnssType::RANGE) && !type.isInList(types))
                  types.push_back(type);
              if(!types.size())
                continue;

              // transmitted types
              std::vector<GnssType> typesTrans;
              Matrix T;
              recv->signalComposition(NULLINDEX/*idEpoch*/, types, typesTrans, T);

              // design matrix
              Matrix A(types.size(), parameterCount);
              if(idxClockRecv.at(recv->idRecv()) != NULLINDEX)    // clock recv
                for(UInt i=0; i<types.size(); i++)
                  A(i, idxClockRecv.at(recv->idRecv())) = 1.;
              if(idxClockTrans.at(trans->idTrans()) != NULLINDEX) // clock trans
                for(UInt i=0; i<types.size(); i++)
                  A(i, idxClockTrans.at(trans->idTrans())) = -1.;
              if(idxBiasRecv.at(recv->idRecv()) != NULLINDEX)     // bias recv
                for(UInt i=0; i<types.size(); i++)
                  copy(paraRecv.at(recv->idRecv())->Bias.row(GnssType::index(recv->signalBias.types, types.at(i))),
                      A.slice(i, idxBiasRecv.at(recv->idRecv()), 1, paraRecv.at(recv->idRecv())->Bias.columns()));
              if(idxBiasTrans.at(trans->idTrans()) != NULLINDEX)  // bias trans
                for(UInt k=0; k<typesTrans.size(); k++)
                  matMult(1., T.column(k), paraTrans.at(trans->idTrans())->Bias.row(GnssType::index(trans->signalBias.types, typesTrans.at(k))),
                          A.column(idxBiasTrans.at(trans->idTrans()), paraTrans.at(trans->idTrans())->Bias.columns()));
              // eliminate STEC
              Vector STEC(types.size());
              for(UInt i=0; i<types.size(); i++)
                STEC(i) = types.at(i).ionosphericFactor();
              eliminationParameter(STEC, {A});
              rankKUpdate(1., A, N);
            }
        Parallel::reduceSum(N, 0, normalEquationInfo.comm);

        if(Parallel::isMaster(normalEquationInfo.comm))
        {
          // add zero mean of clocks and eliminate clocks
          Matrix N11 = N.slice(idxClocks, idxClocks, parameterCount-idxClocks, parameterCount-idxClocks);
          for(UInt i=0; i<N11.rows(); i++)
            if(N11(i,i) == 0)
              N11(i,i) = 1.;
          rankKUpdate(1., Matrix(1, N11.rows(), 1.), N11); // add zero mean
          cholesky(N11);
          triangularSolve(1., N11.trans(), N.slice(0, idxClocks, idxClocks, N11.rows()).trans());
          rankKUpdate(-1., N.slice(0, idxClocks, idxClocks, N11.rows()).trans(), N.slice(0, 0, idxClocks, idxClocks));
          N = N.slice(0, 0, idxClocks, idxClocks);

          zeroMeanDesign = nullSpace(N, FALSE).trans(); // null space defines the constraint equations
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

void GnssParametrizationCodeBiases::aprioriParameter(const GnssNormalEquationInfo &normalEquationInfo, MatrixSliceRef x0) const
{
  try
  {
    if(Parallel::isMaster(normalEquationInfo.comm))
      for(auto para : paraTrans)
        if(para && para->index)
          copy(leastSquares(Matrix(para->Bias), Vector(para->trans->signalBias.biases)), x0.row(normalEquationInfo.index(para->index), para->Bias.columns()));

    for(auto para : paraRecv)
      if(para && para->index && para->recv->isMyRank())
        copy(leastSquares(Matrix(para->Bias), Vector(para->recv->signalBias.biases)), x0.row(normalEquationInfo.index(para->index), para->Bias.columns()));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssParametrizationCodeBiases::designMatrix(const GnssNormalEquationInfo &/*normalEquationInfo*/, const GnssObservationEquation &eqn, GnssDesignMatrix &A) const
{
  try
  {
    auto paraTrans = this->paraTrans.at(eqn.transmitter->idTrans());
    if(paraTrans && paraTrans->index)
    {
      UInt idx;
      MatrixSlice Design(A.column(paraTrans->index));
      for(UInt idType=0; idType<eqn.typesTransmitted.size(); idType++)
        if((eqn.typesTransmitted.at(idType) == GnssType::RANGE) && eqn.typesTransmitted.at(idType).isInList(eqn.transmitter->signalBias.types, idx))
          matMult(1., eqn.A.column(GnssObservationEquation::idxUnit + eqn.types.size() + idType), paraTrans->Bias.row(idx), Design);
    }

    auto paraRecv = this->paraRecv.at(eqn.receiver->idRecv());
    if(paraRecv && paraRecv->index)
    {
      UInt idx;
      MatrixSlice Design(A.column(paraRecv->index));
      for(UInt idType=0; idType<eqn.types.size(); idType++)
        if((eqn.types.at(idType) == GnssType::RANGE) && eqn.types.at(idType).isInList(eqn.receiver->signalBias.types, idx))
          matMult(1., eqn.A.column(GnssObservationEquation::idxUnit + idType), paraRecv->Bias.row(idx), Design);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssParametrizationCodeBiases::constraints(const GnssNormalEquationInfo &normalEquationInfo, MatrixDistributed &normals, std::vector<Matrix> &n, Double &lPl, UInt &obsCount) const
{
  try
  {
    if(Parallel::isMaster(normalEquationInfo.comm) && applyConstraint)
    {
      logStatus<<"apply "<<zeroMeanDesign.rows()<<" zero mean equations for code bias parameters"<<Log::endl;
      if(zeroMeanDesign.size())
      {
        GnssDesignMatrix A(normalEquationInfo, Vector(zeroMeanDesign.rows()));
        for(auto para : paraTrans)
          if(para && para->index && (idxBiasTrans.at(para->trans->idTrans()) != NULLINDEX))
            axpy(1./sigmaZeroMean, zeroMeanDesign.column(idxBiasTrans.at(para->trans->idTrans()), para->Bias.columns()), A.column(para->index));
        for(auto para : paraRecv)
          if(para && para->index && (idxBiasRecv.at(para->recv->idRecv()) != NULLINDEX))
            axpy(1./sigmaZeroMean, zeroMeanDesign.column(idxBiasRecv.at(para->recv->idRecv()), para->Bias.columns()), A.column(para->index));
        A.accumulateNormals(normals, n, lPl, obsCount);
      }
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Double GnssParametrizationCodeBiases::updateParameter(const GnssNormalEquationInfo &normalEquationInfo, const_MatrixSliceRef x, const_MatrixSliceRef /*Wz*/)
{
  try
  {
    Double maxChange = 0;
    Gnss::InfoParameterChange infoTrans("mm");
    for(auto para : paraTrans)
      if(para && para->index)
      {
        const Vector dBias = para->Bias * x.row(normalEquationInfo.index(para->index), para->Bias.columns());
        for(UInt idType=0; idType<dBias.size(); idType++)
          para->trans->signalBias.biases.at(idType) += dBias(idType);
        for(UInt idType=0; idType<dBias.size(); idType++)
          if(infoTrans.update(1e3*dBias(idType)))
            infoTrans.info = "code bias transmitter ("+para->trans->signalBias.types.at(idType).str()+")";
      }
    infoTrans.synchronizeAndPrint(normalEquationInfo.comm, 1e-3, maxChange);

    Gnss::InfoParameterChange infoRecv("mm");
    for(auto para : paraRecv)
      if(para && para->index)
      {
        const Vector dBias = para->Bias * x.row(normalEquationInfo.index(para->index), para->Bias.columns());
        for(UInt idType=0; idType<dBias.size(); idType++)
          para->recv->signalBias.biases.at(idType) += dBias(idType);
        for(UInt idType=0; idType<dBias.size(); idType++)
          if(infoRecv.update(1e3*dBias(idType)))
            infoRecv.info = "code bias receiver ("+para->recv->name()+", "+para->recv->signalBias.types.at(idType).str()+")";
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
