/***********************************************/
/**
* @file gnssParametrizationTecBiases.cpp
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
#include "gnss/gnss.h"
#include "gnss/gnssTransceiverSelector/gnssTransceiverSelector.h"
#include "gnss/gnssParametrization/gnssParametrizationTecBiases.h"

/***********************************************/

GnssParametrizationTecBiases::GnssParametrizationTecBiases(Config &config)
{
  try
  {
    readConfig(config, "name",                    name,               Config::OPTIONAL, "parameter.tecBiases", "used for parameter selection");
    readConfig(config, "selectTransmitters",      selectTransmitters, Config::DEFAULT,  R"(["all"])", "");
    readConfig(config, "selectReceivers",         selectReceivers,    Config::DEFAULT,  R"(["all"])", "");
    readConfig(config, "linearGlonassBias",       isLinearBias,       Config::DEFAULT,  "0", "phase or code biases depend linear on frequency channel number");
    readConfig(config, "nameConstraint",          nameConstraint,     Config::OPTIONAL, "constraint.tecBiases", "used for parameter selection");
    readConfig(config, "sigmaZeroMeanConstraint", sigmaZeroMean,      Config::DEFAULT,  "0.0001", "(0 = unconstrained) sigma [m] for null space constraint");
    if(isCreateSchema(config)) return;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

GnssParametrizationTecBiases::~GnssParametrizationTecBiases()
{
  for(auto para : paraTrans)
    delete para;
  for(auto para : paraRecv)
    delete para;
}

/***********************************************/

void GnssParametrizationTecBiases::init(Gnss *gnss, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    this->gnss = gnss;

    auto selectedTransmitters = selectTransmitters->select(gnss->transmitters);
    paraTrans.resize(gnss->transmitters.size(), nullptr);
    for(UInt idTrans=0; idTrans<gnss->transmitters.size(); idTrans++)
      if(selectedTransmitters.at(idTrans) && gnss->transmitters.at(idTrans)->useable())
      {
        auto para = new Parameter();
        paraTrans.at(idTrans) = para;
        para->trans = gnss->transmitters.at(idTrans);
      }

    auto selectedReceivers = selectReceivers->select(gnss->receivers);
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

void GnssParametrizationTecBiases::initParameter(GnssNormalEquationInfo &normalEquationInfo)
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

    // transmitter parameters
    // ----------------------
    UInt countParaTrans = 0;
    if(!normalEquationInfo.isEachReceiverSeparately)
      for(auto para : paraTrans)
        if(para && para->trans->useable())
        {
          // determine nullspace
          Matrix N(para->trans->signalBias.types.size(), Matrix::SYMMETRIC);
          for(UInt idRecv=0; idRecv<gnss->receivers.size(); idRecv++)
          {
            std::vector<GnssType> types;
            for(GnssType type : gnss->typesRecvTrans.at(idRecv).at(para->trans->idTrans()))
              if(((type == GnssType::PHASE) || (type == GnssType::RANGE)) && !type.isInList(types))
                types.push_back(type);
            if(!types.size())
              continue;

            // Composed signals (e.g. C2DG)
            std::vector<GnssType> typesTrans;
            Matrix T;
            gnss->receivers.at(idRecv)->signalComposition(NULLINDEX, types, typesTrans, T);

            Matrix A(types.size(), para->trans->signalBias.types.size());
            UInt idx;
            for(UInt i=0; i<typesTrans.size(); i++)
              if(typesTrans.at(i).isInList(para->trans->signalBias.types, idx))
                axpy(1., T.column(i), A.column(idx));
            Vector STEC(types.size());
            for(UInt i=0; i<types.size(); i++)
              STEC(i) = types.at(i).ionosphericFactor();
            eliminationParameter(STEC, {A});
            rankKUpdate(1., A, N);
          }

          // determine eigen values
          Vector eigen = eigenValueDecomposition(N, TRUE);
          eigen *= (eigen(eigen.rows()-1) > 1e-4) ? 1./eigen(eigen.rows()-1) : 0.;
          UInt countZeros = 0;
          while((countZeros < eigen.rows()) && (eigen(countZeros) < 1e-8))
            countZeros++;
          para->Bias = N.column(0, countZeros);
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
            parameterNames.push_back(ParameterName(para->trans->name(), "tecBias"+(i+1)%"%02i("s+typeStr+")"));
          }
          para->index = normalEquationInfo.parameterNamesTransmitter(para->trans->idTrans(), parameterNames);
          countParaTrans += parameterNames.size();
        }
    if(countParaTrans)
      logInfo<<countParaTrans%"%9i transmitter TEC bias parameters"s<<Log::endl;

    // receiver parameters
    // -------------------
    UInt countParaRecv = 0;
    for(auto para : paraRecv)
      if(para && para->recv->useable() && normalEquationInfo.estimateReceiver.at(para->recv->idRecv()))
      {
        std::vector<GnssType> biasTypes;
        std::vector<GnssType> trendTypes;
        for(GnssType type : para->recv->signalBias.types)
          if((type == GnssType::RANGE) || (type == GnssType::PHASE))
          {
            if(!type.isInList(biasTypes))
              biasTypes.push_back(type & (isLinearBias ? ~GnssType::FREQ_NO : GnssType::ALL));
            else if(!type.isInList(trendTypes))
              trendTypes.push_back(type & ~GnssType::FREQ_NO);
          }

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
        Matrix N(T.columns(), Matrix::SYMMETRIC);
        for(UInt idTrans=0; idTrans<gnss->transmitters.size(); idTrans++)
        {
          std::vector<GnssType> types;
          for(GnssType type : gnss->typesRecvTrans.at(para->recv->idRecv()).at(idTrans))
            if(((type == GnssType::PHASE) || (type == GnssType::RANGE)) && !type.isInList(types))
              types.push_back(type);
          if(!types.size())
            continue;

          Matrix A(types.size(), T.columns());
          for(UInt i=0; i<types.size(); i++)
            copy(T.row(GnssType::index(para->recv->signalBias.types, types.at(i))), A.row(i));
          Vector STEC(types.size());
          for(UInt i=0; i<types.size(); i++)
            STEC(i) = types.at(i).ionosphericFactor();
          eliminationParameter(STEC, {A});
          rankKUpdate(1., A, N);
        }

        // determine eigen values
        Vector eigen = eigenValueDecomposition(N, TRUE);
        eigen *= (eigen(eigen.rows()-1) > 1e-4) ? 1./eigen(eigen.rows()-1) : 0.;
        UInt countZeros = 0;
        while((countZeros < eigen.rows()) && (eigen(countZeros) < 1e-8))
          countZeros++;
        para->Bias = T * N.column(0, countZeros);
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
          parameterNames.push_back(ParameterName(para->recv->name(), "tecBias"+(i+1)%"%02i("s+typeStr+")"));
        }
        para->index = normalEquationInfo.parameterNamesReceiver(para->recv->idRecv(), parameterNames);
        countParaRecv += parameterNames.size();
      }
    if(countParaRecv)
      logInfo<<countParaRecv%"%9i receiver TEC bias parameters"s<<Log::endl;

    applyConstraint = isEnabled(normalEquationInfo, nameConstraint) && sigmaZeroMean
                    && (countParaTrans+countParaRecv) &&  !normalEquationInfo.isEachReceiverSeparately;

    // calculate constraint equations (zero mean of signal biases)
    // -----------------------------------------------------------
    if(applyConstraint)
    {
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
              if(idxBiasRecv.at(recv->idRecv()) != NULLINDEX)     // bias recv
                for(UInt i=0; i<types.size(); i++)
                  copy(paraRecv.at(recv->idRecv())->Bias.row(GnssType::index(recv->signalBias.types, types.at(i))),
                      A.slice(i, idxBiasRecv.at(recv->idRecv()), 1, paraRecv.at(recv->idRecv())->Bias.columns()));
              if(idxBiasTrans.at(trans->idTrans()) != NULLINDEX)  // bias recv
                for(UInt k=0; k<typesTrans.size(); k++)
                  matMult(1., T.column(k), paraTrans.at(trans->idTrans())->Bias.row(GnssType::index(trans->signalBias.types, typesTrans.at(k))),
                          A.column(idxBiasTrans.at(trans->idTrans()), paraTrans.at(trans->idTrans())->Bias.columns()));
              rankKUpdate(1., A, N);
            }
        Parallel::reduceSum(N, 0, normalEquationInfo.comm);

        if(Parallel::isMaster(normalEquationInfo.comm))
        {
          // determine eigen values
          Vector eigen = eigenValueDecomposition(N, TRUE);
          eigen *= (eigen(eigen.rows()-1) > 1e-4) ? 1./eigen(eigen.rows()-1) : 0.;
          UInt countZeros = 0;
          while((countZeros < eigen.rows()) && (eigen(countZeros) < 1e-8))
            countZeros++;
          zeroMeanDesign = N.column(0, countZeros).trans(); // null space defines the constraint equations
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

void GnssParametrizationTecBiases::aprioriParameter(const GnssNormalEquationInfo &normalEquationInfo, MatrixSliceRef x0) const
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

void GnssParametrizationTecBiases::designMatrix(const GnssNormalEquationInfo &/*normalEquationInfo*/, const GnssObservationEquation &eqn, GnssDesignMatrix &A) const
{
  try
  {
    auto paraTrans = this->paraTrans.at(eqn.transmitter->idTrans());
    if(paraTrans && paraTrans->index)
    {
      MatrixSlice Design(A.column(paraTrans->index));
      UInt idx;
      for(UInt idType=0; idType<eqn.typesTransmitted.size(); idType++)
        if(eqn.typesTransmitted.at(idType).isInList(eqn.transmitter->signalBias.types, idx))
          matMult(1., eqn.A.column(GnssObservationEquation::idxUnit + eqn.types.size() + idType), paraTrans->Bias.row(idx), Design);
    }

    auto paraRecv = this->paraRecv.at(eqn.receiver->idRecv());
    if(paraRecv && paraRecv->index)
    {
      MatrixSlice Design(A.column(paraRecv->index));
      UInt idx;
      for(UInt idType=0; idType<eqn.types.size(); idType++)
        if(eqn.types.at(idType).isInList(eqn.receiver->signalBias.types, idx))
          matMult(1., eqn.A.column(GnssObservationEquation::idxUnit + idType), paraRecv->Bias.row(idx), Design);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssParametrizationTecBiases::constraints(const GnssNormalEquationInfo &normalEquationInfo, MatrixDistributed &normals, std::vector<Matrix> &n, Double &lPl, UInt &obsCount) const
{
  try
  {
    if(Parallel::isMaster(normalEquationInfo.comm) && applyConstraint)
    {
      logStatus<<"apply "<<zeroMeanDesign.rows()<<" zero mean equations for tec bias parameters"<<Log::endl;
      if(zeroMeanDesign.rows())
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

Double GnssParametrizationTecBiases::updateParameter(const GnssNormalEquationInfo &normalEquationInfo, const_MatrixSliceRef x, const_MatrixSliceRef /*Wz*/)
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
            infoTrans.info = "tec bias transmitter ("+para->trans->signalBias.types.at(idType).str()+")";
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
            infoRecv.info = "tec bias receiver ("+para->recv->name()+", "+para->recv->signalBias.types.at(idType).str()+")";
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
