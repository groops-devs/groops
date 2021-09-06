/***********************************************/
/**
* @file gnssParametrizationAmbiguities.cpp
*
* @brief integer and float ambiguities.
*
* @author Torsten Mayer-Guerr
* @author Sebastian Strasser
* @date 2013-06-24
*
*/
/***********************************************/

#include "base/import.h"
#include "config/config.h"
#include "inputOutput/logging.h"
#include "files/fileMatrix.h"
#include "gnss/gnssLambda.h"
#include "gnss/gnssReceiver.h"
#include "gnss/gnssParametrization/gnssParametrizationAmbiguities.h"


/***********************************************/

GnssParametrizationAmbiguities::GnssParametrizationAmbiguities(Config &config)
{
  try
  {
    std::string choice;

    readConfig(config, "name",                           name,                    Config::OPTIONAL, "parameter.ambiguities", "used for parameter selection");
    readConfig(config, "estimateTransmitterPhaseBias",   selectTransmitters,      Config::DEFAULT,  "[\"all\"]", "");
    readConfig(config, "estimateReceiverPhaseBias",      selectReceivers,         Config::DEFAULT,  "[\"all\"]", "");
    readConfig(config, "linearGlonassBias",              isLinearBias,            Config::DEFAULT,  "0", "bias depends linear on frequency channel number");
    if(isCreateSchema(config)) return;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

GnssParametrizationAmbiguities::~GnssParametrizationAmbiguities()
{
  for(auto para : paraTrans)
    delete para;
  for(auto para : paraRecv)
    delete para;
}

/***********************************************/

void GnssParametrizationAmbiguities::init(Gnss *gnss, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    this->gnss = gnss;

    // float bias at transmitter
    auto selectedTransmitters = selectTransmitters->select(gnss->transmitters);
    paraTrans.resize(gnss->transmitters.size(), nullptr);
    for(UInt idTrans=0; idTrans<gnss->transmitters.size(); idTrans++)
      if(selectedTransmitters.at(idTrans) && gnss->transmitters.at(idTrans)->useable())
      {
        auto para = new ParameterTrans();
        paraTrans.at(idTrans) = para;
        para->trans = gnss->transmitters.at(idTrans);
      }

    // float bias at receiver
    auto selectedReceivers = selectReceivers->select(gnss->receivers);
    paraRecv.resize(gnss->receivers.size(), nullptr);
    for(UInt idRecv=0; idRecv<gnss->receivers.size(); idRecv++)
      if(selectedReceivers.at(idRecv) && gnss->receivers.at(idRecv)->useable())
      {
        auto para = new ParameterRecv();
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

std::vector<GnssParametrizationAmbiguities::Ambiguity*> GnssParametrizationAmbiguities::getAmbiguities() const
{
  try
  {
    std::vector<Ambiguity*> ambiguities;
    for(auto recv : gnss->receivers)
      if(recv->isMyRank())
        for(auto &track : recv->tracks)
          if(track->ambiguity)
            ambiguities.push_back(dynamic_cast<Ambiguity*>(track->ambiguity));
    return ambiguities;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/

// to share information between processes
class GnssParametrizationAmbiguities::AmbiguityInfo
{
public:
  Ambiguity            *ambi;
  UInt                  idRecv, idTrans;
  UInt                  idEpochStart, idEpochEnd;
  UInt                  parameterCount;
  std::vector<GnssType> typesAmbiguity, typesTrack;
  std::vector<GnssType> typesNew;
  Bool                  processed;

  AmbiguityInfo() : ambi(nullptr) {}

  AmbiguityInfo(Ambiguity *ambi) :
    ambi(ambi), idRecv(ambi->track->receiver->idRecv()), idTrans(ambi->track->transmitter->idTrans()),
    idEpochStart(ambi->track->idEpochStart), idEpochEnd(ambi->track->idEpochEnd),
    parameterCount(ambi->value.rows()), typesAmbiguity(ambi->types), typesTrack(ambi->track->types) {}

  void save(OutArchive &oa) const
  {
    oa<<nameValue("idRecv",         idRecv);
    oa<<nameValue("idTrans",        idTrans);
    oa<<nameValue("idEpochStart",   idEpochStart);
    oa<<nameValue("idEpochEnd",     idEpochEnd);
    oa<<nameValue("parameterCount", parameterCount);
    oa<<nameValue("typesAmbiguity", typesAmbiguity);
    oa<<nameValue("typesTrack",     typesTrack);
  }

  void load(InArchive  &ia)
  {
    ia>>nameValue("idRecv",         idRecv);
    ia>>nameValue("idTrans",        idTrans);
    ia>>nameValue("idEpochStart",   idEpochStart);
    ia>>nameValue("idEpochEnd",     idEpochEnd);
    ia>>nameValue("parameterCount", parameterCount);
    ia>>nameValue("typesAmbiguity", typesAmbiguity);
    ia>>nameValue("typesTrack",     typesTrack);
  }
};

/***********************************************/
/***********************************************/

void GnssParametrizationAmbiguities::initParameter(GnssNormalEquationInfo &normalEquationInfo)
{
  try
  {
    // setup new ambiguities
    // ---------------------
    for(auto recv : gnss->receivers)
      if(recv->isMyRank())
      {
        recv->deleteEmptyTracks();
        for(auto &track : recv->tracks)
          if(!track->ambiguity)
            new Ambiguity(track.get()); // track is owner of ambiguity
      }

    // reset parameter indices
    // -----------------------
    for(auto para : paraTrans)
      if(para)
      {
        para->index = GnssParameterIndex();
        para->types.clear();
      }
    for(auto para : paraRecv)
      if(para)
      {
        para->index = GnssParameterIndex();
        para->types.clear();
      }
    for(auto ambi : getAmbiguities())
      ambi->index = GnssParameterIndex();

    if(!isEnabled(normalEquationInfo, name))
      return;

    // distribute ambiguities to all nodes
    // -----------------------------------
    std::vector<AmbiguityInfo> ambiguityInfos;
    for(UInt idProcess=0; idProcess<Parallel::size(normalEquationInfo.comm); idProcess++)
    {
      std::vector<AmbiguityInfo> ambiguityInfosLocal;
      if(Parallel::myRank(normalEquationInfo.comm) == idProcess)
        for(auto ambi : getAmbiguities())
          if(ambi->track->receiver->useable() && normalEquationInfo.estimateReceiver.at(ambi->track->receiver->idRecv()))
            ambiguityInfosLocal.emplace_back(ambi);
      // synchronize information between processes
      Parallel::broadCast(ambiguityInfosLocal, idProcess, normalEquationInfo.comm);
      ambiguityInfos.insert(ambiguityInfos.end(), ambiguityInfosLocal.begin(), ambiguityInfosLocal.end());
    }

    // float biases at receivers
    // -------------------------
    for(auto &info : ambiguityInfos)
      for(GnssType &type : info.typesTrack)
        if(paraRecv.at(info.idRecv) && !type.isInList(paraRecv.at(info.idRecv)->types))
          paraRecv.at(info.idRecv)->types.push_back(type & ~GnssType::PRN);

    std::vector<std::vector<GnssType>> typesRecv(paraRecv.size()); // float ambiguity at receiver
    UInt countParaRecv = 0;
    for(auto para : paraRecv)
      if(para)
      {
        auto &types = para->types;
        std::sort(types.begin(), types.end());
        std::vector<GnssType> typesBias;
        std::vector<GnssType> typesTrend;
        for(UInt idType=0; idType<types.size(); idType++)
        {
          // count same types with different frequency numbers
          UInt count = 1;
          while((idType+count < types.size()) && ((types.at(idType+count) & ~GnssType::FREQ_NO) == types.at(idType)))
            count++;
          // can type be represented by a trend? (need two adjacent freq. numbers)
          Bool isTrend = FALSE;
          if(isLinearBias && (count > 2))
            for(UInt i=0; i<count-1; i++)
              if((types.at(idType+i+1).frequencyNumber() - types.at(idType+i).frequencyNumber()) == 1)
              {
                isTrend = TRUE;
                typesRecv.at(para->recv->idRecv()).push_back(types.at(idType+i));
                typesRecv.at(para->recv->idRecv()).push_back(types.at(idType+i+1));
                break;
              }
          if(isTrend)
          {
            typesBias.push_back (types.at(idType) & ~GnssType::FREQ_NO);
            if(!types.at(idType).isInList(typesTrend))
              typesTrend.push_back(types.at(idType) & ~(GnssType::ATTRIBUTE+GnssType::FREQ_NO));
          }
          else
            for(UInt i=0; i<count; i++)
            {
              typesBias.push_back(types.at(idType+i));
              typesRecv.at(para->recv->idRecv()).push_back(types.at(idType+i));
            }
          idType += count-1;
        } // for(idType)

        // transformation matrix
        para->Bias = Matrix(types.size(), typesBias.size()+typesTrend.size());
        UInt idx;
        for(UInt i=0; i<types.size(); i++)
          para->Bias(i, GnssType::index(typesBias, types.at(i))) = LIGHT_VELOCITY/types.at(i).frequency();
        std::sort(typesTrend.begin(), typesTrend.end());
        for(UInt i=0; i<types.size(); i++)
          if(types.at(i).isInList(typesTrend, idx))
            para->Bias(i, typesBias.size()+idx) = LIGHT_VELOCITY/types.at(i).frequency() * types.at(i).frequencyNumber();

        // parameter names
        std::vector<ParameterName> parameterNames(typesBias.size()+typesTrend.size());
        for(UInt i=0; i<typesBias.size(); i++)
          parameterNames.at(i) = ParameterName(para->recv->name(), "phaseBias("+typesBias.at(i).str()+")");
        for(UInt i=0; i<typesTrend.size(); i++)
          parameterNames.at(typesBias.size()+i) = ParameterName(para->recv->name(), "phaseBiasTrend("+typesTrend.at(i).str()+")");
        para->index = normalEquationInfo.parameterNamesReceiver(para->recv->idRecv(), parameterNames);
        countParaRecv += parameterNames.size();
      } // for(paraRecv)

    // determine integer ambiguities (skip zero/single diff float ambiguities)
    // -----------------------------------------------------------------------
    // sort with decreasing track length
    std::stable_sort(ambiguityInfos.begin(), ambiguityInfos.end(), [](const AmbiguityInfo &a, const AmbiguityInfo &b)
                    {return (a.idEpochEnd-a.idEpochStart) > (b.idEpochEnd-b.idEpochStart);});
    // remove attributes
    for(auto &info : ambiguityInfos)
      info.typesTrack = GnssType::replaceCompositeSignals(info.typesTrack);

    std::vector<GnssType> types;
    for(auto &info : ambiguityInfos)
      for(GnssType &type : info.typesTrack)
        if(!type.isInList(types))
          types.push_back(type & ~GnssType::PRN);

    std::vector<GnssType> typesZero; // float ambiguity at receiver and transmitter
    for(GnssType type : types)
    {
      std::vector<Bool> floatTrans(paraTrans.size(), FALSE);
      std::vector<Bool> floatRecv(paraRecv.size(), FALSE);
      Bool              floatBoth = FALSE;
      for(auto &info : ambiguityInfos)
        info.processed = !type.isInList(info.typesTrack);
      for(;;)
      {
        Bool restart = FALSE;
        UInt idx;
        for(auto &info : ambiguityInfos)
          if(!info.processed && type.isInList(info.typesTrack, idx))
          {
            const Bool isFloatRecv  = !floatRecv.at(info.idRecv)   && type.isInList(typesRecv.at(info.idRecv));
            const Bool isFloatTrans = !floatTrans.at(info.idTrans) && paraTrans.at(info.idTrans) && !normalEquationInfo.isEachReceiverSeparately;
            if(isFloatTrans && isFloatRecv && floatBoth)                                        // only one zero difference ambiguity is allowed for each type
              continue;                                                                         // try again later when new single float ambiguities are set up
            if(!isFloatTrans && !isFloatRecv) info.typesNew.push_back(info.typesTrack.at(idx)); // integer double difference
            if(isFloatRecv)                   floatRecv.at(info.idRecv)   = TRUE;               // float receiver single difference
            if(isFloatTrans)                  floatTrans.at(info.idTrans) = TRUE;               // float transmitter single difference
            if(isFloatTrans && isFloatRecv)   floatBoth = TRUE;                                 // zero difference
            if(isFloatTrans && isFloatRecv)   typesZero.push_back(info.typesTrack.at(idx));     // zero difference
            info.processed = TRUE;
            restart = (isFloatTrans || isFloatRecv);                                            // new float ambiguity -> restart searching
            if(restart)
              break;
          }
        if(!restart)
          break;
      }

      for(auto &info : ambiguityInfos)
        if(!info.processed)
          throw(Exception("Cannot setup ambiguities. Network is separated into independent parts!?"));
    } // for(type)

    // ambiguity parameters names
    // --------------------------
    // sort in temporal order
    std::stable_sort(ambiguityInfos.begin(), ambiguityInfos.end(), [](const AmbiguityInfo &a, const AmbiguityInfo &b) {return (a.idEpochStart+a.idEpochEnd) < (b.idEpochStart+b.idEpochEnd);});
    UInt countParaAmbi = 0;
    for(auto &info : ambiguityInfos)
    {
      std::sort(info.typesNew.begin(), info.typesNew.end());
      if(info.typesNew != info.typesAmbiguity)
        info.parameterCount = info.typesNew.size();

      // parameter names
      std::vector<ParameterName> parameterNames;
      ParameterName name(gnss->receivers.at(info.idRecv)->name()+"."+gnss->transmitters.at(info.idTrans)->name(),
                         "ambiguity", "", gnss->times.at(info.idEpochStart), gnss->times.at(info.idEpochEnd));
      std::string typeStr;
      for(auto type : info.typesNew)
        typeStr += type.str().substr(0,2);
      for(UInt i=0; i<info.parameterCount; i++)
        parameterNames.push_back(ParameterName(name.object, name.type+(i+1)%"%iof"s+info.parameterCount%"%i("s+typeStr+")", name.temporal, name.interval));
      GnssParameterIndex index = normalEquationInfo.parameterNamesAmbiguity((info.idEpochStart+info.idEpochEnd)/2, info.idRecv, info.idTrans, parameterNames);
      countParaAmbi += parameterNames.size();

      if(info.ambi)
      {
        info.ambi->index     = index;
        info.ambi->isInteger = gnss->receivers.at(info.idRecv)->integerAmbiguities;
        if(info.ambi->types != info.typesNew)
        {
          info.ambi->types = info.typesNew;
          info.ambi->T     = GnssLambda::phaseDecorrelation(info.ambi->types, gnss->receivers.at(info.idRecv)->wavelengthFactor);
          info.ambi->value = Vector(info.ambi->T.columns());
        }
      }
    }

    // float biases at transmitters
    // ----------------------------
    if(!normalEquationInfo.isEachReceiverSeparately)
      for(auto &info : ambiguityInfos)
        for(GnssType &type : info.typesTrack)
          if(paraTrans.at(info.idTrans) && !type.isInList(paraTrans.at(info.idTrans)->types) && !type.isInList(typesZero))
            paraTrans.at(info.idTrans)->types.push_back(type);

    UInt countParaTrans = 0;
    for(auto para : paraTrans)
      if(para && para->types.size())
      {
        std::sort(para->types.begin(), para->types.end());
        std::vector<ParameterName> parameterNames(para->types.size());
        for(UInt i=0; i<para->types.size(); i++)
          parameterNames.at(i) = ParameterName(para->trans->name(), "phaseBias("+para->types.at(i).str()+")");
        para->index = normalEquationInfo.parameterNamesTransmitter(para->trans->idTrans(), parameterNames);
        countParaTrans += parameterNames.size();
      }

    if(countParaTrans) logInfo<<countParaTrans%"%9i transmitter phase bias parameters"s<<Log::endl;
    if(countParaRecv)  logInfo<<countParaRecv% "%9i receiver phase bias parameters"s<<Log::endl;
    if(countParaAmbi)  logInfo<<countParaAmbi% "%9i ambiguity parameters"s<<Log::endl;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssParametrizationAmbiguities::aprioriParameter(const GnssNormalEquationInfo &normalEquationInfo, MatrixSliceRef x0) const
{
  try
  {
    // float bias at transmitter
    if(Parallel::isMaster(normalEquationInfo.comm))
      for(auto para : paraTrans)
        if(para && para->index)
          copy(para->trans->signalBias.compute(para->types), x0.row(normalEquationInfo.index(para->index), para->types.size()));

    // float bias at receiver
    for(auto para : paraRecv)
      if(para && para->index && para->recv->isMyRank())
        copy(leastSquares(Matrix(para->Bias), para->recv->signalBias.compute(para->types)), x0.row(normalEquationInfo.index(para->index), para->Bias.columns()));

    // ambiguities
    for(auto ambi : getAmbiguities())
      if(ambi->index)
        copy(ambi->value, x0.row(normalEquationInfo.index(ambi->index), ambi->value.rows()));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssParametrizationAmbiguities::designMatrix(const GnssNormalEquationInfo &/*normalEquationInfo*/, const GnssObservationEquation &eqn, GnssDesignMatrix &A) const
{
  try
  {
    // float bias at transmitter
    auto paraTrans = this->paraTrans.at(eqn.transmitter->idTrans());
    if(paraTrans && paraTrans->index)
    {
      UInt idx;
      MatrixSlice Design(A.column(paraTrans->index));
      for(UInt idType=0; idType<eqn.typesTransmitted.size(); idType++)
        if(eqn.typesTransmitted.at(idType).isInList(paraTrans->types, idx))
          copy(eqn.A.column(GnssObservationEquation::idxUnit + eqn.types.size() + idType), Design.column(idx));
    }

    // float bias at receiver
    auto paraRecv = this->paraRecv.at(eqn.receiver->idRecv());
    if(paraRecv && paraRecv->index)
    {
      UInt idx;
      MatrixSlice Design(A.column(paraRecv->index));
      for(UInt idType=0; idType<eqn.types.size(); idType++)
        if(eqn.types.at(idType).isInList(paraRecv->types, idx))
          matMult(1., eqn.A.column(GnssObservationEquation::idxUnit + idType), paraRecv->Bias.row(idx), Design);
    }

    // ambiguities
    if(eqn.track && eqn.track->ambiguity)
    {
      const Ambiguity *ambiguity = dynamic_cast<Ambiguity*>(eqn.track->ambiguity);
      if(ambiguity->index)
      {
        UInt idx;
        MatrixSlice Design(A.column(ambiguity->index));
        for(UInt idType=0; idType<eqn.typesTransmitted.size(); idType++)
          if(eqn.typesTransmitted.at(idType).isInList(ambiguity->types, idx))
            matMult(1., eqn.A.column(GnssObservationEquation::idxUnit + eqn.types.size() + idType), ambiguity->T.row(idx), Design);
      }
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Double GnssParametrizationAmbiguities::ambiguityResolve(const GnssNormalEquationInfo &normalEquationInfo, MatrixDistributed &normals, std::vector<Matrix> &n, Double &lPl, UInt &obsCount,
                                                        const std::function<Vector(const_MatrixSliceRef xFloat, MatrixSliceRef W, const_MatrixSliceRef d, Vector &isNotFixed, Double &sigma)> &searchInteger)
{
  try
  {
    const UInt block0 = normalEquationInfo.blockAmbiguity();
    const UInt index0 = normalEquationInfo.blockIndex(block0);
    const UInt dim    = normalEquationInfo.parameterCount()-index0;
    Double sigmaFloat = sqrt(lPl/obsCount);
    if(dim == 0)
      return sigmaFloat;

    // apriori float ambiguities x0
    // ----------------------------
    Vector x0(dim);
    Vector isInteger(dim);
    for(auto ambi : getAmbiguities())
      if(ambi->index && ambi->value.size())
      {
        const UInt index = normalEquationInfo.index(ambi->index) - index0;
        copy(ambi->value, x0.row(index, ambi->value.rows()));
        if(ambi->isInteger)
          isInteger.row(index, ambi->value.rows()).fill(1);
      }
    Parallel::reduceSum(x0,        0, normalEquationInfo.comm);
    Parallel::reduceSum(isInteger, 0, normalEquationInfo.comm);
    Parallel::broadCast(isInteger, 0, normalEquationInfo.comm);

    Vector isNotResolved(dim, 1.);
    if(!isStrictlyZero(isInteger))
    {
      // right hand side and diagonal
      Vector rhs(dim);
      Vector diagonal(dim);
      for(UInt i=block0; i<normals.blockCount(); i++)
      {
        copy(n.at(i), rhs.row(normals.blockIndex(i)-index0, normals.blockSize(i)));
        if(normals.isMyRank(i,i))
        {
          Matrix &N = normals.N(i,i);
          for(UInt k=0; k<N.rows(); k++)
            diagonal(normals.blockIndex(i)-index0+k) = N(k,k);
        }
      }
      Parallel::reduceSum(diagonal, 0, normalEquationInfo.comm);
      Parallel::broadCast(diagonal, 0, normalEquationInfo.comm);

      // sort float ambiguites to begin and accurate integer ambiguities to end
      // ----------------------------------------------------------------------
      GnssLambda::Transformation Z(dim);
      std::vector<UInt> index(dim);
      std::iota(index.begin(), index.end(), 0);
      for(UInt i=0; i<dim; i++)
        if(isInteger(i))
        {
          UInt idxMin = i;
          for(UInt k=i+1; k<dim; k++)
            if(!isInteger(k) || (diagonal(k) < diagonal(idxMin)))
              idxMin = k;
          if(i == idxMin)
            continue;
          std::swap(diagonal.at(i), diagonal.at(idxMin));
          std::swap(isInteger(i),   isInteger(idxMin));
          std::swap(index.at(i),    index.at(idxMin));
          Z.swap(i, idxMin);
        }

      logStatus<<"- reorder normals of integer part"<<Log::endl;
      const UInt countInteger = static_cast<UInt>(sum(isInteger));
      const UInt startInteger = dim - countInteger;
      MatrixDistributed normalsAmbi = normals;
      normalsAmbi.eraseBlocks(0, block0);
      normalsAmbi.reorder(index, {0, startInteger, dim}, [](UInt, UInt, UInt){return 0;});

      // perform cholesky for float part of the normals
      // ----------------------------------------------
      logStatus<<"- cholesky normals of integer part"<<Log::endl;
      if(startInteger > 0)
        normalsAmbi.cholesky(FALSE/*timing*/, 0, 1, TRUE/*collect*/);
      if(normalsAmbi.isMyRank(1,1))
        GnssLambda::choleskyReversePivot(normalsAmbi.N(1,1), Z, TRUE/*timing*/);  // Z is now only valid at master

      // float solution
      // --------------
      rhs = Z.transform(rhs); // only valid at master
      x0  = Z.transform(x0);  // only valid at master
      normalsAmbi.triangularTransSolve(rhs);
      Vector dxFloat = rhs;
      sigmaFloat = std::sqrt((lPl-quadsum(rhs))/(obsCount-dim));
      normalsAmbi.triangularSolve(dxFloat);

      // inverse of cholesky
      // -------------------
      normalsAmbi.choleskyInverse(FALSE/*timing*/, 1, normalsAmbi.blockCount()-1);

      // xInteger = LAMBDA(W, xFloat)
      // ============================
      Vector dxInteger(dim);
      if(Parallel::isMaster(normalEquationInfo.comm))
      {
        // decorrelate integer ambiguities (may take a while)
        // -------------------------------------------------
        logStatus<<"- decorrelate ambiguities (may take a while)"<<Log::endl;
        Matrix &W(normalsAmbi.N(1,1));
        GnssLambda::Transformation T(W.rows());
        const Vector d = GnssLambda::choleskyTransform(W, T, TRUE/*timing*/);

        // check for NANs
        for(UInt i=0; i<W.rows(); i++)
          for(UInt k=i; k<W.rows(); k++)
            if(std::isnan(W(i,k)))
              throw(Exception("Cholesky matrix contains NANs"));
        for(UInt i=0; i<d.rows(); i++)
          if(std::isnan(d(i)))
            throw(Exception("diagonal vector d contains NANs"));

        // resolve ambiguities
        // -------------------
        logStatus<<"- search integer (blocked algorithm)"<<Log::endl;
        Vector notResolved;
        Double sigma;
        Vector xInteger = searchInteger(T.transform((dxFloat+x0).row(startInteger, countInteger)), W, d, notResolved, sigma);
        logInfo<<"  Ambiguity fixing norm = "<<sigma%"%3.2f"s<<Log::endl;
        logInfo<<"  "<<sum(notResolved)%"%4i of "s<<notResolved.rows()%"%i transformed ambiguities not resolved ("s<<(100.*sum(notResolved)/notResolved.rows())%"%.1f%%)"s<<Log::endl;
        logInfo<<"  "<<sum(T.distributeBack(notResolved))%"%4i remaining float ambiguities"s<<Log::endl;

        copy(T.distributeBack(notResolved), isNotResolved.row(startInteger, countInteger));
        copy(T.transformBack(xInteger) - x0.row(startInteger, countInteger), dxInteger.row(startInteger, countInteger));
        for(UInt i=0; i<dxInteger.rows(); i++)
          if(isNotResolved(i))
            dxInteger(i) = 0.;
      } // if(isMaster)
      // =================================================

      dxInteger     = Z.transformBack(dxInteger);
      isNotResolved = Z.distributeBack(isNotResolved);
      Parallel::broadCast(dxInteger,     0, normalEquationInfo.comm);
      Parallel::broadCast(isNotResolved, 0, normalEquationInfo.comm);

      // mark resolved ambiguities
      for(auto ambi : getAmbiguities())
        if(ambi->index && ambi->value.size() && ambi->isInteger)
        {
          ambi->resolved = Vector(ambi->value.rows(), 1.);
          const UInt index = normalEquationInfo.index(ambi->index) - index0;
          for(UInt i=0; i<ambi->resolved.rows(); i++)
            if(isNotResolved(index+i))
              ambi->resolved(i) = 0;
        }

      // Remove solved integer ambiguities from the right hand side n -= N*(x_integer-x0)
      // --------------------------------------------------------------------------------
      Vector dn(dxInteger.rows());
      for(UInt i=block0; i<normals.blockCount(); i++)
        if(normals.isMyRank(i,i))
          matMult(1., normals.N(i,i), dxInteger.row(normals.blockIndex(i)-index0, normals.blockSize(i)), dn.row(normals.blockIndex(i)-index0, normals.blockSize(i)));
      for(UInt i=block0; i<normals.blockCount(); i++)
        for(UInt k=i+1; k<normals.blockCount(); k++)
          if(normals.isMyRank(i,k))
          {
            matMult(1., normals.N(i,k),         dxInteger.row(normals.blockIndex(k)-index0, normals.blockSize(k)), dn.row(normals.blockIndex(i)-index0, normals.blockSize(i)));
            matMult(1., normals.N(i,k).trans(), dxInteger.row(normals.blockIndex(i)-index0, normals.blockSize(i)), dn.row(normals.blockIndex(k)-index0, normals.blockSize(k)));
          }
      Parallel::reduceSum(dn, 0, normalEquationInfo.comm);
      if(Parallel::isMaster(normalEquationInfo.comm))
        for(UInt i=block0; i<normals.blockCount(); i++)
        {
          lPl     += inner(dxInteger.row(normals.blockIndex(i)-index0, normals.blockSize(i)), dn.row(normals.blockIndex(i)-index0, normals.blockSize(i)))
                     -2.*inner(dxInteger.row(normals.blockIndex(i)-index0, normals.blockSize(i)), n.at(i));
          n.at(i) -= dn.row(normals.blockIndex(i)-index0, normals.blockSize(i));
        }

      // Integer ambiguities are estimated indepedently with high weight
      // ---------------------------------------------------------------
      const Double weightInteger = std::pow(1./0.001, 2); // [1/cycles^2]
      for(UInt i=block0; i<normals.blockCount(); i++)
        for(UInt z=0; z<normals.blockSize(i); z++)
          if(isNotResolved(normals.blockIndex(i)+z-index0) == 0.)
          {
            for(UInt k=block0; k<=i; k++)
              if(normals.isMyRank(k,i))
                normals.N(k,i).column(z).setNull(); // set column 0
            for(UInt k=i; k<normals.blockCount(); k++)
              if(normals.isMyRank(i,k))
                normals.N(i,k).row(z).setNull();    // set row 0
            if(normals.isMyRank(i,i))
              normals.N(i,i)(z,z) = weightInteger;
            if(Parallel::isMaster(normalEquationInfo.comm))
              n.at(i)(z,0) = weightInteger * dxInteger(normals.blockIndex(i)+z-index0);
          }
    } // if(some resolvable)

    // solve (forward step)
    // --------------------
    normals.cholesky(FALSE/*timing*/, block0, normals.blockCount()-block0, TRUE/*collect*/);
    normals.triangularTransSolve(n, block0, normals.blockCount()-block0, TRUE/*collect*/); // forward step
    if(Parallel::isMaster(normalEquationInfo.comm))
      for(UInt i=block0; i<normals.blockCount(); i++)
        for(UInt z=0; z<normals.blockSize(i); z++)
          if(isNotResolved(normals.blockIndex(i)+z-index0))
          {
            lPl      -= std::pow(n.at(i)(z,0), 2); // lPl = lPl - n1' N1^(-1) n1
            obsCount -= 1;
          }

    return sigmaFloat;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Double GnssParametrizationAmbiguities::updateParameter(const GnssNormalEquationInfo &normalEquationInfo, const_MatrixSliceRef x, const_MatrixSliceRef /*Wz*/)
{
  try
  {
    // float bias at transmitter
    // -------------------------
    Double maxChange = 0;
    Gnss::InfoParameterChange infoTrans("mm");
    for(auto para : paraTrans)
      if(para && para->index)
      {
        const Vector dBias = x.row(normalEquationInfo.index(para->index), para->types.size());
        for(UInt idType=0; idType<para->types.size(); idType++)
          para->trans->signalBias.biases.at(GnssType::index(para->trans->signalBias.types, para->types.at(idType))) += dBias(idType);
        for(UInt idType=0; idType<para->types.size(); idType++)
          if(infoTrans.update(1e3*dBias(idType)))
            infoTrans.info = "phase bias transmitter ("+para->types.at(idType).str()+")";
      }
    infoTrans.synchronizeAndPrint(normalEquationInfo.comm, 1e-3, maxChange);

    // float bias at receiver
    // ----------------------
    Gnss::InfoParameterChange infoRecv("mm");
    for(auto para : paraRecv)
      if(para && para->index)
      {
        const Vector dBias = para->Bias * x.row(normalEquationInfo.index(para->index), para->Bias.columns());
        for(UInt idType=0; idType<dBias.size(); idType++)
          para->recv->signalBias.biases.at(GnssType::index(para->recv->signalBias.types, para->types.at(idType))) += dBias(idType);
        for(UInt idType=0; idType<para->types.size(); idType++)
          if(infoRecv.update(1e3*dBias(idType)))
            infoRecv.info = "phase bias receiver ("+para->recv->name()+", "+para->types.at(idType).str()+")";
      }
    infoRecv.synchronizeAndPrint(normalEquationInfo.comm, 1e-3, maxChange);

    // ambiguities
    // -----------
    Gnss::InfoParameterChange info("cyc");
    for(auto ambi : getAmbiguities())
      if(ambi->index)
      {
        const Vector dx = x.row(normalEquationInfo.index(ambi->index), ambi->value.rows());
        ambi->value += dx;
        for(UInt k=0; k<dx.rows(); k++)
          if(info.update(dx(k)))
            info.info = "ambiguity ("+normalEquationInfo.parameterNames().at(normalEquationInfo.index(ambi->index)+k).str()+")";
      }
    info.synchronizeAndPrint(normalEquationInfo.comm, 0.19, maxChange); // cycles -> meter

    // remove integer part from observations
    // -------------------------------------
    for(auto ambi : getAmbiguities())
      if(ambi->value.size())
      {
        Vector x(ambi->value.rows());
        for(UInt i=0; i<x.rows(); i++)
          x(i) = std::round(ambi->value(i));
        ambi->value -= x;
        ambi->track->removeAmbiguitiesFromObservations(ambi->types, Vector(ambi->T*x));
      }

    // remove resolved ambiguities
    // ---------------------------
    for(auto ambi : getAmbiguities())
      if(ambi->resolved.size())
      {
        const UInt remove = static_cast<UInt>(sum(ambi->resolved));
        Matrix T(ambi->T.rows(), ambi->T.columns() - remove);
        Vector value(ambi->value.rows() - remove);
        UInt idx = 0;
        for(UInt i=0; i<ambi->value.rows(); i++)
          if(ambi->resolved(i) == 0.)
          {
            copy(ambi->T.column(i), T.column(idx));
            value(idx) = ambi->value(i);
            idx++;
          }
        ambi->T        = T;
        ambi->value    = value;
        ambi->resolved = Vector();
      }

    return maxChange;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/

Vector GnssParametrizationAmbiguities::Ambiguity::ambiguities(const std::vector<GnssType> &types) const
{
  try
  {
    Vector value(types.size());
    UInt idx;
    if(this->value.size())
      for(UInt idType=0; idType<types.size(); idType++)
        if(types.at(idType).isInList(this->types, idx))
          value(idType) = (T.row(idx) * this->value)(0,0);
    return value;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
