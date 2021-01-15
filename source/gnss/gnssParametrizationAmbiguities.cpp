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

#define DOCSTRING_GnssParametrizationAmbiguities

#include "base/import.h"
#include "parallel/matrixDistributed.h"
#include "parser/dataVariables.h"
#include "config/config.h"
#include "config/configRegister.h"
#include "inputOutput/logging.h"
#include "files/fileMatrix.h"
#include "files/fileParameterName.h"
#include "gnss/gnss.h"
#include "gnss/gnssObservation.h"
#include "gnss/gnssDesignMatrix.h"
#include "gnss/gnssTransmitter.h"
#include "gnss/gnssReceiver.h"
#include "gnss/gnssParametrizationAmbiguities.h"

/***********************************************/

GROOPS_REGISTER_CLASS_WITHOUT_SUBS(GnssParametrizationAmbiguities, "gnssParametrizationAmbiguitiesType")
GROOPS_READCONFIG_CLASS(GnssParametrizationAmbiguities, "gnssParametrizationAmbiguitiesType")

/***********************************************/

GnssParametrizationAmbiguities::GnssParametrizationAmbiguities(Config &config, const std::string &name)
{
  try
  {
    std::string choice;
    sigmaMaxResolve  = 1e99;
    maxSearchSteps   = MAX_UINT;

    readConfigSequence(config, name, Config::MUSTSET, "", "");
    readConfig(config, "outputfileAmbiguities", fileNameAmbiguities, Config::OPTIONAL, "",    "resolved ambiguities, station name is appended to file name for single receivers");
    readConfig(config, "sigmaMaxResolve",       sigmaMaxResolve,     Config::OPTIONAL, "0.2", "max. allowed std. dev. of ambiguity to resolve [cycles]");
    readConfig(config, "searchBlockSize",       searchBlockSize,     Config::DEFAULT,  "200", "block size for blocked integer search");
    readConfig(config, "maxSearchSteps",        maxSearchSteps,      Config::OPTIONAL, "200000000", "max. steps of integer search for each block");
    if(readConfigChoice(config, "incompleteAction", choice, Config::MUSTSET, "shrinkBlockSize", "if not all solutions tested after maxSearchSteps"))
    {
      if(readConfigChoiceElement(config, "stop",            choice, "stop searching, ambiguities remain float in this block")) incompleteAction = IncompleteAction::STOP;
      if(readConfigChoiceElement(config, "resolve",         choice, "use best integer solution found so far"))                 incompleteAction = IncompleteAction::IGNORE;
      if(readConfigChoiceElement(config, "shrinkBlockSize", choice, "try again with half block size"))                         incompleteAction = IncompleteAction::SHRINKBLOCKSIZE;
      if(readConfigChoiceElement(config, "throwException",  choice, "stop and throw an exception"))                            incompleteAction = IncompleteAction::EXCEPTION;
      endChoice(config);
    }
    endSequence(config);
    if(isCreateSchema(config)) return;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssParametrizationAmbiguities::initIntervalLate(Gnss::AnalysisType /*analysisType*/, const std::vector<Time> &/*times*/, const Time &/*timeMargin*/, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    ambiguities.clear();
    for(auto &recv : gnss().receiver)
      if(recv->useable())
        for(auto &track : recv->track)
        {
          AmbiguityPtr ambiguity  = std::make_shared<Ambiguity>();
          ambiguity->track        = track.get();
          ambiguity->types        = Gnss::replaceCompositeSignals(track->types);
          ambiguity->T            = phaseDecorrelation(ambiguity->types, recv->wavelengthFactor());
          ambiguity->value        = Vector(ambiguity->T.columns());
          ambiguity->isInteger    = FALSE;
          track->ambiguity = ambiguity.get();
          ambiguities.push_back(ambiguity);
        }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

// to share information between processes
class GnssParametrizationAmbiguities::AmbiguityInfo
{
public:
  UInt                  idProcess, idAmbi;
  UInt                  idRecv, idTrans;
  UInt                  idEpochStart, idEpochEnd;
  UInt                  parameterCount;
  std::vector<GnssType> typesAmbiguity, typesTrack;
  std::vector<GnssType> typesNew;

  AmbiguityInfo() {}

  void save(OutArchive &oa) const
  {
    oa<<nameValue("idProcess",      idProcess);
    oa<<nameValue("idAmbi",         idAmbi);
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
    ia>>nameValue("idProcess",      idProcess);
    ia>>nameValue("idAmbi",         idAmbi);
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

void GnssParametrizationAmbiguities::initParameter(Gnss::NormalEquationInfo &normalEquationInfo)
{
  try
  {
    // reset parameter indices
    // -----------------------
    for(auto ambi : ambiguities)
      ambi->parameterIndex = Gnss::ParameterIndex();
    receiverName.clear();

    // Code only -> no ambiguities
    if(!((normalEquationInfo.analysisType & Gnss::ANALYSIS_PHASE) && (normalEquationInfo.estimationType & Gnss::NormalEquationInfo::ESTIMATE_AMBIGUITIES)))
      return;

    // receiver name of a single receiver solution
    if(std::count(normalEquationInfo.estimateReceiver.begin(), normalEquationInfo.estimateReceiver.end(), TRUE) == 1) // only single station
      receiverName = gnss().receiver.at(std::distance(normalEquationInfo.estimateReceiver.begin(), std::find(normalEquationInfo.estimateReceiver.begin(), normalEquationInfo.estimateReceiver.end(), TRUE)))->name();

    // distribute ambiguities to all nodes
    // -----------------------------------
    std::vector<AmbiguityInfo> ambiguityInfos;
    for(UInt idProcess=0; idProcess<Parallel::size(normalEquationInfo.comm); idProcess++)
    {
      std::vector<AmbiguityInfo> ambiguityInfosLocal;
      if(Parallel::myRank(normalEquationInfo.comm) == idProcess)
        for(auto ambi : ambiguities)
          if(ambi->track->receiver->useable() && normalEquationInfo.estimateReceiver.at(ambi->track->receiver->idRecv()))
          {
            AmbiguityInfo info;
            info.idProcess      = idProcess;
            info.idAmbi         = std::distance(ambiguities.begin(), std::find(ambiguities.begin(), ambiguities.end(), ambi));
            info.idRecv         = ambi->track->receiver->idRecv();
            info.idTrans        = ambi->track->transmitter->idTrans();
            info.idEpochStart   = ambi->track->idEpochStart;
            info.idEpochEnd     = ambi->track->idEpochEnd;
            info.parameterCount = ambi->value.rows();
            info.typesAmbiguity = ambi->types;
            info.typesTrack     = Gnss::replaceCompositeSignals(ambi->track->types);
            ambiguityInfosLocal.push_back(info);
          }
      // synchronize information between processes
      Parallel::broadCast(ambiguityInfosLocal, idProcess, normalEquationInfo.comm);
      ambiguityInfos.insert(ambiguityInfos.end(), ambiguityInfosLocal.begin(), ambiguityInfosLocal.end());
    }

    // sort with decreasing track length
    std::stable_sort(ambiguityInfos.begin(), ambiguityInfos.end(), [](const AmbiguityInfo &a, const AmbiguityInfo &b)
                    {
                      if((a.typesTrack.size()-a.typesAmbiguity.size()) == (b.typesTrack.size()-b.typesAmbiguity.size()))
                        return (a.idEpochEnd-a.idEpochStart) > (b.idEpochEnd-b.idEpochStart);
                      return (a.typesTrack.size()-a.typesAmbiguity.size()) > (b.typesTrack.size()-b.typesAmbiguity.size());
                    });

    // determine integer ambiguities (skip zero/single diff float ambiguities)
    // -----------------------------------------------------------------------
    std::vector<std::vector<GnssType>> typesTransmitter(gnss().transmitter.size()); // float ambiguity at transmitter
    std::vector<std::vector<GnssType>> typesReceiver(gnss().receiver.size());       // float ambiguity at receiver
    std::vector<GnssType>              typesBoth;                                   // float ambiguity at receiver and transmitter
    Bool restart = TRUE;
    while(restart)
    {
      restart = FALSE;
      for(auto &info : ambiguityInfos)
      {
        for(GnssType &type : info.typesTrack)
          if(type == GnssType::PHASE)
          {
            type &= ~GnssType::PRN; // ambiguity depends not on PRN
            const Bool isFloatTrans = gnss().transmitter.at(info.idTrans)->isPhaseBiasEstimated(normalEquationInfo) && (GnssType::index(typesTransmitter.at(info.idTrans), type) == NULLINDEX);
            const Bool isFloatRecv  = gnss().receiver.at(info.idRecv)->isPhaseBiasEstimated(normalEquationInfo)     && (GnssType::index(typesReceiver.at(info.idRecv),     type) == NULLINDEX);

            // only one zero difference ambiguity is allowed for each type
            if(isFloatTrans && isFloatRecv && (GnssType::index(typesBoth, type) != NULLINDEX))
              continue; // try again later when new single float ambiguities are set up

            if(isFloatTrans && isFloatRecv)
              typesBoth.push_back(type);                         // zero difference
            if(isFloatTrans)
              typesTransmitter.at(info.idTrans).push_back(type); // float transmitter single difference
            if(isFloatRecv)
              typesReceiver.at(info.idRecv).push_back(type);     // float receiver single difference
            if(!isFloatTrans && !isFloatRecv)
              info.typesNew.push_back(type);                     // integer double difference
            type = GnssType::RANGE;                              // to indicate type is already processed
            restart = (isFloatTrans || isFloatRecv);             // new float ambiguity -> restart searching
            if(restart)
              break;
          }
        if(restart)
          break;
      } // for(ambiguityInfos)
    } // while(restart)

    for(auto &info : ambiguityInfos)
      if(GnssType::index(info.typesTrack, GnssType::PHASE) != NULLINDEX)
        throw(Exception("Cannot setup ambiguities. Network is separated into independent parts?"));

    // calculate parameter names and indicies
    // --------------------------------------
    // sort in temporal order
    std::stable_sort(ambiguityInfos.begin(), ambiguityInfos.end(), [](const AmbiguityInfo &a, const AmbiguityInfo &b) {return (a.idEpochStart+a.idEpochEnd) < (b.idEpochStart+b.idEpochEnd);});

    for(auto &info : ambiguityInfos)
    {
      if(info.typesNew != info.typesAmbiguity)
        info.parameterCount = info.typesNew.size();

      // parameter names
      std::vector<ParameterName> parameterNames;
      ParameterName name(gnss().receiver.at(info.idRecv)->name()+"."+gnss().transmitter.at(info.idTrans)->name(),
                         "ambiguity", "", gnss().times.at(info.idEpochStart), gnss().times.at(info.idEpochEnd));
      std::string typeStr;
      for(auto type : info.typesAmbiguity)
        typeStr += type.str().substr(0,2);
      for(UInt i=0; i<info.parameterCount; i++)
        parameterNames.push_back(ParameterName(name.object, name.type+(i+1)%"%iof"s+info.parameterCount%"%i("s+typeStr+")", name.temporal, name.interval));
      Gnss::ParameterIndex parameterIndex = normalEquationInfo.parameterNamesAmbiguity((info.idEpochStart+info.idEpochEnd)/2, info.idRecv, info.idTrans, parameterNames);

      if(Parallel::myRank(normalEquationInfo.comm) == info.idProcess)
      {
        AmbiguityPtr ambi  = ambiguities.at(info.idAmbi);
        ambi->isInteger    = ambi->track->receiver->supportsIntegerAmbiguities(normalEquationInfo) &&
                             ambi->track->transmitter->supportsIntegerAmbiguities(normalEquationInfo);
        if(info.typesNew != ambi->types)
        {
          ambi->types = info.typesNew;
          ambi->T     = phaseDecorrelation(ambi->types, gnss().receiver.at(info.idRecv)->wavelengthFactor());
          ambi->value = Vector(ambi->T.columns());
        }
        ambi->parameterIndex = parameterIndex;
      }
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssParametrizationAmbiguities::aprioriParameter(const Gnss::NormalEquationInfo &normalEquationInfo, MatrixSliceRef x0) const
{
  try
  {
    for(auto ambi : ambiguities)
      if(ambi->parameterIndex)
        copy(ambi->value, x0.row(normalEquationInfo.index(ambi->parameterIndex), ambi->value.rows()));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Bool GnssParametrizationAmbiguities::isDesignMatrix(const Gnss::NormalEquationInfo &/*normalEquationInfo*/, UInt idRecv, UInt idTrans, UInt idEpoch) const
{
  try
  {
    const Gnss::Track *track = gnss().receiver.at(idRecv)->observation(idTrans, idEpoch)->track;
    if(!track || !track->ambiguity)
      return FALSE;
    const Ambiguity *ambiguity = dynamic_cast<Ambiguity*>(track->ambiguity);
    return (ambiguity->parameterIndex) && ambiguity->value.size();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssParametrizationAmbiguities::designMatrix(const Gnss::NormalEquationInfo &/*normalEquationInfo*/, const Gnss::ObservationEquation &eqn, Gnss::DesignMatrix &A) const
{
  try
  {
    if((!eqn.track) || (!eqn.track->ambiguity))
      return;
    const Ambiguity *ambiguity = dynamic_cast<Ambiguity*>(eqn.track->ambiguity);
    if(ambiguity->parameterIndex)
    {
      MatrixSlice Design(A.column(ambiguity->parameterIndex));
      for(UInt idType=0; idType<eqn.typesTransmitted.size(); idType++)
      {
        const UInt idx = GnssType::index(ambiguity->types, eqn.typesTransmitted.at(idType));
        if(idx != NULLINDEX)
          matMult(1., eqn.A.column(Gnss::ObservationEquation::idxUnit + eqn.types.size() + idType), ambiguity->T.row(idx), Design);
      }
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Double GnssParametrizationAmbiguities::updateParameter(const Gnss::NormalEquationInfo &normalEquationInfo, const_MatrixSliceRef x, const_MatrixSliceRef /*Wz*/, Bool printStatistics)
{
  try
  {
    Double      maxChange = 0;
    std::string infoMaxChange;

    for(auto ambi : ambiguities)
      if(ambi->parameterIndex)
      {
        const Vector dx = x.row(normalEquationInfo.index(ambi->parameterIndex), ambi->value.rows());
        ambi->value += dx;

        for(UInt k=0; k<dx.rows(); k++)
          if(fabs(dx(k)) > maxChange)
          {
            maxChange = fabs(dx(k));
            infoMaxChange = "  "+normalEquationInfo.parameterNames().at(normalEquationInfo.index(ambi->parameterIndex)+k).str()+" change = "+dx(k)%"%8.3f"s+" cycles";
          }
      }

    // remove integer part from observations
    // -------------------------------------
    for(auto ambi : ambiguities)
      if(ambi->value.size())
      {
        Vector x(ambi->value.rows());
        for(UInt i=0; i<x.rows(); i++)
          x(i) = std::round(ambi->value(i));
        ambi->value -= x;
        ambi->track->removeAmbiguitiesFromObservations(ambi->types, Vector(ambi->T*x));
      }

    // adjust resolved ambiguities
    // ---------------------------
    for(auto ambi : ambiguities)
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

    Gnss::checkMaxChange(maxChange, infoMaxChange, printStatistics, normalEquationInfo.comm);
    return maxChange * 0.19; // cycles -> meter
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Double GnssParametrizationAmbiguities::ambiguityResolve(const Gnss::NormalEquationInfo &normalEquationInfo, MatrixDistributed &normals, std::vector<Matrix> &n, Double &lPl, UInt &obsCount, Bool dryRun)
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
    for(auto ambi : ambiguities)
      if(ambi->parameterIndex && ambi->value.size())
      {
        const UInt index = normalEquationInfo.index(ambi->parameterIndex) - index0;
        copy(ambi->value, x0.row(index, ambi->value.rows()));
        if(ambi->isInteger)
          isInteger.row(index, ambi->value.rows()).fill(1);
      }
    Parallel::reduceSum(x0, 0, normalEquationInfo.comm);
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
      Transformation Z(dim);
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
        choleskyReversePivot(normalsAmbi.N(1,1), Z);  // Z is now only valid at master

      // float solution
      // --------------
      rhs = Z.transform(rhs); // only valid at master
      x0  = Z.transform(x0);  // only valid at master
      normalsAmbi.triangularTransSolve(rhs);
      Vector dxFloat = rhs;
      sigmaFloat = sqrt((lPl-quadsum(rhs))/(obsCount-dim));
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
        Transformation T(W.rows());
        const Vector d = choleskyTransform(W, T);

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
        Vector xInteger;
        Vector notResolved;
        Matrix solutionSteps;
        Double sigma;
        searchIntegerBlocked(T.transform((dxFloat+x0).row(startInteger, countInteger)), W, d,
                             xInteger, notResolved, solutionSteps, sigma);
        logInfo<<"  Ambiguity fixing norm = "<<sigma%"%3.2f"s<<Log::endl;
        logInfo<<"  "<<sum(notResolved)%"%4i of "s<<notResolved.rows()%"%i transformed ambiguities not resolved ("s<<(100.*sum(notResolved)/notResolved.rows())%"%.1f%%)"s<<Log::endl;
        logInfo<<"  "<<sum(T.distributeBack(notResolved))%"%4i remaining float ambiguities"s<<Log::endl;

        copy(T.distributeBack(notResolved), isNotResolved.row(startInteger, countInteger));
        copy(T.transformBack(xInteger) - x0.row(startInteger, countInteger), dxInteger.row(startInteger, countInteger));
        for(UInt i=0; i<dxInteger.rows(); i++)
          if(isNotResolved(i))
            dxInteger(i) = 0.;

        if(!fileNameAmbiguities.empty())
        {
          this->solutionSteps.push_back(solutionSteps);
          this->solutionNames.push_back(receiverName);
        }
      } // if(isMaster)
      // =================================================

      dxInteger     = Z.transformBack(dxInteger);
      isNotResolved = Z.distributeBack(isNotResolved);
      Parallel::broadCast(dxInteger,     0, normalEquationInfo.comm);
      Parallel::broadCast(isNotResolved, 0, normalEquationInfo.comm);

      // mark resolved ambiguities
      if(!dryRun)
        for(auto ambi : ambiguities)
          if(ambi->parameterIndex && ambi->value.size() && ambi->isInteger)
          {
            ambi->resolved = Vector(ambi->value.rows(), 1.);
            const UInt index = normalEquationInfo.index(ambi->parameterIndex) - index0;
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

void GnssParametrizationAmbiguities::writeResults(const Gnss::NormalEquationInfo &/*normalEquationInfo*/, const std::string &suffix)
{
  try
  {
    VariableList fileNameVariableList;
    addTimeVariables(fileNameVariableList);
    evaluateTimeVariables(0, gnss().times.at(0), gnss().times.back(), fileNameVariableList);

    if(!fileNameAmbiguities.empty())
    {
      logStatus<<"write ambiguities to file(s) <"<<fileNameAmbiguities(fileNameVariableList).appendBaseName(suffix)<<">"<<Log::endl;
      for(UInt i=0; i<solutionSteps.size(); i++)
        writeFileMatrix(fileNameAmbiguities(fileNameVariableList).appendBaseName((solutionNames.at(i).empty() ? "" : "."+solutionNames.at(i))+suffix), solutionSteps.at(i));
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/

Matrix GnssParametrizationAmbiguities::phaseDecorrelation(const std::vector<GnssType> &types, Double wavelengthFactor)
{
  try
  {
    if(types.size() == 0)
      return Matrix();
    if(types.size() == 1)
      return Matrix(1, 1, wavelengthFactor*LIGHT_VELOCITY/types.at(0).frequency());

    // design matrix
    // assume for every phase observation an additional range observation
    // parameters: range, TEC, ambiguities
    const UInt dim = types.size();
    Matrix A(2*dim, 2+dim);
    for(UInt i=0; i<dim; i++)
    {
      // phase observations:
      A(i, 0)   = 1.;                                                      // range
      A(i, 1)   = -Ionosphere::Ap/pow(types.at(i).frequency(), 2);         // TEC
      A(i, 2+i) = wavelengthFactor*LIGHT_VELOCITY/types.at(i).frequency(); // ambiguity
      // range observations (100 times less accurate):
      A(i+dim, 0) = 1./100.;       // range
      A(i+dim, 1) = -A(i, 1)/100.; // TEC
    }

    // solve & deccorelate
    QR_decomposition(A);
    inverse(A.slice(2, 2, dim, dim));
    Transformation Z(dim);
    choleskyTransform(A.slice(2, 2, dim, dim), Z);

    // Transformation matrix cycles -> meter
    Matrix T = Z.transformBack(identityMatrix(types.size()));
    for(UInt idType=0; idType<types.size(); idType++)
      T.row(idType) *= wavelengthFactor*LIGHT_VELOCITY/types.at(idType).frequency();

    return T;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssParametrizationAmbiguities::choleskyReversePivot(Matrix &N, Transformation &Z)
{
  try
  {
    Vector tmp(N.rows());
    logTimerStart;
    for(UInt i=0; i<N.rows(); i++)
    {
      logTimerLoop(i, N.rows());

      // find minimium
      UInt   k    = i;
      Double minN = N(i,i)-tmp(i);
      for(UInt j=i+1; j<N.rows(); j++)
        if(N(j,j)-tmp(j) < minN)
        {
          k    = j;
          minN = N(k,k)-tmp(k);
        }
      // swap
      if(i != k)
      {
        Z.swap(i, k);
        std::swap(tmp(i), tmp(k));
        std::swap(N(i,i), N(k,k));
        if(i>0)          swap(N.slice(0, i,i, 1),               N.slice(0,   k,   i,     1));
        if(k<N.rows()-1) swap(N.slice(i, k+1, 1, N.rows()-k-1), N.slice(k,   k+1, 1,     N.rows()-k-1));
        if(k>i+1)        swap(N.slice(i, i+1, 1, k-i-1),        N.slice(i+1, k,   k-i-1, 1).trans());
      }
      // cholesky
      N(i,i) -= tmp(i);
      N(i,i)  = std::sqrt(N(i,i));
      if(i+1<N.rows())
      {
        if(i > 0)
          matMult(-1., N.slice(0, i, i, 1).trans(), N.slice(0, i+1, i, N.rows()-1-i), N.slice(i, i+1, 1, N.rows()-1-i));
        N.slice(i,i+1,1,N.rows()-1-i) *= 1./N(i,i);
        for(UInt k=i+1; k<N.rows(); k++)
          tmp(k) += std::pow(N(i,k), 2);
      }
    }
    logTimerLoopEnd(N.rows());
    N.setType(Matrix::TRIANGULAR);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

// swap element i and i+1 of cholesky of covariance matrix
inline void GnssParametrizationAmbiguities::choleskySwap(UInt i, Double delta, MatrixSliceRef W, std::vector<Double> &d, Transformation &Z)
{
  try
  {
    const UInt   dim    = W.rows();
    const Double eta    = d[i]/delta;
    const Double lambda = d[i+1]*W(i,i+1)/delta;
    d[i]   = eta*d[i+1];
    d[i+1] = delta;

    if(i>0)
      copy(W.slice(0,i,i,2)*Matrix({{-W(i,i+1),  eta}, {1., lambda}}), W.slice(0,i,i,2));
    W(i,i+1) = lambda;
    if(i+2<W.rows())
      swap(W.slice(i,i+2,1,dim-2-i), W.slice(i+1,i+2,1,dim-2-i));

    // store applied transformation
    Z.swap(i, i+1);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

// reduce element (i, k) (reduce row k from i)
inline void GnssParametrizationAmbiguities::choleskyReduce(UInt i, UInt k, MatrixSliceRef W, Transformation &Z)
{
  try
  {
    const Double alpha = std::round(W(i,k));
    if(alpha == 0.)
      return;
    axpy(-alpha, W.slice(k,k,1,W.rows()-k), W.slice(i,k,1,W.rows()-k));
    Z.reduce(alpha, k, i);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Vector GnssParametrizationAmbiguities::choleskyTransform(MatrixSliceRef W, Transformation &Z)
{
  try
  {
    const UInt dim = W.rows();

    // Covariance = W D W^T
    // D: diagonal matrix, W: upper triangular
    std::vector<Double> d(dim);
    for(UInt i=0; i<dim; i++)
      d[i] = W(i,i)*W(i,i);
    for(UInt i=0; i<dim; i++)
      W.column(i) *= 1./W(i,i);

    // decorrelate
    for(UInt i=dim; i-->0;)
      for(UInt k=i+1; k<dim; k++)
        choleskyReduce(i, k, W, Z);

    std::vector<Double> delta(dim-1);
    for(UInt i=0; i<dim-1; i++)
      delta[i] = d[i] + std::pow(W(i,i+1), 2) * d[i+1];

    std::multimap<Double, UInt> ratio;
    for(UInt i=0; i<dim-1; i++)
      if(delta[i] < d[i+1])
        ratio.insert(std::make_pair(delta[i]/d[i+1], i));

    if(dim>100) logTimerStart;
    UInt iter = 0;
    while(!ratio.empty())
    {
      if((dim>100) && (iter%1000 == 0)) logTimerLoop(iter, iter+1);
      iter++;

      // find maximum change
      const UInt k = ratio.begin()->second;

      // erase ratios that will be changed
      ratio.erase(ratio.begin());
      if((k > 0) && (delta[k-1] < d[k]))
      {
        auto iterpair = ratio.equal_range(delta[k-1]/d[k]);
        ratio.erase(std::find_if(iterpair.first, iterpair.second, [k](const std::pair<Double, UInt> &x){return x.second == k-1;}));
      }
      if((k+1 < dim-1) && (delta[k+1] < d[k+2]))
      {
        auto iterpair = ratio.equal_range(delta[k+1]/d[k+2]);
        ratio.erase(std::find_if(iterpair.first, iterpair.second, [k](const std::pair<Double, UInt> &x){return x.second == k+1;}));
      }

      // swap k and k+1
      choleskySwap(k, delta[k], W, d, Z);

      // decorrelate, update delta and ratio
      for(UInt i=std::min(k+2,dim-1); i-->std::max(k,UInt(1))-1;)
      {
        // mathematically only choleskyReduce(i, i+1, W, Z) is needed
        for(UInt k=i+1; k<dim; k++)
          choleskyReduce(i, k, W, Z);
        delta[i] = d[i] + std::pow(W(i,i+1), 2) * d[i+1];
        if(delta[i] < d[i+1])
          ratio.insert(std::make_pair(delta[i]/d[i+1], i));
      }
    }
    if(dim>100) logTimerLoopEnd(iter);

    // decorrelate rest of the triangle
    for(UInt i=dim; i-->0;)
      for(UInt k=i+1; k<dim; k++)
        choleskyReduce(i, k, W, Z);

    return d;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Bool GnssParametrizationAmbiguities::searchInteger(const_MatrixSliceRef xFloat, MatrixSliceRef W, const_MatrixSliceRef d,
                                                   UInt maxSearchSteps, Vector &solution, Double &minNorm)
{
  try
  {
    const UInt dim = W.rows();
    std::vector<Double> xInt(dim, 0);
    std::vector<Double> norm(dim, 0);
    std::vector<Double> step(dim, 0);
    std::vector<Double> dx  (dim, 0);

    // store xBar in lower triangle of W
    // xBar vector of step i starts at W(dim-1-i, dim-1-i), reverse order
    MatrixSlice xBar(W);
    copy(xFloat, xBar.column(0));

    UInt i  = dim-1;
    xInt[i] = std::round(xBar(dim-1, dim-1-i));
    dx[i]   = xInt[i] - xBar(dim-1, dim-1-i);
    step[i] = (dx[i] < 0) ? 1 : -1;
    minNorm = 1e99;
    Double newNorm = dx[i]*dx[i]/d(i,0);

    UInt iter=0;
    for(; iter<maxSearchSteps; iter++)
    {
      // move down
      const UInt lastValidXBar = i;
      while((newNorm < minNorm) && (i-- > 0))
      {
        // compute new xBar
        for(UInt k=lastValidXBar; k-->i;)
          xBar(i+dim-1-k, dim-1-k) = xBar(i+dim-2-k, dim-2-k) + dx[k+1] * W(i,k+1);

        xInt[i]  = std::round(xBar(dim-1, dim-1-i));
        dx[i]    = xInt[i] - xBar(dim-1, dim-1-i);
        step[i]  = (dx[i] < 0) ? 1 : -1;
        norm[i]  = newNorm;
        newNorm += dx[i]*dx[i]/d(i,0);
      }

      // new solution?
      if(newNorm < minNorm)
      {
        i = 0;
        solution = Vector(xInt);
        minNorm  = newNorm;
        newNorm  = 2e99;
      }

      // move up
      while((newNorm >= minNorm) && (i++ < dim-1))
      {
        xInt[i] += step[i];
        dx[i]   += step[i];
        step[i]  = (step[i]>0) ? (-step[i]-1) : (-step[i]+1); // zig-zag search
        newNorm  = norm[i] + dx[i]*dx[i]/d(i,0);
      }

      if(newNorm >= minNorm)
        break;
    } // for(iter)

    // restore diagonal
    for(UInt i=0; i<dim; i++)
      W(i,i) = 1.;

    return (iter < maxSearchSteps);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssParametrizationAmbiguities::searchIntegerBlocked(const_MatrixSliceRef xFloat, MatrixSliceRef W, const_MatrixSliceRef d,
                                           Vector &xInt, Vector &isNotFixed, Matrix &solutionSteps, Double &sigma) const
{
  try
  {
    const UInt dim = W.rows();
    xInt           = Vector(dim);
    isNotFixed     = Vector(dim, 1.);
    solutionSteps  = Matrix();
    sigma          = 0.;

    // find values to be solved
    UInt minIndex = dim;
    while((minIndex>0) && (d(minIndex-1,0) < sigmaMaxResolve*sigmaMaxResolve))
      minIndex--;
    if(minIndex >= dim)
      return;
    UInt idxSolved = dim;

    // float solution
    std::vector<Vector> solutions(1, d);
    auto storeSolution = [&](Vector xBar, UInt idxSolved)
    {
      for(UInt i=dim; i-->idxSolved;)
        axpy(-xBar(i), W.slice(0,i,i,1), xBar.row(0,i));
      for(UInt i=0; i<idxSolved; i++)
        xBar(i) = std::fabs(std::remainder(xBar(i), 1));
      for(UInt i=idxSolved; i<dim; i++)
        xBar(i) = -std::fabs(xBar(i));
      solutions.push_back(xBar);
    };
    storeSolution(xFloat, idxSolved);

    UInt defaultBlockSize = std::min(dim-minIndex, searchBlockSize);
    UInt blockSize     = defaultBlockSize;
    UInt blockStart    = dim-blockSize;
    UInt blockStartOld = dim;
    UInt iter=0;
    logTimerStart;
    for(;;)
    {
      logTimerLoop(iter++, (dim-minIndex)/(defaultBlockSize/2));

      // compute xBar for found solution so far
      Vector xBar = xFloat-xInt;
      for(UInt i=dim; i-->blockStart+blockSize;)
        axpy(-xBar(i), W.slice(0,i,i,1), xBar.row(0,i));

      Vector dxInt;
      Double ePe;
      Bool completed = searchInteger(xBar.row(blockStart, blockSize),
                                     W.slice(blockStart, blockStart, blockSize, blockSize),
                                     d.row(blockStart, blockSize), maxSearchSteps, dxInt, ePe);
      axpy(1., dxInt, xInt.row(blockStart, blockSize));
      copy(Vector(blockSize, 0), isNotFixed.row(blockStart, blockSize));
      idxSolved = blockStart;

      if(!completed)
      {
        if(incompleteAction == IncompleteAction::STOP)
        {
          logWarning<<receiverName<<": searchInteger("<<dim<<").slice("<<blockStart<<", "<<blockSize<<"): cannot find a solution -> stop searching"<<Log::endl;
          idxSolved = blockStart+blockSize;
          break;
        }
        else if(incompleteAction == IncompleteAction::SHRINKBLOCKSIZE)
        {
          logWarning<<receiverName<<": searchInteger("<<dim<<").slice("<<blockStart<<", "<<blockSize<<"): cannot find a solution -> restart with smaller block size"<<Log::endl;
          defaultBlockSize = (defaultBlockSize+1)/2;
          blockStart += (blockSize+1)/2;
          blockSize  -= (blockSize+1)/2;
          iter--;
          continue;
        }
        else if(incompleteAction == IncompleteAction::EXCEPTION)
        {
          logWarning<<receiverName<<": searchInteger("<<dim<<").slice("<<blockStart<<", "<<blockSize<<"): cannot find a solution"<<Log::endl;
          throw(Exception("Ambiguity resolution failed."));
        }
        // incompleteAction == IncompleteAction::IGNORE)
      }

      storeSolution(xFloat-xInt, idxSolved);

      // check consistence with solution of old block, if not restart with larger block
      if(quadsum(dxInt.row(blockStartOld-blockStart, blockStart+blockSize-blockStartOld)))
      {
        logWarning<<receiverName<<": searchInteger("<<dim<<").slice("<<blockStart<<", "<<blockSize<<"): not consistent -> restart"<<Log::endl;
        blockStartOld = blockStart+blockSize;
        blockSize    += std::min(defaultBlockSize/2, dim-blockStart-blockSize);
        iter--;
        continue;
      }

      blockStartOld = blockStart;
      blockSize     = defaultBlockSize;
      if(blockStart == minIndex) // at beginning?
        break;
      blockStart = std::max(blockStart, minIndex+defaultBlockSize/2) - defaultBlockSize/2;
    } // for(blocks)
    logTimerLoopEnd(iter);

    xInt.row(0, idxSolved).fill(0);
    isNotFixed.row(0, idxSolved).fill(1);

    // compute norm
    Vector xBar = xFloat - xInt;
    triangularSolve(1., W, xBar);
    sigma = 0;
    for(UInt i=idxSolved; i<dim; i++)
      sigma += xBar(i)*xBar(i)/d(i,0);
    sigma = std::sqrt(sigma/(dim-idxSolved));

    // copy solutions into one matrix
    solutionSteps = Matrix(dim, solutions.size());
    for(UInt i=0; i<solutions.size(); i++)
      copy(solutions.at(i), solutionSteps.column(i));
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
    if(this->value.size())
      for(UInt idType=0; idType<types.size(); idType++)
      {
        const UInt i = GnssType::index(this->types, types[idType]);
        if(i != NULLINDEX)
          value(idType) = (T.row(i) * this->value)(0,0);
      }
    return value;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/

GnssParametrizationAmbiguities::Transformation::Transformation(UInt dim)
  : columnForward(dim), columnBackward(dim)
{
  try
  {
    // unity matrix in both directions
    for(UInt i=0; i<dim; i++)
    {
      columnForward[i].insert (std::make_pair(i, 1.0));
      columnBackward[i].insert(std::make_pair(i, 1.0));
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

// row(k) -= alpha * row(i)
inline void GnssParametrizationAmbiguities::Transformation::reduce(Double alpha, UInt i, UInt k)
{
  try
  {
    if(alpha == 0.)
      return;

    // --- lambda ---
    auto reduce = [](Double alpha, std::map<UInt, Double> &column1, std::map<UInt, Double> &column2)
    {
      for(auto &col : column1)
      {
        auto element = column2.find(col.first);
        if(element != column2.end())
        {
          element->second -= alpha * col.second;
          if(std::fabs(element->second) < 1e-9)
            column2.erase(element);
        }
        else
          column2.insert(std::make_pair(col.first, -alpha * col.second));
      }
    };
    // -----------

    reduce(+alpha, columnForward[i],  columnForward[k]);
    reduce(-alpha, columnBackward[k], columnBackward[i]);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline void GnssParametrizationAmbiguities::Transformation::swap(UInt i, UInt k)
{
  std::swap(columnForward.at(i),  columnForward.at(k));
  std::swap(columnBackward.at(i), columnBackward.at(k));
}

/***********************************************/

Matrix GnssParametrizationAmbiguities::Transformation::transform(const_MatrixSliceRef x) const
{
  try
  {
    if(x.rows() != columnForward.size())
      throw(Exception("Dimension error"));

    Matrix y(x.rows(), x.columns());
    for(UInt i=0; i<columnForward.size(); i++)
      for(const auto &col : columnForward[i])
        axpy(col.second, x.row(col.first), y.row(i));
    return y;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Matrix GnssParametrizationAmbiguities::Transformation::transformBack(const_MatrixSliceRef x) const
{
  try
  {
    Matrix y(x.rows(), x.columns());
    for(UInt i=0; i<columnBackward.size(); i++)
      for(const auto &col : columnBackward[i])
        axpy(col.second, x.row(i), y.row(col.first));
    return y;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Matrix GnssParametrizationAmbiguities::Transformation::distributeBack(const_MatrixSliceRef x) const
{
  try
  {
    Matrix y(x.rows(), x.columns());
    for(UInt i=0; i<columnBackward.size(); i++)
      for(UInt k=0; k<x.columns(); k++)
        if(x(i, k))
          for(const auto &col : columnBackward[i])
            y(col.first, k) = 1.;
    return y;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/
