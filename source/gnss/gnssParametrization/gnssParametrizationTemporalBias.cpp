/***********************************************/
/**
* @file gnssParametrizationTemporalBias.cpp
*
* @brief Temporal changing signal bias.
* @see GnssParametrization
*
* @author Torsten Mayer-Guerr
* @author Sebastian Strasser
* @date 2021-01-23
*
*/
/***********************************************/

#include "base/import.h"
#include "base/polynomial.h"
#include "config/config.h"
#include "files/fileInstrument.h"
#include "gnss/gnss.h"
#include "classes/parametrizationTemporal/parametrizationTemporal.h"
#include "gnss/gnssTransceiverSelector/gnssTransceiverSelector.h"
#include "gnss/gnssParametrization/gnssParametrizationTemporalBias.h"

/***********************************************/

GnssParametrizationTemporalBias::GnssParametrizationTemporalBias(Config &config)
{
  try
  {
    readConfig(config, "name",                     name,               Config::OPTIONAL, "parameter.temporalBias", "used for parameter selection");
    readConfig(config, "selectTransmitters",       selectTransmitters, Config::MUSTSET,  "",       "");
    readConfig(config, "outputfileBiasTimeSeries", fileNameOut,        Config::OPTIONAL, "",       "variable {prn} available");
    readConfig(config, "inputfileBiasTimeSeries",  fileNameIn,         Config::OPTIONAL, "",       "variable {prn} available");
    readConfig(config, "type",                     type,               Config::MUSTSET,  "L5*G",   "");
    readConfig(config, "parametrizationTemporal",  temporal,           Config::DEFAULT,  "",       "");
    readConfig(config, "nameConstraint",           nameConstraint,     Config::OPTIONAL, "constraint.temporalBias", "used for parameter selection");
    readConfig(config, "sigmaZeroMeanConstraint",  sigmaZeroMean,      Config::DEFAULT,  "0.0001", "(0 = unconstrained) sigma [m] for temporal zero-mean constraint");
    if(isCreateSchema(config)) return;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

GnssParametrizationTemporalBias::~GnssParametrizationTemporalBias()
{
  for(Parameter *para : parameters)
    delete para;
}

/***********************************************/

void GnssParametrizationTemporalBias::init(Gnss *gnss, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    this->gnss = gnss;
    auto selectedTransmitters = selectTransmitters->select(gnss->transmitters);

    VariableList fileNameVariableList;
    addVariable("prn", "***", fileNameVariableList);
    Bool foundAnyFile = FALSE;
    parameters.resize(gnss->transmitters.size(), nullptr);
    for(UInt idTrans=0; idTrans<gnss->transmitters.size(); idTrans++)
      if(selectedTransmitters.at(idTrans) && gnss->transmitters.at(idTrans)->useable())
      {
        auto para = new Parameter();
        parameters.at(idTrans) = para;
        para->trans = gnss->transmitters.at(idTrans);

        para->bias = Vector(gnss->times.size());
        if(!fileNameIn.empty())
        {
          try
          {
            fileNameVariableList["prn"]->setValue(para->trans->name());
            MiscValueArc arc = InstrumentFile::read(fileNameIn(fileNameVariableList));
            Polynomial polynomial(3);
            para->bias = polynomial.interpolate(gnss->times, arc.times(), arc.matrix().column(1));
            foundAnyFile = TRUE;
          }
          catch(std::exception &/*e*/)
          {
          }
        }

        if(temporal->parameterCount())
        {
          Matrix A(gnss->times.size(), temporal->parameterCount());
          for(UInt idEpoch=0; idEpoch<gnss->times.size(); idEpoch++)
            copy(temporal->factors(gnss->times.at(idEpoch)).trans(), A.row(idEpoch));
          para->x = leastSquares(A, Vector(para->bias));
        }
      }

    if(!fileNameIn.empty() && !foundAnyFile)
      logWarningOnce<<"time-variable signal bias files for "<<type.str()<<" not found for any transmitter."<<Log::endl;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssParametrizationTemporalBias::observationCorrections(GnssObservationEquation &eqn) const
{
  try
  {
    auto para = parameters.at(eqn.transmitter->idTrans());
    if(para)
      for(UInt idType=0; idType<eqn.typesTransmitted.size(); idType++)
        if(eqn.typesTransmitted.at(idType) == type)
          axpy(-para->bias(eqn.idEpoch), eqn.A.column(GnssObservationEquation::idxUnit + eqn.types.size() + idType), eqn.l);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssParametrizationTemporalBias::initParameter(GnssNormalEquationInfo &normalEquationInfo)
{
  try
  {
    for(auto para : parameters)
      if(para)
        para->index = GnssParameterIndex();
    applyConstraint = FALSE;
    if(!isEnabled(normalEquationInfo, name) || normalEquationInfo.isEachReceiverSeparately)
      return;

    UInt countPara = 0;
    std::vector<GnssType> types = GnssType::replaceCompositeSignals(gnss->types());
    std::vector<ParameterName> parameterNames;
    temporal->parameterName({ParameterName("", "signalBias."+type.str())}, parameterNames);
    for(auto para : parameters)
      if(para && para->trans->useable() && (para->trans->PRN() == type) && (GnssType::index(types, para->trans->PRN()+type) != NULLINDEX))
      {
        for(auto &name : parameterNames)
          name.object = para->trans->name();
        para->index = normalEquationInfo.parameterNamesTransmitter(para->trans->idTrans(), parameterNames);
        countPara += parameterNames.size();
      }
    if(countPara)
      logInfo<<countPara%"%9i temporal bias ("s<<type.str()<<") parameters"<<Log::endl;

    applyConstraint = isEnabled(normalEquationInfo, nameConstraint)  && sigmaZeroMean && countPara;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssParametrizationTemporalBias::aprioriParameter(const GnssNormalEquationInfo &normalEquationInfo, MatrixSliceRef x0) const
{
  try
  {
    if(Parallel::isMaster(normalEquationInfo.comm))
      for(auto para : parameters)
        if(para && para->index)
          copy(para->x, x0.row(normalEquationInfo.index(para->index), para->x.rows()));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssParametrizationTemporalBias::designMatrix(const GnssNormalEquationInfo &/*normalEquationInfo*/, const GnssObservationEquation &eqn, GnssDesignMatrix &A) const
{
  try
  {
    auto para = parameters.at(eqn.transmitter->idTrans());
    if(para && para->index)
    {
      MatrixSlice Design(A.column(para->index));
      for(UInt idType=0; idType<eqn.typesTransmitted.size(); idType++)
        if(eqn.typesTransmitted.at(idType) == type)
          matMult(1., eqn.A.column(GnssObservationEquation::idxUnit + eqn.types.size() + idType),
                  temporal->factors(std::max(eqn.timeTrans, gnss->times.at(0))).trans(), Design);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssParametrizationTemporalBias::constraints(const GnssNormalEquationInfo &normalEquationInfo, MatrixDistributed &normals, std::vector<Matrix> &n, Double &lPl, UInt &obsCount) const
{
  try
  {
    if(!Parallel::isMaster(normalEquationInfo.comm) || !applyConstraint)
      return;

    logStatus<<"apply temporal zero mean of temporal signal bias of "<<type.str()<<Log::endl;
    Vector mean(temporal->parameterCount());
    for(const Time &time : gnss->times)
      mean += temporal->factors(time);

    for(auto para : parameters)
      if(para && para->index)
      {
        GnssDesignMatrix A(normalEquationInfo, Vector(1));
        A.init(-1./sigmaZeroMean/gnss->times.size() * mean.trans() * para->x); // constrain towards zero (0-x0)
        axpy(1./sigmaZeroMean/gnss->times.size(), mean.trans(), A.column(para->index));
        A.accumulateNormals(normals, n, lPl, obsCount);
      }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Double GnssParametrizationTemporalBias::updateParameter(const GnssNormalEquationInfo &normalEquationInfo, const_MatrixSliceRef x, const_MatrixSliceRef /*Wz*/)
{
  try
  {
    Double maxChange = 0;
    Gnss::InfoParameterChange info("mm");
    for(auto para : parameters)
      if(para && para->index)
      {
        const Vector dx = x.row(normalEquationInfo.index(para->index), para->x.rows());
        para->x += dx;
        for(UInt idEpoch=0; idEpoch<gnss->times.size(); idEpoch++)
        {
          const Double dBias = inner(temporal->factors(gnss->times.at(idEpoch)), dx);
          para->bias(idEpoch) += dBias;
          if(info.update(1e3*dBias))
            info.info = "temporal signal bias ("+type.str()+")";
        }
      }
    info.synchronizeAndPrint(normalEquationInfo.comm, 1e-3, maxChange);
    return maxChange;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssParametrizationTemporalBias::writeResults(const GnssNormalEquationInfo &normalEquationInfo, const std::string &suffix) const
{
  try
  {
    if(!isEnabled(normalEquationInfo, name) || !Parallel::isMaster(normalEquationInfo.comm))
      return;

    if(!fileNameOut.empty() && std::any_of(parameters.begin(), parameters.end(), [](const Parameter *p){return p && p->index;}))
    {
      VariableList fileNameVariableList;
      addVariable("prn", "***", fileNameVariableList);
      logStatus<<"write transmitter time variable bias to files <"<<fileNameOut(fileNameVariableList).appendBaseName(suffix)<<">"<<Log::endl;
      for(auto para : parameters)
        if(para && para->index)
        {
          MiscValueArc arc;
          for(UInt idEpoch : normalEquationInfo.idEpochs)
            if(para->trans->useable(idEpoch))
            {
              MiscValueEpoch epoch;
              epoch.time  = gnss->times.at(idEpoch);
              epoch.value = para->bias(idEpoch);
              arc.push_back(epoch);
            }
          fileNameVariableList["prn"]->setValue(para->trans->name());
          InstrumentFile::write(fileNameOut(fileNameVariableList).appendBaseName(suffix), arc);
        }
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
