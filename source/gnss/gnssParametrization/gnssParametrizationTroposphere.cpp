/***********************************************/
/**
* @file gnssParametrizationTroposphere.cpp
*
* @brief Troposphere.
* @see GnssParametrization
*
* @author Torsten Mayer-Guerr
* @author Sebastian Strasser
* @date 2021-01-23
*
*/
/***********************************************/

#include "base/import.h"
#include "config/config.h"
#include "files/fileMatrix.h"
#include "classes/troposphere/troposphere.h"
#include "classes/parametrizationTemporal/parametrizationTemporal.h"
#include "gnss/gnssTransceiverSelector/gnssTransceiverSelector.h"
#include "gnss/gnssParametrization/gnssParametrizationTroposphere.h"

/***********************************************/

GnssParametrizationTroposphere::GnssParametrizationTroposphere(Config &config)
{
  try
  {
    readConfig(config, "name",                          name,                    Config::OPTIONAL, "parameter.troposphere", "used for parameter selection");
    readConfig(config, "selectReceivers",               selectReceivers,         Config::MUSTSET,  "",  "");
    readConfig(config, "outputfileTroposphere",         fileNameTropo,           Config::OPTIONAL, "output/troposphere_{loopTime:%D}.{station}.txt", "columns: MJD, ZHD, ZWD, dry north gradient, wet north gradient, dry east gradient, wet east gradient");
    readConfig(config, "troposphere",                   troposphere,             Config::MUSTSET,  "",  "a priori troposphere model");
    readConfig(config, "troposphereWetEstimation",      parametrizationWet,      Config::DEFAULT,  "",  "[m] parametrization of zenith wet delays");
    readConfig(config, "troposphereGradientEstimation", parametrizationGradient, Config::DEFAULT,  "",  "[degree] parametrization of north and east gradients");
    if(isCreateSchema(config)) return;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

GnssParametrizationTroposphere::~GnssParametrizationTroposphere()
{
  for(Parameter *para : parameters)
    delete para;
}

/***********************************************/

void GnssParametrizationTroposphere::init(Gnss *gnss, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    this->gnss = gnss;
    auto selectedReceivers = selectReceivers->select(gnss->receivers);
    parameters.resize(gnss->receivers.size(), nullptr);
    UInt idTropo = 0;
    std::vector<Vector3d> positions;
    for(UInt idRecv=0; idRecv<gnss->receivers.size(); idRecv++)
      if(selectedReceivers.at(idRecv) && gnss->receivers.at(idRecv)->useable())
      {
        auto para = new Parameter();
        parameters.at(idRecv) = para;
        para->idRecv = idRecv;
        if(gnss->receivers.at(idRecv)->isMyRank())
        {
          positions.push_back(gnss->receivers.at(idRecv)->position(0));
          para->idTropo   = idTropo++;
          para->xWet      = Vector(parametrizationWet->parameterCount());
          para->xGradient = Vector(2*parametrizationGradient->parameterCount());
          para->zenitDelayWet.resize(gnss->times.size(), 0);
          para->gradientX.resize(gnss->times.size(), 0);
          para->gradientY.resize(gnss->times.size(), 0);
        }
      }

    troposphere->init(positions);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssParametrizationTroposphere::initParameter(GnssNormalEquationInfo &normalEquationInfo)
{
  try
  {
    for(auto para : parameters)
      if(para)
        para->indexWet = para->indexGradient = GnssParameterIndex();
    if(!isEnabled(normalEquationInfo, name))
      return;

    // wet troposphere
    UInt countParaWet = 0;
    if(parametrizationWet->parameterCount())
      for(auto para : parameters)
        if(para && normalEquationInfo.estimateReceiver.at(para->idRecv))
        {
          std::vector<ParameterName> parameterNames;
          parametrizationWet->parameterName({ParameterName(gnss->receivers.at(para->idRecv)->name(), "troposphereWet")}, parameterNames);
          para->indexWet = normalEquationInfo.parameterNamesReceiver(para->idRecv, parameterNames);
          countParaWet += parameterNames.size();
        }
    if(countParaWet)
      logInfo<<countParaWet%"%9i troposphere wet parameters"s<<Log::endl;

    // troposphere gradient
    UInt countParaGradient = 0;
    if(parametrizationGradient->parameterCount())
      for(auto para : parameters)
        if(para && normalEquationInfo.estimateReceiver.at(para->idRecv))
        {
          std::vector<ParameterName> parameterNames;
          std::vector<ParameterName> name({{gnss->receivers.at(para->idRecv)->name(), "troposphereGradient.x"}, {gnss->receivers.at(para->idRecv)->name(), "troposphereGradient.y"}});
          parametrizationGradient->parameterName(name, parameterNames);
          para->indexGradient = normalEquationInfo.parameterNamesReceiver(para->idRecv, parameterNames);
          countParaGradient += parameterNames.size();
        }
    if(countParaGradient)
      logInfo<<countParaGradient%"%9i troposphere gradient parameters"s<<Log::endl;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssParametrizationTroposphere::observationCorrections(GnssObservationEquation &eqn) const
{
  try
  {
    auto para = parameters.at(eqn.receiver->idRecv());
    if(!para)
      return;

    const Time t = std::max(eqn.timeRecv, gnss->times.at(0));
    // apriori value
    Double delay = troposphere->slantDelay(t, para->idTropo, eqn.azimutRecvLocal, eqn.elevationRecvLocal);
    // estimated wet effect
    delay += troposphere->mappingFunctionWet(t, para->idTropo, eqn.azimutRecvLocal, eqn.elevationRecvLocal) * para->zenitDelayWet.at(eqn.idEpoch);
    // estimated gradient
    Double dx, dy;
    troposphere->mappingFunctionGradient(t, para->idTropo, eqn.azimutRecvLocal, eqn.elevationRecvLocal, dx, dy);
    delay += dx * para->gradientX.at(eqn.idEpoch) + dy * para->gradientY.at(eqn.idEpoch);

    for(UInt i=0; i<eqn.types.size(); i++)
      if((eqn.types.at(i) == GnssType::RANGE) || (eqn.types.at(i) == GnssType::PHASE))
        eqn.l(i) -= delay;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssParametrizationTroposphere::aprioriParameter(const GnssNormalEquationInfo &normalEquationInfo, MatrixSliceRef x0) const
{
  try
  {
    // update wet troposphere
    for(auto para : parameters)
      if(para && para->indexWet && gnss->receivers.at(para->idRecv)->isMyRank())
        copy(para->xWet, x0.row(normalEquationInfo.index(para->indexWet), para->xWet.rows()));

    // update troposphere gradient
    for(auto para : parameters)
      if(para && para->indexGradient && gnss->receivers.at(para->idRecv)->isMyRank())
        copy(para->xGradient, x0.row(normalEquationInfo.index(para->indexGradient), para->xGradient.rows()));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssParametrizationTroposphere::designMatrix(const GnssNormalEquationInfo &/*normalEquationInfo*/, const GnssObservationEquation &eqn, GnssDesignMatrix &A) const
{
  try
  {
    auto para = parameters.at(eqn.receiver->idRecv());
    if(!para)
      return;

    // temporal parametrization
    auto designMatrixTemporal = [&](ParametrizationTemporalPtr parametrization, const_MatrixSliceRef B, const GnssParameterIndex &index)
    {
      std::vector<UInt>   idx;
      std::vector<Double> factor;
      parametrization->factors(std::max(eqn.timeRecv, gnss->times.at(0)), idx, factor);
      MatrixSlice Design(A.column(index));
      for(UInt i=0; i<factor.size(); i++)
        axpy(factor.at(i), B, Design.column(B.columns()*idx.at(i), B.columns()));
    };

    // troposphere wet
    if(para->indexWet)
    {
      const Double mappingFunctionWet = troposphere->mappingFunctionWet(std::max(eqn.timeRecv, gnss->times.at(0)), para->idTropo, eqn.azimutRecvLocal, eqn.elevationRecvLocal);
      const Matrix B = mappingFunctionWet * eqn.A.column(GnssObservationEquation::idxRange,1);
      designMatrixTemporal(parametrizationWet, B, para->indexWet);
    }

    // troposphere gradient
    if(para->indexGradient)
    {
      Double dx, dy;
      troposphere->mappingFunctionGradient(std::max(eqn.timeRecv, gnss->times.at(0)), para->idTropo, eqn.azimutRecvLocal, eqn.elevationRecvLocal, dx, dy);
      Matrix B(eqn.A.rows(), 2);
      axpy(dx, eqn.A.column(GnssObservationEquation::idxRange,1), B.column(0));
      axpy(dy, eqn.A.column(GnssObservationEquation::idxRange,1), B.column(1));
      designMatrixTemporal(parametrizationGradient, B, para->indexGradient);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Double GnssParametrizationTroposphere::updateParameter(const GnssNormalEquationInfo &normalEquationInfo, const_MatrixSliceRef x, const_MatrixSliceRef /*Wz*/)
{
  try
  {
    // update wet troposphere
    Double maxChange = 0;
    Gnss::InfoParameterChange infoWet("mm");
    for(auto para : parameters)
      if(para && para->indexWet && gnss->receivers.at(para->idRecv)->isMyRank())
      {
        auto recv = gnss->receivers.at(para->idRecv);
        para->xWet += x.row(normalEquationInfo.index(para->indexWet), parametrizationWet->parameterCount());
        std::vector<UInt>   index;
        std::vector<Double> factor;
        for(UInt idEpoch=0; idEpoch<gnss->times.size(); idEpoch++)
        {
          parametrizationWet->factors(std::max(recv->timeCorrected(idEpoch), gnss->times.at(0)), index, factor);
          Double z = 0;
          for(UInt k=0; k<factor.size(); k++)
            z += factor.at(k) * para->xWet(index.at(k));
          const Double zOld = para->zenitDelayWet.at(idEpoch);
          para->zenitDelayWet.at(idEpoch) = z;
          if(infoWet.update(1e3*(z-zOld)))
            infoWet.info = "troposphere wet ("+recv->name()+", "+gnss->times.at(idEpoch).dateTimeStr()+")";
        }
      }
    infoWet.synchronizeAndPrint(normalEquationInfo.comm, 1e-3, maxChange);

    // update troposphere gradient
    Gnss::InfoParameterChange infoGradient("mm");
    for(auto para : parameters)
      if(para && para->indexGradient && gnss->receivers.at(para->idRecv)->isMyRank())
      {
        auto recv = gnss->receivers.at(para->idRecv);
        para->xGradient += x.row(normalEquationInfo.index(para->indexGradient), 2*parametrizationGradient->parameterCount());
        std::vector<UInt>   index;
        std::vector<Double> factor;
        for(UInt idEpoch=0; idEpoch<gnss->times.size(); idEpoch++)
        {
          parametrizationGradient->factors(std::max(recv->timeCorrected(idEpoch), gnss->times.at(0)), index, factor);
          Double gx=0, gy=0;
          for(UInt k=0; k<factor.size(); k++)
          {
            gx += factor.at(k) * para->xGradient(2*index.at(k)+0);
            gy += factor.at(k) * para->xGradient(2*index.at(k)+1);
          }
          const Double dx = gx - para->gradientX.at(idEpoch);
          const Double dy = gy - para->gradientY.at(idEpoch);
          para->gradientX.at(idEpoch) = gx;
          para->gradientY.at(idEpoch) = gy;
          const Double dxdy = std::sqrt(dx*dx+dy*dy);
          if(infoGradient.update(1e3*dxdy))
            infoGradient.info = "troposphere gradient ("+recv->name()+", "+gnss->times.at(idEpoch).dateTimeStr()+")";
        }
      }
    infoGradient.synchronizeAndPrint(normalEquationInfo.comm, 1e-3, maxChange);

    return maxChange;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssParametrizationTroposphere::writeResults(const GnssNormalEquationInfo &normalEquationInfo, const std::string &suffix) const
{
  try
  {
    if(!isEnabled(normalEquationInfo, name) || fileNameTropo.empty())
      return;

    VariableList fileNameVariableList;
    addVariable("station", "****", fileNameVariableList);
    fileNameVariableList["station"]->setValue("****");
    for(auto para : parameters)
      if(para && normalEquationInfo.estimateReceiver.at(para->idRecv) && gnss->receivers.at(para->idRecv)->isMyRank())
      {
        Matrix A(normalEquationInfo.idEpochs.size(), 12);
        for(UInt i=0; i<normalEquationInfo.idEpochs.size(); i++)
        {
          const UInt idEpoch = normalEquationInfo.idEpochs.at(i);
          Double zenithWetDelay, zenithDryDelay, gradientWetNorth, gradientDryNorth, gradientWetEast, gradientDryEast, aDry, aWet;
          troposphere->getAprioriValues(gnss->times.at(idEpoch), para->idTropo, zenithDryDelay, zenithWetDelay,
                                        gradientDryNorth, gradientWetNorth, gradientDryEast, gradientWetEast, aDry, aWet);
          A(i,  0) = gnss->times.at(idEpoch).mjd();
          A(i,  1) = zenithDryDelay;                                     // tropospheric zenith dry delay [m] (only from model)
          A(i,  2) = zenithWetDelay   + para->zenitDelayWet.at(idEpoch); // tropospheric zenith wet delay [m] (model + delta estimate)
          A(i,  3) = gradientDryNorth + para->gradientX.at(idEpoch);     // tropospheric dry gradient - north direction [m] (model + delta estimate, due to same mapping function)
          A(i,  4) = gradientWetNorth;                                   // tropospheric wet gradient - north direction [m] (only from model)
          A(i,  5) = gradientDryEast  + para->gradientY.at(idEpoch);     // tropospheric dry gradient - east component [m] (model + delta estimate, due to same mapping function)
          A(i,  6) = gradientWetEast;                                    // tropospheric wet gradient - east component [m] (only from model)
          A(i,  7) = para->zenitDelayWet.at(idEpoch);                    // tropospheric zenith wet delay [m] (delta estimate)
          A(i,  8) = para->gradientX.at(idEpoch);                        // tropospheric gradient - north [m] (delta estimate)
          A(i,  9) = para->gradientY.at(idEpoch);                        // tropospheric gradient - east  [m] (delta estimate)
          A(i, 10) = aDry;                                               // dry mapping function coefficient a []
          A(i, 11) = aWet;                                               // wet mapping function coefficient a []
        }

        fileNameVariableList["station"]->setValue(gnss->receivers.at(para->idRecv)->name());
        writeFileMatrix(fileNameTropo(fileNameVariableList).appendBaseName(suffix), A);
      }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
