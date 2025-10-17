/***********************************************/
/**
* @file gnssParametrizationIonosphereVTEC.cpp
*
* @brief IonosphereVTEC.
* @see GnssParametrization
*
* @author Torsten Mayer-Guerr
* @date 2021-01-23
*
*/
/***********************************************/

#include "base/import.h"
#include "base/planets.h"
#include "config/config.h"
#include "files/fileInstrument.h"
#include "gnss/gnssParametrization/gnssParametrizationIonosphereVTEC.h"

/***********************************************/

GnssParametrizationIonosphereVTEC::GnssParametrizationIonosphereVTEC(Config &config)
{
  try
  {
    readConfig(config, "name",            name,                    Config::OPTIONAL, "parameter.VTEC", "");
    readConfig(config, "selectReceivers", selectReceivers,         Config::MUSTSET,   R"(["all"])", "");
    readConfig(config, "outputfileVTEC",  fileNameVTEC,            Config::OPTIONAL, "output/vtec_{loopTime:%D}.{station}.dat", "variable {station} available");
    readConfig(config, "mapR",            mapR,                    Config::DEFAULT,  "6371e3",     "constant of MSLM mapping function");
    readConfig(config, "mapH",            mapH,                    Config::DEFAULT,  "506.7e3",    "constant of MSLM mapping function");
    readConfig(config, "mapAlpha",        mapAlpha,                Config::DEFAULT,  "0.9782",     "constant of MSLM mapping function");
    readConfig(config, "gradient",        parametrizationGradient, Config::DEFAULT,  "",           "parametrization of north and east gradients");
    if(isCreateSchema(config)) return;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssParametrizationIonosphereVTEC::init(Gnss *gnss, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    this->gnss = gnss;
    selectedReceivers = gnss->selectReceivers(selectReceivers);
    indexVTEC.resize(gnss->receivers.size());
    VTEC.resize(gnss->receivers.size());

    indexGradient.resize(gnss->receivers.size());
    xGradient.resize(gnss->receivers.size());
    gradientX.resize(gnss->receivers.size());
    gradientY.resize(gnss->receivers.size());

    for(auto recv : gnss->receivers)
    {
      const UInt idRecv = recv->idRecv();
      if(recv->useable() && recv->isMyRank())
      {
        VTEC.at(idRecv).resize(gnss->times.size(), 0);
        xGradient.at(idRecv) = Vector(2*parametrizationGradient->parameterCount());
        gradientX.at(idRecv).resize(gnss->times.size(), 0);
        gradientY.at(idRecv).resize(gnss->times.size(), 0);
      }
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssParametrizationIonosphereVTEC::initParameter(GnssNormalEquationInfo &normalEquationInfo)
{
  try
  {
    indexVTEC.clear();
    indexVTEC.resize(gnss->receivers.size());
    indexGradient.clear();
    indexGradient.resize(gnss->receivers.size(),GnssParameterIndex());
    if(!isEnabled(normalEquationInfo, name))
      return;

    UInt countPara = 0;
    for(auto recv : gnss->receivers)
    {
      const UInt idRecv = recv->idRecv();
      if(recv->useable() && normalEquationInfo.estimateReceiver.at(idRecv) && selectedReceivers.at(idRecv))
      {
        indexVTEC.at(idRecv).resize(gnss->times.size());
        if(recv->isMyRank())
          VTEC.at(idRecv).resize(gnss->times.size(), 0);
        for(UInt idEpoch : normalEquationInfo.idEpochs)
          if(recv->useable(idEpoch))
          {
            indexVTEC.at(idRecv).at(idEpoch) = normalEquationInfo.parameterNamesEpochReceiver(idEpoch, idRecv, {ParameterName(recv->name(), "VTEC", "", gnss->times.at(idEpoch))});
            countPara++;
          }
      }
    }
    if(countPara)
      logInfo<<countPara%"%9i VTEC epoch parameters"s<<Log::endl;

    // VTEC gradient
    UInt countParaGradient = 0;
    if(parametrizationGradient->parameterCount())
      for(auto recv : gnss->receivers)
      {
        const UInt idRecv = recv->idRecv();
        if(recv->useable() && normalEquationInfo.estimateReceiver.at(idRecv) && selectedReceivers.at(idRecv))
        {
          std::vector<ParameterName> parameterNames;
          std::vector<ParameterName> name({{gnss->receivers.at(idRecv)->name(), "VTECGradient.x"}, {gnss->receivers.at(idRecv)->name(), "VTECGradient.y"}});
          parametrizationGradient->parameterName(name, parameterNames);
          indexGradient.at(idRecv) = normalEquationInfo.parameterNamesReceiver(idRecv, parameterNames);
          countParaGradient += parameterNames.size();
        }
      }
    if(countParaGradient)
      logInfo<<countParaGradient%"%9i VTEC gradient parameters"s<<Log::endl;

  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssParametrizationIonosphereVTEC::aprioriParameter(const GnssNormalEquationInfo &normalEquationInfo, MatrixSliceRef x0) const
{
  try
  {
    for(auto recv : gnss->receivers)
      if(recv->isMyRank())
      {
        const UInt idRecv = recv->idRecv();
        for(UInt idEpoch : normalEquationInfo.idEpochs)
          if(indexVTEC.at(idRecv).size() && indexVTEC.at(idRecv).at(idEpoch))
            x0(normalEquationInfo.index(indexVTEC.at(idRecv).at(idEpoch)), 0) = VTEC.at(idRecv).at(idEpoch);
        if(indexGradient.at(idRecv))
          copy(xGradient.at(idRecv), x0.row(normalEquationInfo.index(indexGradient.at(idRecv)), xGradient.at(idRecv).rows()));
      }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

// Mapping function
Double GnssParametrizationIonosphereVTEC::mapping(Angle elevation) const
{
  return 1./std::cos(std::asin(mapR/(mapR+mapH) * std::sin(mapAlpha*(PI/2-elevation))));
}

/***********************************************/

void GnssParametrizationIonosphereVTEC::designMatrix(const GnssNormalEquationInfo &/*normalEquationInfo*/, const GnssObservationEquation &eqn, GnssDesignMatrix &A) const
{
  try
  {
    const Double map = mapping(eqn.elevationRecvLocal);

    // VTEC at station per epoch
    if(indexVTEC.at(eqn.receiver->idRecv()).size() && indexVTEC.at(eqn.receiver->idRecv()).at(eqn.idEpoch))
      axpy(map, eqn.A.column(GnssObservationEquation::idxSTEC), A.column(indexVTEC.at(eqn.receiver->idRecv()).at(eqn.idEpoch)));

    // VTEC gradient at station
    if(indexGradient.at(eqn.receiver->idRecv()))
    {
      Matrix B(eqn.A.rows(), 2);
      axpy(map*std::cos(eqn.azimutRecvLocal), eqn.A.column(GnssObservationEquation::idxSTEC), B.column(0));
      axpy(map*std::sin(eqn.azimutRecvLocal), eqn.A.column(GnssObservationEquation::idxSTEC), B.column(1));
      // temporal parametrization
      std::vector<UInt>   idx;
      std::vector<Double> factor;
      parametrizationGradient->factors(std::max(eqn.timeRecv, gnss->times.at(0)), idx, factor);
      MatrixSlice Design(A.column(indexGradient.at(eqn.receiver->idRecv())));
      for(UInt i=0; i<factor.size(); i++)
        axpy(factor.at(i), B, Design.column(B.columns()*idx.at(i), B.columns()));
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Double GnssParametrizationIonosphereVTEC::updateParameter(const GnssNormalEquationInfo &normalEquationInfo, const_MatrixSliceRef x, const_MatrixSliceRef /*Wz*/)
{
  try
  {
    if(!isEnabled(normalEquationInfo, name))
      return 0;

    Double minVTEC   = 1e+99;
    Double maxVTEC   = 1e-99;
    Double meanVTEC  = 0;
    Double stdVTEC   = 0;
    UInt   countVTEC = 0;

    Gnss::InfoParameterChange infoVTEC("tec");
    Gnss::InfoParameterChange infoGradient("tec");
    for(auto recv : gnss->receivers)
      if(recv->isMyRank())
      {
        const UInt idRecv = recv->idRecv();
        Vector dxGradient;
        if(indexGradient.at(idRecv))
        {
          dxGradient = x.row(normalEquationInfo.index(indexGradient.at(idRecv)), 2*parametrizationGradient->parameterCount());
          xGradient.at(idRecv) += dxGradient;
        }

        for(UInt idEpoch=0; idEpoch<gnss->times.size(); idEpoch++)
        {
          // update VTEC
          Double dVTEC = 0;
          if(indexVTEC.at(idRecv).size() && indexVTEC.at(idRecv).at(idEpoch))
          {
            dVTEC = x(normalEquationInfo.index(indexVTEC.at(idRecv).at(idEpoch)), 0);
            VTEC.at(idRecv).at(idEpoch) += dVTEC;
            if(infoVTEC.update(dVTEC))
              infoVTEC.info = "VTEC ("+recv->name()+", "+gnss->times.at(idEpoch).dateTimeStr()+")";

            minVTEC   = std::min(VTEC.at(idRecv).at(idEpoch), minVTEC);
            maxVTEC   = std::max(VTEC.at(idRecv).at(idEpoch), maxVTEC);
            meanVTEC += VTEC.at(idRecv).at(idEpoch);
            stdVTEC  += VTEC.at(idRecv).at(idEpoch)*VTEC.at(idRecv).at(idEpoch);
            countVTEC++;
          }

          // update gradient
          Double dgx=0, dgy=0;
          if(indexGradient.at(idRecv))
          {
            std::vector<UInt>   index;
            std::vector<Double> factor;
            parametrizationGradient->factors(std::max(recv->timeCorrected(idEpoch), gnss->times.at(0)), index, factor);
            for(UInt k=0; k<factor.size(); k++)
            {
              dgx += factor.at(k) * dxGradient(2*index.at(k)+0);
              dgy += factor.at(k) * dxGradient(2*index.at(k)+1);
            }
            gradientX.at(idRecv).at(idEpoch) += dgx;
            gradientY.at(idRecv).at(idEpoch) += dgy;
            if(infoGradient.update(std::sqrt(dgx*dgx+dgy*dgy)))
              infoGradient.info = "VTEC gradient ("+recv->name()+", "+gnss->times.at(idEpoch).dateTimeStr()+")";
          }

          // update STEC
          if(dVTEC || dgx || dgy)
            for(UInt idTrans=0; idTrans<recv->idTransmitterSize(idEpoch); idTrans++)
              if(recv->observation(idTrans, idEpoch))
              {
                GnssObservationEquation eqn(*recv->observation(idTrans, idEpoch), *recv, *gnss->transmitters.at(idTrans), gnss->funcRotationCrf2Trf,
                                            nullptr/*reduceModels*/, idEpoch, FALSE/*homogenize*/, {}/*types*/);
                recv->observation(idTrans, idEpoch)->STEC += mapping(eqn.elevationRecvLocal) * (dVTEC + std::cos(eqn.azimutRecvLocal)*dgx + std::sin(eqn.azimutRecvLocal)*dgy);
              }
        } // for(idEpoch)
      } // for(recv)

    // VTEC statistics
    Parallel::reduceMin(minVTEC,   0, normalEquationInfo.comm);
    Parallel::reduceMax(maxVTEC,   0, normalEquationInfo.comm);
    Parallel::reduceSum(meanVTEC,  0, normalEquationInfo.comm);
    Parallel::reduceSum(stdVTEC,   0, normalEquationInfo.comm);
    Parallel::reduceSum(countVTEC, 0, normalEquationInfo.comm);
    stdVTEC   = std::sqrt((stdVTEC-meanVTEC*meanVTEC/countVTEC)/(countVTEC-1));
    meanVTEC /= countVTEC;
    std::string infoStr = " (total: "+meanVTEC%"%.2f +- "s+stdVTEC%"%.2f ["s+minVTEC%"%.2f -- "s+maxVTEC%"%.2f])"s;
    Parallel::broadCast(infoStr, 0, normalEquationInfo.comm);
    infoVTEC.info += infoStr;

    Double maxChange = 0;
    infoVTEC.synchronizeAndPrint    (normalEquationInfo.comm, 0, maxChange);
    infoGradient.synchronizeAndPrint(normalEquationInfo.comm, 0, maxChange);
    return maxChange;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssParametrizationIonosphereVTEC::writeResults(const GnssNormalEquationInfo &normalEquationInfo, const std::string &suffix) const
{
  try
  {
    if(!isEnabled(normalEquationInfo, name))
      return;

    if(!fileNameVTEC.empty() && VTEC.size())
    {
      VariableList fileNameVariableList;
      fileNameVariableList.setVariable("station", "****");
      logStatus<<"write VTEC to files <"<<fileNameVTEC(fileNameVariableList).appendBaseName(suffix)<<">"<<Log::endl;

      for(auto &recv : gnss->receivers)
        if(normalEquationInfo.estimateReceiver.at(recv->idRecv()) && recv->isMyRank())
        {
          MiscValuesArc arc;
          for(UInt idEpoch : normalEquationInfo.idEpochs)
          {
            if(indexVTEC.at(recv->idRecv()).size() && indexVTEC.at(recv->idRecv()).at(idEpoch))
            {
              MiscValuesEpoch epoch(3);
              epoch.time      = gnss->times.at(idEpoch);
              epoch.values(0) = VTEC.at(recv->idRecv()).at(idEpoch);
              epoch.values(1) = gradientX.at(recv->idRecv()).at(idEpoch);
              epoch.values(2) = gradientY.at(recv->idRecv()).at(idEpoch);
              arc.push_back(epoch);
            }
          }
          fileNameVariableList.setVariable("station", recv->name());
          if(arc.size())
            InstrumentFile::write(fileNameVTEC(fileNameVariableList).appendBaseName(suffix), arc);
        }
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
