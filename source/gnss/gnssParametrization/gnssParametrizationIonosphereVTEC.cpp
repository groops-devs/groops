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
#include "classes/magnetosphere/magnetosphere.h"
#include "gnss/gnssParametrization/gnssParametrizationIonosphereVTEC.h"

/***********************************************/

GnssParametrizationIonosphereVTEC::GnssParametrizationIonosphereVTEC(Config &config)
{
  try
  {
    readConfig(config, "name",              name,              Config::OPTIONAL, "parameter.VTEC", "");
    readConfig(config, "selectReceivers",   selectReceivers,   Config::MUSTSET,   R"(["all"])", "");
    readConfig(config, "outputfileVTEC",    fileNameVTEC,      Config::OPTIONAL, "output/vtec_{loopTime:%D}.{station}.dat", "variable {station} available");
    readConfig(config, "mapR",              mapR,              Config::DEFAULT,  "6371e3",     "constant of MSLM mapping function");
    readConfig(config, "mapH",              mapH,              Config::DEFAULT,  "506.7e3",    "constant of MSLM mapping function");
    readConfig(config, "mapAlpha",          mapAlpha,          Config::DEFAULT,  "0.9782",     "constant of MSLM mapping function");
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
    index.resize(gnss->receivers.size());
    VTEC.resize(gnss->receivers.size());
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
    index.clear();
    index.resize(gnss->receivers.size());
    if(!isEnabled(normalEquationInfo, name))
      return;

    UInt countPara = 0;
    for(auto recv : gnss->receivers)
    {
      const UInt idRecv = recv->idRecv();
      if(recv->useable() && normalEquationInfo.estimateReceiver.at(idRecv) && selectedReceivers.at(idRecv))
      {
        index.at(idRecv).resize(gnss->times.size());
        if(recv->isMyRank())
          VTEC.at(idRecv).resize(gnss->times.size(), 0);
        for(UInt idEpoch : normalEquationInfo.idEpochs)
          if(recv->useable(idEpoch))
          {
            index.at(idRecv).at(idEpoch) = normalEquationInfo.parameterNamesEpochReceiver(idEpoch, idRecv, {ParameterName(recv->name(), "VTEC", "", gnss->times.at(idEpoch))});
            countPara++;
          }
      }
    }
    if(countPara)
      logInfo<<countPara%"%9i VTEC epoch parameters"s<<Log::endl;
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
          if(index.at(idRecv).size() && index.at(idRecv).at(idEpoch))
            x0(normalEquationInfo.index(index.at(idRecv).at(idEpoch)), 0) = VTEC.at(idRecv).at(idEpoch);
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
    if(!index.at(eqn.receiver->idRecv()).size() || !index.at(eqn.receiver->idRecv()).at(eqn.idEpoch))
      return;

    // VTEC at station per epoch
    axpy(mapping(eqn.elevationRecvLocal), eqn.A.column(GnssObservationEquation::idxSTEC), A.column(index.at(eqn.receiver->idRecv()).at(eqn.idEpoch)));
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

    Gnss::InfoParameterChange info("tec");
    for(auto recv : gnss->receivers)
      if(recv->isMyRank())
      {
        const UInt idRecv = recv->idRecv();
        for(UInt idEpoch : normalEquationInfo.idEpochs)
          if(index.at(idRecv).size() && index.at(idRecv).at(idEpoch))
          {
            // update VTEC
            const Double dVTEC = x(normalEquationInfo.index(index.at(idRecv).at(idEpoch)), 0);
            VTEC.at(idRecv).at(idEpoch) += dVTEC;
            if(info.update(dVTEC))
              info.info = "VTEC ("+recv->name()+", "+gnss->times.at(idEpoch).dateTimeStr()+")";

            minVTEC   = std::min(VTEC.at(idRecv).at(idEpoch), minVTEC);
            maxVTEC   = std::max(VTEC.at(idRecv).at(idEpoch), maxVTEC);
            meanVTEC += VTEC.at(idRecv).at(idEpoch);
            stdVTEC  += VTEC.at(idRecv).at(idEpoch)*VTEC.at(idRecv).at(idEpoch);
            countVTEC++;

            // update STEC
            for(UInt idTrans=0; idTrans<recv->idTransmitterSize(idEpoch); idTrans++)
              if(recv->observation(idTrans, idEpoch))
              {
                GnssObservationEquation eqn(*recv->observation(idTrans, idEpoch), *recv, *gnss->transmitters.at(idTrans), gnss->funcRotationCrf2Trf,
                                            nullptr/*reduceModels*/, idEpoch, FALSE/*decorrelate*/, {}/*types*/);
                recv->observation(idTrans, idEpoch)->STEC += mapping(eqn.elevationRecvLocal) * dVTEC;
              }
          }
      }

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
    info.info += infoStr;

    Double maxChange = 0;
    info.synchronizeAndPrint(normalEquationInfo.comm, 0, maxChange);
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
          MiscValueArc arc;
          for(UInt idEpoch : normalEquationInfo.idEpochs)
            if(index.at(recv->idRecv()).size() && index.at(recv->idRecv()).at(idEpoch))
            {
              MiscValueEpoch epoch;
              epoch.time  = gnss->times.at(idEpoch);
              epoch.value = VTEC.at(recv->idRecv()).at(idEpoch);
              arc.push_back(epoch);
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
