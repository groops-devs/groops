/***********************************************/
/**
* @file gnssParametrizationKinematicPositions.cpp
*
* @brief Position estimation each epoch.
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
#include "gnss/gnssTransceiverSelector/gnssTransceiverSelector.h"
#include "gnss/gnssParametrization/gnssParametrization.h"
#include "gnss/gnssParametrization/gnssParametrizationKinematicPositions.h"

/***********************************************/

GnssParametrizationKinematicPositions::GnssParametrizationKinematicPositions(Config &config)
{
  try
  {
    readConfig(config, "name",                      name,               Config::OPTIONAL, "parameter.kinematicPositions", "used for parameter selection");
    readConfig(config, "selectReceivers",           selectReceivers,    Config::MUSTSET,  "", "");
    readConfig(config, "outputfilePositions",       fileNamePositions,  Config::OPTIONAL, "output/positions_{loopTime:%D}.{station}.dat",          "variable {station} available, estimated kinematic positions/orbit");
    readConfig(config, "outputfileCovarianceEpoch", fileNameCovariance, Config::OPTIONAL, "output/covariance3x3Epoch_{loopTime:%D}.{station}.dat", "variable {station} available, 3x3 epoch covariances");
    if(isCreateSchema(config)) return;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssParametrizationKinematicPositions::init(Gnss *gnss, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    this->gnss = gnss;
    selectedReceivers = selectReceivers->select(gnss->receivers);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssParametrizationKinematicPositions::requirements(GnssNormalEquationInfo &normalEquationInfo,
                                                         std::vector<UInt> &/*transCount*/, std::vector<UInt> &/*transCountEpoch*/,
                                                         std::vector<UInt> &/*recvCount*/,  std::vector<UInt> &recvCountEpoch)
{
  try
  {
    if(isEnabled(normalEquationInfo, name))
      for(auto recv : gnss->receivers)
        if(normalEquationInfo.estimateReceiver.at(recv->idRecv()) && selectedReceivers.at(recv->idRecv()))
          recvCountEpoch.at(recv->idRecv()) += 3;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssParametrizationKinematicPositions::initParameter(GnssNormalEquationInfo &normalEquationInfo)
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
      if(recv->useable() && selectedReceivers.at(idRecv) && normalEquationInfo.estimateReceiver.at(idRecv))
      {
        index.at(idRecv).resize(gnss->times.size());
        for(UInt idEpoch : normalEquationInfo.idEpochs)
          if(recv->useable(idEpoch))
          {
            std::vector<ParameterName> parameterNames;
            parameterNames.push_back(ParameterName(recv->name(), "position.x", "", recv->timeCorrected(idEpoch)));
            parameterNames.push_back(ParameterName(recv->name(), "position.y", "", recv->timeCorrected(idEpoch)));
            parameterNames.push_back(ParameterName(recv->name(), "position.z", "", recv->timeCorrected(idEpoch)));
            index.at(idRecv).at(idEpoch) = normalEquationInfo.parameterNamesEpochReceiver(idEpoch, idRecv, parameterNames);
            countPara += parameterNames.size();
          }
      }
    }
    if(countPara)
      logInfo<<countPara%"%9i kinematic position parameters"s<<Log::endl;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssParametrizationKinematicPositions::aprioriParameter(const GnssNormalEquationInfo &normalEquationInfo, MatrixSliceRef x0) const
{
  try
  {
    for(UInt idRecv=0; idRecv<index.size(); idRecv++)
      if(gnss->receivers.at(idRecv)->isMyRank())
        for(UInt idEpoch=0; idEpoch<index.at(idRecv).size(); idEpoch++)
          if(index.at(idRecv).at(idEpoch))
            copy(gnss->receivers.at(idRecv)->pos.at(idEpoch).vector(), x0.row(normalEquationInfo.index(index.at(idRecv).at(idEpoch)), 3));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssParametrizationKinematicPositions::designMatrix(const GnssNormalEquationInfo &/*normalEquationInfo*/, const GnssObservationEquation &eqn, GnssDesignMatrix &A) const
{
  try
  {
    if(index.at(eqn.receiver->idRecv()).size() && index.at(eqn.receiver->idRecv()).at(eqn.idEpoch))
    {
      if(eqn.receiver->isEarthFixed())
        matMult(1., eqn.A.column(GnssObservationEquation::idxPosRecv, 3),
                gnss->rotationCrf2Trf(eqn.timeRecv).matrix().trans(),
                A.column(index.at(eqn.receiver->idRecv()).at(eqn.idEpoch)));
      else
        copy(eqn.A.column(GnssObservationEquation::idxPosRecv, 3), A.column(index.at(eqn.receiver->idRecv()).at(eqn.idEpoch)));
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Double GnssParametrizationKinematicPositions::updateParameter(const GnssNormalEquationInfo &normalEquationInfo, const_MatrixSliceRef x, const_MatrixSliceRef /*Wz*/)
{
  try
  {
    Double maxChange = 0;
    Gnss::InfoParameterChange info("mm");
    for(UInt idRecv=0; idRecv<index.size(); idRecv++)
      if(gnss->receivers.at(idRecv)->isMyRank())
        for(UInt idEpoch=0; idEpoch<index.at(idRecv).size(); idEpoch++)
          if(index.at(idRecv).at(idEpoch))
          {
            const Vector3d dp(x.row(normalEquationInfo.index(index.at(idRecv).at(idEpoch)), 3));
            gnss->receivers.at(idRecv)->pos.at(idEpoch) += dp;
            if(info.update(1e3*dp.r()))
              info.info = "kinematic position ("+gnss->receivers.at(idRecv)->name()+", "+gnss->times.at(idEpoch).dateTimeStr()+")";
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

void GnssParametrizationKinematicPositions::updateCovariance(const GnssNormalEquationInfo &normalEquationInfo, const MatrixDistributed &covariance)
{
  try
  {
    // 3x3 epoch covariance matrix
    // ---------------------------
    cov.resize(index.size());
    for(UInt idRecv=0; idRecv<index.size(); idRecv++)
      if(gnss->receivers.at(idRecv)->isMyRank())
      {
        cov.at(idRecv).resize(index.at(idRecv).size());
        for(UInt idEpoch=0; idEpoch<index.at(idRecv).size(); idEpoch++)
          if(index.at(idRecv).at(idEpoch))
          {
            const UInt idBlock = normalEquationInfo.block(index.at(idRecv).at(idEpoch));
            const UInt idx     = normalEquationInfo.index(index.at(idRecv).at(idEpoch)) - covariance.blockIndex(idBlock);
            cov.at(idRecv).at(idEpoch).xx() = covariance.N(idBlock, idBlock)(idx+0, idx+0);
            cov.at(idRecv).at(idEpoch).yy() = covariance.N(idBlock, idBlock)(idx+1, idx+1);
            cov.at(idRecv).at(idEpoch).zz() = covariance.N(idBlock, idBlock)(idx+2, idx+2);
            cov.at(idRecv).at(idEpoch).xy() = covariance.N(idBlock, idBlock)(idx+0, idx+1);
            cov.at(idRecv).at(idEpoch).xz() = covariance.N(idBlock, idBlock)(idx+0, idx+2);
            cov.at(idRecv).at(idEpoch).yz() = covariance.N(idBlock, idBlock)(idx+1, idx+2);
          }
      }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssParametrizationKinematicPositions::writeResults(const GnssNormalEquationInfo &normalEquationInfo, const std::string &suffix) const
{
  try
  {
    if(!isEnabled(normalEquationInfo, name))
      return;

    VariableList fileNameVariableList;
    addVariable("station", fileNameVariableList);

    // write kinematic orbits
    // ----------------------
    if(!fileNamePositions.empty())
    {
      fileNameVariableList["station"]->setValue("****");
      logStatus<<"write kinematic position data <"<<fileNamePositions(fileNameVariableList).appendBaseName(suffix)<<">"<<Log::endl;
      for(UInt idRecv=0; idRecv<index.size(); idRecv++)
        if(gnss->receivers.at(idRecv)->isMyRank())
        {
          Arc arc;
          for(UInt idEpoch=0; idEpoch<index.at(idRecv).size(); idEpoch++)
            if(index.at(idRecv).at(idEpoch))
            {
              if(gnss->receivers.at(idRecv)->isEarthFixed())
              {
                Vector3dEpoch epoch;
                epoch.time     = gnss->receivers.at(idRecv)->timeCorrected(idEpoch);
                epoch.vector3d = gnss->receivers.at(idRecv)->pos.at(idEpoch);
                arc.push_back(epoch);
              }
              else
              {
                OrbitEpoch epoch;
                epoch.time     = gnss->receivers.at(idRecv)->timeCorrected(idEpoch);
                epoch.position = gnss->receivers.at(idRecv)->pos.at(idEpoch);
                arc.push_back(epoch);
              }
            }
          fileNameVariableList["station"]->setValue(gnss->receivers.at(idRecv)->name());
          InstrumentFile::write(fileNamePositions(fileNameVariableList).appendBaseName(suffix), arc);
        }
    }

    // write epoch covariance
    // ----------------------
    if(!fileNameCovariance.empty() && cov.size())
    {
      fileNameVariableList["station"]->setValue("****");
      logStatus<<"write epoch covariance data <"<<fileNameCovariance(fileNameVariableList).appendBaseName(suffix)<<">"<<Log::endl;
       for(UInt idRecv=0; idRecv<index.size(); idRecv++)
        if(gnss->receivers.at(idRecv)->isMyRank() && cov.at(idRecv).size())
        {
          Covariance3dArc arc;
          for(UInt idEpoch=0; idEpoch<index.at(idRecv).size(); idEpoch++)
            if(index.at(idRecv).at(idEpoch))
            {
              Covariance3dEpoch epoch;
              epoch.time       = gnss->receivers.at(idRecv)->timeCorrected(idEpoch);
              epoch.covariance = cov.at(idRecv).at(idEpoch);
              arc.push_back(epoch);
            }
          fileNameVariableList["station"]->setValue(gnss->receivers.at(idRecv)->name());
          InstrumentFile::write(fileNameCovariance(fileNameVariableList).appendBaseName(suffix), arc);
        }
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
