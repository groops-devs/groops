/***********************************************/
/**
* @file gnssParametrizationIonosphereMap.cpp
*
* @brief IonosphereMap.
* @see GnssParametrization
*
* @author Torsten Mayer-Guerr
* @date 2021-01-23
*
*/
/***********************************************/

#include "base/import.h"
#include "base/planets.h"
#include "base/sphericalHarmonics.h"
#include "config/config.h"
// #include "files/fileGnssIonosphereMaps.h"
#include "classes/magnetosphere/magnetosphere.h"
#include "gnss/gnssParametrization/gnssParametrizationIonosphereMap.h"

/***********************************************/

GnssParametrizationIonosphereMap::GnssParametrizationIonosphereMap(Config &config)
{
  try
  {
    readConfig(config, "name",            name,            Config::OPTIONAL, "parameter.mapVTEC", "");
    readConfig(config, "selectReceivers", selectReceivers, Config::MUSTSET,  R"(["all"])", "");
//     readConfig(config, "outputfileMap",   fileNameOut,     Config::OPTIONAL, "",        "");
//     readConfig(config, "inputfileMap",    fileNameIn,      Config::OPTIONAL, "",        "");
    readConfig(config, "maxDegree",       maxDegree,       Config::MUSTSET,  "15",      "spherical harmonics");
    readConfig(config, "temporal",        temporal,        Config::DEFAULT,  "",        "temporal evolution of TEC values");
    readConfig(config, "mapR",            mapR,            Config::DEFAULT,  "6371e3",  "[m] constant of MSLM mapping function");
    readConfig(config, "mapH",            mapH,            Config::DEFAULT,  "506.7e3", "[m] constant of MSLM mapping function");
    readConfig(config, "mapAlpha",        mapAlpha,        Config::DEFAULT,  "0.9782",  "constant of MSLM mapping function");
    readConfig(config, "magnetosphere",   magnetosphere,   Config::MUSTSET,  "",        "");
    if(isCreateSchema(config)) return;

  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssParametrizationIonosphereMap::init(Gnss *gnss, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    this->gnss = gnss;
    selectedReceivers = gnss->selectReceivers(selectReceivers);
    x.resize(temporal->parameterCount(), Vector((maxDegree+1)*(maxDegree+1)));

    // apriori model
    // -------------
//     if(!fileNameIn.empty())
//     {
//       GnssIonosphereMapsFile file(fileNameIn, maxDegree);
//
//       // aprioriParameter by least squares adjustment
//       Matrix A(gnss->times.size(), x.size());
//       if(x.size())
//       {
//         for(UInt idEpoch=0; idEpoch<gnss->times.size(); idEpoch++)
//           copy(temporal->factors(gnss->times.at(idEpoch)).trans(), A.row(idEpoch));
//         // A' := (A'A)^-1A'
//         const Vector tau = QR_decomposition(A);
//         const Matrix W = A.row(0, A.columns());
//         generateQ(A, tau);
//         triangularSolve(1., W, A.trans());
//       }
//
//       for(UInt idEpoch=0; idEpoch<gnss->times.size(); idEpoch++)
//       {
//         const Vector xt = file.sphericalHarmonics(gnss->times.at(idEpoch)).x();
//         for(UInt i=0; i<A.columns(); i++)
//           axpy(A(idEpoch, i), xt, x.at(i));
//
//         // update STEC in observations
//         for(auto recv : gnss->receivers)
//           if(recv->isMyRank() && recv->useable(idEpoch) && selectedReceivers.at(recv->idRecv()))
//             for(UInt idTrans=0; idTrans<recv->idTransmitterSize(idEpoch); idTrans++)
//               if(recv->observation(idTrans, idEpoch))
//               {
//                 GnssObservationEquation eqn(*recv->observation(idTrans, idEpoch), *recv, *gnss->transmitters.at(idTrans), gnss->funcRotationCrf2Trf,
//                                             nullptr/*reduceModels*/, idEpoch, FALSE/*decorrelate*/, {}/*types*/);
//                 // spatial representation
//                 Matrix Cnm, Snm;
//                 const Vector3d point = file.rotary(eqn.timeRecv).rotate(intersection(eqn.posRecv, eqn.posTrans, eqn.elevationRecvLocal));
//                 SphericalHarmonics::CnmSnm(normalize(point), maxDegree, Cnm, Snm);
//                 Double VTEC = 0;
//                 UInt count = 0;
//                 for(UInt n=0; n<=maxDegree; n++)
//                 {
//                   VTEC += Cnm(n,0) * xt(count++);
//                   for(UInt m=1; m<=n; m++)
//                   {
//                     VTEC += Cnm(n,m) * xt(count++);
//                     VTEC += Snm(n,m) * xt(count++);
//                   }
//                 }
//                 // VTEC -> STEC
//                 recv->observation(idTrans, idEpoch)->STEC += mapping(eqn.elevationRecvLocal) * VTEC;
//               } // for(recv, trans)
//       } // for(idEpoch)
//     } // if(fileNameIn)
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssParametrizationIonosphereMap::initParameter(GnssNormalEquationInfo &normalEquationInfo)
{
  try
  {
    index.clear();
    if(!isEnabled(normalEquationInfo, name) || normalEquationInfo.isEachReceiverSeparately)
      return;

    std::vector<ParameterName> baseNames;
    for(UInt n=0; n<=maxDegree; n++)
    {
      baseNames.push_back(ParameterName("VTEC", "sphericalHarmonics.c_"+n%"%i"s+"_"+0%"%i"s));
      for(UInt m=1; m<=n; m++)
      {
        baseNames.push_back(ParameterName("VTEC", "sphericalHarmonics.c_"+n%"%i"s+"_"+m%"%i"s));
        baseNames.push_back(ParameterName("VTEC", "sphericalHarmonics.s_"+n%"%i"s+"_"+m%"%i"s));
      }
    }

    std::vector<ParameterName> temporalNames;
    temporal->parameterName(temporalNames);
    for(UInt i=0; i<temporalNames.size(); i++)
    {
      std::vector<ParameterName> parameterNames = baseNames;
      for(ParameterName &name : parameterNames)
        name.combine(temporalNames.at(i));
      index.push_back(normalEquationInfo.parameterNamesOther(parameterNames));
    }

    logInfo<<(temporalNames.size()*baseNames.size())%"%9i VTEC map parameters"s<<Log::endl;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssParametrizationIonosphereMap::aprioriParameter(const GnssNormalEquationInfo &normalEquationInfo, MatrixSliceRef x0) const
{
  try
  {
    if(Parallel::isMaster(normalEquationInfo.comm))
      for(UInt i=0; i<index.size(); i++)
        copy(x.at(i), x0.row(normalEquationInfo.index(index.at(i)), x.at(i).rows()));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

// intersection point in ionosphere height
// classic approach is to use ~450/500 km, not applicable for satellites flying possibly higher
// Hence, satellite position + 50 km, no evidence if true, but LEOs are in principle flying
// IN the ionosphere and so pierce point is set slightly higher
Vector3d GnssParametrizationIonosphereMap::intersection(const Vector3d &posRecv, const Vector3d &posTrans, Angle elevation) const
{
  try
  {
    const Double r = posRecv.r();
    const Double H = std::max(mapR+mapH, r + 50e3);  // single layer ionosphere in 450 km or 50 km above satellite
    return posRecv - normalize(posRecv-posTrans) * ((H * std::sin(PI/2 - elevation - std::asin((r/H)*std::cos(elevation)))) / std::cos(elevation));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

// Mapping function
Double GnssParametrizationIonosphereMap::mapping(Angle elevation) const
{
  return 1./std::cos(std::asin(mapR/(mapR+mapH) * std::sin(mapAlpha*(PI/2-elevation))));
}

/***********************************************/

void GnssParametrizationIonosphereMap::designMatrix(const GnssNormalEquationInfo &/*normalEquationInfo*/, const GnssObservationEquation &eqn, GnssDesignMatrix &A) const
{
  try
  {
    if(!index.size() || !selectedReceivers.at(eqn.receiver->idRecv()))
      return;

    // spatial representation
    Matrix Cnm, Snm;
    const Rotary3d rot   = magnetosphere->rotaryCelestial2SolarGeomagneticFrame(eqn.timeRecv);
    const Vector3d point = rot.rotate(intersection(eqn.posRecv, eqn.posTrans, eqn.elevationRecvLocal));
    SphericalHarmonics::CnmSnm(normalize(point), maxDegree, Cnm, Snm);
    Matrix Ynm(1, (maxDegree+1)*(maxDegree+1));
    UInt count = 0;
    for(UInt n=0; n<=maxDegree; n++)
    {
      Ynm(0, count++) = Cnm(n,0);
      for(UInt m=1; m<=n; m++)
      {
        Ynm(0, count++) = Cnm(n,m);
        Ynm(0, count++) = Snm(n,m);
      }
    }

    // VTEC -> STEC
    Matrix B(eqn.l.rows(), Ynm.columns());
    matMult(mapping(eqn.elevationRecvLocal), eqn.A.column(GnssObservationEquation::idxSTEC), Ynm, B);

    // temporal representation
    std::vector<UInt>   idx;
    std::vector<Double> factor;
    temporal->factors(eqn.timeRecv, idx, factor);
    for(UInt i=0; i<factor.size(); i++)
      axpy(factor.at(i), B, A.column(index.at(idx.at(i))));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Double GnssParametrizationIonosphereMap::updateParameter(const GnssNormalEquationInfo &normalEquationInfo, const_MatrixSliceRef x, const_MatrixSliceRef /*Wz*/)
{
  try
  {
    if(!index.size())
      return 0;

    // update parameters
    // -----------------
    for(UInt i=0; i<index.size(); i++)
      this->x.at(i) += x.row(normalEquationInfo.index(index.at(i)), this->x.at(i).rows());

    // update STEC
    // -----------
    Gnss::InfoParameterChange info("tec");
    for(UInt idEpoch=0; idEpoch<gnss->times.size(); idEpoch++)
    {
      // temporal representation
      std::vector<UInt>   idx;
      std::vector<Double> factor;
      Vector dx((maxDegree+1)*(maxDegree+1));
      temporal->factors(gnss->times.at(idEpoch), idx, factor);
      for(UInt i=0; i<factor.size(); i++)
        axpy(factor.at(i), x.row(normalEquationInfo.index(index.at(idx.at(i))), dx.rows()), dx);

      for(auto recv : gnss->receivers)
        if(recv->isMyRank() && recv->useable(idEpoch) && selectedReceivers.at(recv->idRecv()) && normalEquationInfo.estimateReceiver.at(recv->idRecv()))
          for(UInt idTrans=0; idTrans<recv->idTransmitterSize(idEpoch); idTrans++)
            if(recv->observation(idTrans, idEpoch))
            {
              GnssObservationEquation eqn(*recv->observation(idTrans, idEpoch), *recv, *gnss->transmitters.at(idTrans), gnss->funcRotationCrf2Trf,
                                          nullptr/*reduceModels*/, idEpoch, FALSE/*decorrelate*/, {}/*types*/);

              // spatial representation
              Matrix Cnm, Snm;
              const Rotary3d rot   = magnetosphere->rotaryCelestial2SolarGeomagneticFrame(eqn.timeRecv);
              const Vector3d point = rot.rotate(intersection(eqn.posRecv, eqn.posTrans, eqn.elevationRecvLocal));
              SphericalHarmonics::CnmSnm(normalize(point), maxDegree, Cnm, Snm);
              Double dVTEC = 0;
              UInt count = 0;
              for(UInt n=0; n<=maxDegree; n++)
              {
                dVTEC += Cnm(n,0) * dx(count++);
                for(UInt m=1; m<=n; m++)
                {
                   dVTEC += Cnm(n,m) * dx(count++);
                   dVTEC += Snm(n,m) * dx(count++);
                }
              }

              recv->observation(idTrans, idEpoch)->STEC += mapping(eqn.elevationRecvLocal) * dVTEC;
              if(info.update(dVTEC))
                info.info = "VTEC map ("+recv->name()+", "+gnss->times.at(idEpoch).dateTimeStr()+")";
            } // for(recv, trans)
    } // for(idEpoch)

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

void GnssParametrizationIonosphereMap::writeResults(const GnssNormalEquationInfo &normalEquationInfo, const std::string &suffix) const
{
  try
  {
    if(!isEnabled(normalEquationInfo, name) || !Parallel::isMaster(normalEquationInfo.comm) || fileNameOut.empty() || !index.size())
      return;

    logStatus<<"write ionosphere maps to files <"<<fileNameOut.appendBaseName(suffix)<<">"<<Log::endl;

//     // determine sampling times
//     const Time timeStart = gnss->times.front();
//     const Time timeEnd   = gnss->times.back();
//     std::vector<Time> times(temporal->parameterCount(), timeStart);
//     for(UInt i=1; i<temporal->parameterCount(); i++)
//       times.at(i) += i/(times.size()-1.) * (timeEnd-timeStart);
//
//     // interpolate temporal spherical harmonics
//     std::vector<Vector3d> magneticNorthPoles(times.size());
//     std::vector<Matrix>   cnm(times.size(), Matrix(maxDegree+1, Matrix::TRIANGULAR, Matrix::LOWER));
//     std::vector<Matrix>   snm(times.size(), Matrix(maxDegree+1, Matrix::TRIANGULAR, Matrix::LOWER));
//     for(UInt idEpoch=0; idEpoch<times.size(); idEpoch++)
//     {
//       magneticNorthPoles.at(idEpoch) = magnetosphere->geomagneticNorthPole(times.at(idEpoch));
//       // temporal representation
//       Vector xt((maxDegree+1)*(maxDegree+1));
//       std::vector<UInt>   idx;
//       std::vector<Double> factor;
//       temporal->factors(std::min(times.at(idEpoch), times.back()-seconds2time(0.1)), idx, factor);
//       for(UInt i=0; i<factor.size(); i++)
//         axpy(factor.at(i), x.at(idx.at(i)), xt);
//
//       // sort into spherical harmonics
//       UInt count = 0;
//       for(UInt n=0; n<=maxDegree; n++)
//       {
//         cnm.at(idEpoch)(n,0) = xt(count++);
//         for(UInt m=1; m<=n; m++)
//         {
//           cnm.at(idEpoch)(n,m) = xt(count++);
//           snm.at(idEpoch)(n,m) = xt(count++);
//         }
//       }
//     } // for(idEpoch)
//
//     writeFileGnssIonosphereMaps(fileNameOut.appendBaseName(suffix), times, magneticNorthPoles, cnm, snm);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
