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
#include "base/basisSplines.h"
#include "base/sphericalHarmonics.h"
#include "config/config.h"
#include "files/fileGriddedDataTimeSeries.h"
#include "classes/grid/grid.h"
#include "classes/timeSeries/timeSeries.h"
#include "classes/magnetosphere/magnetosphere.h"
#include "gnss/gnssParametrization/gnssParametrizationIonosphereMap.h"

/***********************************************/

GnssParametrizationIonosphereMap::GnssParametrizationIonosphereMap(Config &config)
{
  try
  {
    GridPtr       gridPtr;
    TimeSeriesPtr timeSeriesPtr;

    readConfig(config, "name",                            name,            Config::OPTIONAL, "parameter.mapVTEC", "");
    readConfig(config, "selectReceivers",                 selectReceivers, Config::MUSTSET,  R"(["all"])",        "");
    readConfig(config, "outputfileGriddedDataTimeSeries", fileNameOut,     Config::OPTIONAL, "",                  "single layer VTEC [TECU]");
    readConfig(config, "outputGrid",                      gridPtr,         Config::DEFAULT,  R"({"geograph":{"deltaLambda":"5","deltaPhi":"2.5","height":"450e3","R":"6371e3","inverseFlattening":0}})", "");
    readConfig(config, "outputTimeSeries",                timeSeriesPtr,   Config::DEFAULT,  R"({"uniformSampling":{"sampling":"2/24"}})", "");
    readConfig(config, "inputfileGriddedDataTimeSeries",  fileNameIn,      Config::OPTIONAL, "",                  "single layer VTEC [TECU]");
    readConfig(config, "maxDegree",                       maxDegree,       Config::MUSTSET,  "15",                "spherical harmonics parametrization");
    readConfig(config, "temporal",                        temporal,        Config::DEFAULT,  "",                  "temporal evolution of VTEC values");
    readConfig(config, "radiusIonosphericLayer",          radiusIono,      Config::DEFAULT,  "6371e3+450e3",      "[m] radius of ionospheric single layer");
    readConfig(config, "mapR",                            mapR,            Config::DEFAULT,  "6371e3",            "[m] constant of MSLM mapping function");
    readConfig(config, "mapH",                            mapH,            Config::DEFAULT,  "506.7e3",           "[m] constant of MSLM mapping function");
    readConfig(config, "mapAlpha",                        mapAlpha,        Config::DEFAULT,  "0.9782",            "constant of MSLM mapping function");
    readConfig(config, "magnetosphere",                   magnetosphere,   Config::MUSTSET,  "",                  "");
    if(isCreateSchema(config)) return;

    gridOut  = GriddedData(Ellipsoid(mapR, 0), gridPtr->points(), gridPtr->areas(), {});
    timesOut = timeSeriesPtr->times();
    if(!fileNameOut.empty() && !(gridOut.points.size() && timesOut.size()))
      throw(Exception("outputfileGriddedDataTimeSeries needs outputGrid and outputTimeSeries"));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

static Double interpolateGrid(const Vector3d &point, const GriddedDataRectangular &grid, const Vector &data)
{
  try
  {
    // latitude
    Double tauLat  = Double(point.phi()-grid.latitudes.at(0))/Double(grid.latitudes.at(1)-grid.latitudes.at(0));
    const UInt idxLat1 = std::min(static_cast<UInt>(std::max(std::floor(tauLat), 0.)), grid.latitudes.size()-1);
    const UInt idxLat2 = std::min(idxLat1+1, grid.latitudes.size()-1);
    tauLat -= std::floor(tauLat);
    // longitude
    Double tauLon  = Double(point.lambda()-grid.longitudes.at(0))/Double(grid.longitudes.at(1)-grid.longitudes.at(0));
    const UInt idxLon1 = static_cast<UInt>(std::floor(tauLon)+grid.longitudes.size())%grid.longitudes.size();
    const UInt idxLon2 = (idxLon1+1)%grid.longitudes.size();
    tauLon -= std::floor(tauLon);

    return (1-tauLon) * (1-tauLat) * data(idxLon1 + grid.longitudes.size() * idxLat1)
         + (1-tauLon) *  (tauLat)  * data(idxLon1 + grid.longitudes.size() * idxLat2)
         +  (tauLon)  * (1-tauLat) * data(idxLon2 + grid.longitudes.size() * idxLat1)
         +  (tauLon)  *  (tauLat)  * data(idxLon2 + grid.longitudes.size() * idxLat2);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssParametrizationIonosphereMap::init(Gnss *gnss, Parallel::CommunicatorPtr comm)
{
  try
  {
    this->gnss = gnss;
    selectedReceivers = gnss->selectReceivers(selectReceivers);
    x.resize(temporal->parameterCount(), Vector((maxDegree+1)*(maxDegree+1)));

    // apriori model
    // -------------
    if(!fileNameIn.empty())
    {
      logStatus<<"read apriori ionosphere <"<<fileNameIn<<">"<<Log::endl;
      InFileGriddedDataTimeSeries file(fileNameIn);
      GriddedDataRectangular grid;
      if(!grid.init(file.grid()))
        throw(Exception(fileNameIn.str()+" must contain a rectangular grid"));
      if((gnss->times.front() < file.times().front()) || (gnss->times.back() > file.times().back()))
        throw(Exception("wrong time interval ["+file.times().front().dateTimeStr()+", "+file.times().back().dateTimeStr()+"] in <"+fileNameIn.str()+">"));

      UInt indexData = file.nodeCount();
      std::vector<Vector>   griddedVTEC(file.splineDegree()+1);
      std::vector<Rotary3d> rotSmf2Trf(file.splineDegree()+1); // SolarGeomagneticFrame ->  Terrestrial

      Log::Timer timer(gnss->times.size());
      for(UInt idEpoch=0; idEpoch<gnss->times.size(); idEpoch++)
      {
        timer.loopStep(idEpoch);

        // find time interval and read missing nodes
        const UInt idx = std::min(static_cast<UInt>(std::distance(file.times().begin(),
                                                                  std::upper_bound(file.times().begin(), file.times().end(), gnss->times.at(idEpoch)))), file.times().size()-1)-1;
        const UInt start = (idx > indexData) ? std::max(idx, indexData+griddedVTEC.size()) : idx;
        const UInt end   = (idx > indexData) ? idx+griddedVTEC.size() : std::min(idx+griddedVTEC.size(), indexData);
        for(UInt i=start; i<end; i++)
        {
          griddedVTEC.at(i%griddedVTEC.size()) = file.data(i).column(0);
          rotSmf2Trf.at(i%griddedVTEC.size())  = gnss->rotationCrf2Trf(file.times().at(i)) * inverse(magnetosphere->rotaryCelestial2SolarGeomagneticFrame(file.times().at(i)));
        }
        indexData = idx;

        // prepare spline interpolation
        const Double t     = (gnss->times.at(idEpoch)-file.times().at(idx)).mjd()/(file.times().at(idx+1)-file.times().at(idx)).mjd();
        const Vector coeff = BasisSplines::compute(t, file.splineDegree());

        // Celestial -> SolarGeomagneticFrame
        const Rotary3d rotCrf2Smf = magnetosphere->rotaryCelestial2SolarGeomagneticFrame(gnss->times.at(idEpoch));

        for(auto recv : gnss->receivers)
          if(recv->isMyRank() && recv->useable(idEpoch) && selectedReceivers.at(recv->idRecv()))
            for(UInt idTrans=0; idTrans<recv->idTransmitterSize(idEpoch); idTrans++)
              if(recv->observation(idTrans, idEpoch))
              {

                GnssObservationEquation eqn(*recv->observation(idTrans, idEpoch), *recv, *gnss->transmitters.at(idTrans), gnss->funcRotationCrf2Trf,
                                            nullptr/*reduceModels*/, idEpoch, FALSE/*homogenize*/, {}/*types*/);
                // pierce point in SolarGeomagneticFrame
                const Vector3d point = rotCrf2Smf.rotate(intersection(radiusIono, eqn.posRecv, eqn.posTrans));
                // spline interpolation
                Double VTEC = 0;
                for(UInt i=0; i<coeff.rows(); i++)
                  VTEC += coeff(i) * interpolateGrid(rotSmf2Trf.at((indexData+i)%griddedVTEC.size()).rotate(point), grid, griddedVTEC.at((indexData+i)%griddedVTEC.size()));
                // mapping VTEC -> STEC
                recv->observation(idTrans, idEpoch)->STEC += mapping(eqn.elevationRecvLocal) * VTEC;
              } // for(recv, trans)
      } // for(idEpoch)
      Parallel::barrier(comm);
      timer.loopEnd();
    } // if(fileNameIn)
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

    if(temporalNames.size() && baseNames.size())
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
Vector3d GnssParametrizationIonosphereMap::intersection(const Double radiusIono, const Vector3d &posRecv, const Vector3d &posTrans) const
{
  try
  {
    const Double   rRecv = std::min(posRecv.r(), radiusIono); // LEO satellites flying possibly higher
    const Vector3d k     = normalize(posRecv-posTrans);       // direction from transmitter
    const Double   rk    = inner(posRecv, k);
    return posRecv - (std::sqrt(rk*rk+radiusIono*radiusIono-rRecv*rRecv)+rk) * k;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Double GnssParametrizationIonosphereMap::mapping(Angle elevation) const
{
  return 1./std::cos(std::asin(mapR/(mapR+mapH) * std::sin(mapAlpha*(PI/2-elevation))));
}

/***********************************************/

Double GnssParametrizationIonosphereMap::sphericalHarmonicSynthesis(const Vector3d &point, const Vector &x) const
{
  try
  {
    Matrix Cnm, Snm;
    SphericalHarmonics::CnmSnm(normalize(point), maxDegree, Cnm, Snm);
    Double VTEC  = 0;
    UInt   count = 0;
    for(UInt n=0; n<=maxDegree; n++)
    {
      VTEC += Cnm(n,0) * x(count++);
      for(UInt m=1; m<=n; m++)
      {
        VTEC += Cnm(n,m) * x(count++);
        VTEC += Snm(n,m) * x(count++);
      }
    }
    return VTEC;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
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
    const Vector3d point = rot.rotate(intersection(radiusIono, eqn.posRecv, eqn.posTrans));
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
                                          nullptr/*reduceModels*/, idEpoch, FALSE/*homogenize*/, {}/*types*/);
              // spatial representation
              const Rotary3d rot   = magnetosphere->rotaryCelestial2SolarGeomagneticFrame(eqn.timeRecv);
              const Double   dVTEC = sphericalHarmonicSynthesis(rot.rotate(intersection(radiusIono, eqn.posRecv, eqn.posTrans)), dx);
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
    if(!isEnabled(normalEquationInfo, name) || !Parallel::isMaster(normalEquationInfo.comm) || fileNameOut.empty())
      return;

    logStatus<<"write ionosphere gridded time series to file <"<<fileNameOut.appendBaseName(suffix)<<">"<<Log::endl;
    std::vector<Matrix> data(timesOut.size(), Vector(gridOut.points.size()));

    // apriori model
    // -------------
    if(!fileNameIn.empty())
    {
      InFileGriddedDataTimeSeries file(fileNameIn);
      GriddedDataRectangular gridFile;
      if(!gridFile.init(file.grid()))
        throw(Exception(fileNameIn.str()+" must contain a rectangular grid"));
      if((timesOut.front() < file.times().front()) || (timesOut.back() > file.times().back()))
        throw(Exception("wrong time interval ["+file.times().front().dateTimeStr()+", "+file.times().back().dateTimeStr()+"] in <"+fileNameIn.str()+">"));

      UInt indexData = file.nodeCount();
      std::vector<Vector>   griddedVTEC(file.splineDegree()+1);
      std::vector<Rotary3d> rotSmf2Trf(file.splineDegree()+1); // SolarGeomagneticFrame -> Terrestrial

      for(UInt idEpoch=0; idEpoch<timesOut.size(); idEpoch++)
      {
        // find time interval and read missing nodes
        const UInt idx = std::min(static_cast<UInt>(std::distance(file.times().begin(),
                                                                  std::upper_bound(file.times().begin(), file.times().end(), timesOut.at(idEpoch)))), file.times().size()-1)-1;
        const UInt start = (idx > indexData) ? std::max(idx, indexData+griddedVTEC.size()) : idx;
        const UInt end   = (idx > indexData) ? idx+griddedVTEC.size() : std::min(idx+griddedVTEC.size(), indexData);
        for(UInt i=start; i<end; i++)
        {
          griddedVTEC.at(i%griddedVTEC.size()) = file.data(i).column(0);
          rotSmf2Trf.at(i%griddedVTEC.size())  = gnss->rotationCrf2Trf(file.times().at(i)) * inverse(magnetosphere->rotaryCelestial2SolarGeomagneticFrame(file.times().at(i)));
        }
        indexData = idx;

        // prepare spline interpolation
        const Double t     = (timesOut.at(idEpoch)-file.times().at(idx)).mjd()/(file.times().at(idx+1)-file.times().at(idx)).mjd();
        const Vector coeff = BasisSplines::compute(t, file.splineDegree());

        // Terrestrial -> SolarGeomagneticFrame
        const Rotary3d rotTrf2Smf = magnetosphere->rotaryCelestial2SolarGeomagneticFrame(timesOut.at(idEpoch)) * inverse(gnss->rotationCrf2Trf(timesOut.at(idEpoch)));

        for(UInt i=0; i<gridOut.points.size(); i++)
        {
          Vector3d point = rotTrf2Smf.rotate(gridOut.points.at(i));
          Double VTEC = 0;
          for(UInt i=0; i<coeff.rows(); i++)
            VTEC += coeff(i) * interpolateGrid(rotSmf2Trf.at((indexData+i)%griddedVTEC.size()).rotate(point), gridFile, griddedVTEC.at((indexData+i)%griddedVTEC.size()));
          data.at(idEpoch)(i, 0) += VTEC;
        }
      }
    } // if(!fileNameIn.empty())

    // estimated VTEC maps
    // -------------------
    if(x.size())
      for(UInt idEpoch=0; idEpoch<timesOut.size(); idEpoch++)
      {
        // temporal representation
        std::vector<UInt>   idx;
        std::vector<Double> factor;
        Vector xt((maxDegree+1)*(maxDegree+1));
        temporal->factors(timesOut.at(idEpoch), idx, factor);
        for(UInt i=0; i<factor.size(); i++)
          axpy(factor.at(i), x.at(idx.at(i)), xt);
        const Rotary3d rot = magnetosphere->rotaryCelestial2SolarGeomagneticFrame(timesOut.at(idEpoch))
                           * inverse(gnss->rotationCrf2Trf(timesOut.at(idEpoch))); // mag frame -> TRF
        for(UInt i=0; i<gridOut.points.size(); i++)
          data.at(idEpoch)(i, 0) += sphericalHarmonicSynthesis(rot.rotate(gridOut.points.at(i)), xt);
      } // for(idEpoch)

    writeFileGriddedDataTimeSeries(fileNameOut.appendBaseName(suffix), 1, timesOut, gridOut, data);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
