/***********************************************/
/**
* @file gnssParametrizationEarthRotation.cpp
*
* @brief EarthRotation.
* @see GnssParametrization
*
* @author Sebastian Strasser
* @author Torsten Mayer-Guerr
* @date 2021-01-23
*
*/
/***********************************************/

#include "base/import.h"
#include "config/config.h"
#include "base/planets.h"
#include "files/fileInstrument.h"
#include "classes/parametrizationTemporal/parametrizationTemporal.h"
#include "gnss/gnssParametrization/gnssParametrizationEarthRotation.h"

/***** CLASS ***********************************/

/** @brief EarthRotation.
* @ingroup gnssParametrizationGroup
* @see GnssParametrization */

/***********************************************/

GnssParametrizationEarthRotation::GnssParametrizationEarthRotation(Config &config)
{
  try
  {
    readConfig(config, "name",             name,                    Config::OPTIONAL, "parameter.earthRotation", "used for parameter selection");
    readConfig(config, "outputfileEOP",    fileNameEOP,             Config::OPTIONAL, "", "EOP time series (mjd, xp, yp, sp, dUT1, LOD, X, Y, S)");
    readConfig(config, "estimatePole",     parametrizationPole,     Config::DEFAULT,  "", "xp, yp [mas]");
    readConfig(config, "estimateUT1",      parametrizationUT1,      Config::DEFAULT,  "", "rotation angle [ms]");
    readConfig(config, "estimateNutation", parametrizationNutation, Config::DEFAULT,  "", "dX, dY [mas]");
    if(isCreateSchema(config)) return;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssParametrizationEarthRotation::init(Gnss *gnss, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    this->gnss = gnss;
    xPole     = Vector(2*parametrizationPole->parameterCount());
    xUT1      = Vector(parametrizationUT1->parameterCount());
    xNutation = Vector(2*parametrizationNutation->parameterCount());
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssParametrizationEarthRotation::initParameter(GnssNormalEquationInfo &normalEquationInfo)
{
  try
  {
    indexParameterPole     = GnssParameterIndex();
    indexParameterUT1      = GnssParameterIndex();
    indexParameterNutation = GnssParameterIndex();

    if(!isEnabled(normalEquationInfo, name) || normalEquationInfo.isEachReceiverSeparately)
      return;
    UInt countPara = 0;

    // Polar motion
    if(parametrizationPole->parameterCount())
    {
      std::vector<ParameterName> parameterNames;
      parametrizationPole->parameterName({ParameterName("earth", "polarMotion.xp"), ParameterName("earth", "polarMotion.yp")}, parameterNames);
      indexParameterPole = normalEquationInfo.parameterNamesOther(parameterNames);
      countPara += parameterNames.size();
    }

    // Earth rotation angle (UT1)
    if(parametrizationUT1->parameterCount())
    {
      std::vector<ParameterName> parameterNames;
      parametrizationUT1->parameterName({ParameterName("earth", "UT1")}, parameterNames);
      indexParameterUT1 = normalEquationInfo.parameterNamesOther(parameterNames);
      countPara += parameterNames.size();
    }

    // Nutation
    if(parametrizationNutation->parameterCount())
    {
      std::vector<ParameterName> parameterNames;
      parametrizationNutation->parameterName({ParameterName("earth", "nutation.X"), ParameterName("earth", "nutation.Y")}, parameterNames);
      indexParameterNutation = normalEquationInfo.parameterNamesOther(parameterNames);
      countPara += parameterNames.size();
    }

    if(countPara)
      logInfo<<countPara%"%9i Earth rotation parameters"s<<Log::endl;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssParametrizationEarthRotation::aprioriParameter(const GnssNormalEquationInfo &normalEquationInfo, MatrixSliceRef x0) const
{
  try
  {
    if(!Parallel::isMaster(normalEquationInfo.comm))
      return;
    if(indexParameterPole)
      copy(xPole, x0.row(normalEquationInfo.index(indexParameterPole), xPole.rows()));
    if(indexParameterUT1)
      copy(xUT1, x0.row(normalEquationInfo.index(indexParameterUT1), xUT1.rows()));
    if(indexParameterNutation)
      copy(xNutation, x0.row(normalEquationInfo.index(indexParameterNutation), xNutation.rows()));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssParametrizationEarthRotation::designMatrix(const GnssNormalEquationInfo &/*normalEquationInfo*/, const GnssObservationEquation &eqn, GnssDesignMatrix &A) const
{
  try
  {
    if(!eqn.receiver->isEarthFixed())
      return;
    if(!(indexParameterPole || indexParameterUT1 || indexParameterNutation))
      return;

    const Double xp      = gnss->eop(eqn.idEpoch, 0);
    const Double yp      = gnss->eop(eqn.idEpoch, 1);
    const Double sp      = gnss->eop(eqn.idEpoch, 2);
    const Double deltaUT = gnss->eop(eqn.idEpoch, 3) + (eqn.timeRecv-timeGPS2UTC(eqn.timeRecv)).seconds();
    const Double X       = gnss->eop(eqn.idEpoch, 5);
    const Double Y       = gnss->eop(eqn.idEpoch, 6);
    const Double S       = gnss->eop(eqn.idEpoch, 7);
    const Double ERA     = Planets::ERA(timeGPS2UTC(eqn.timeRecv) + seconds2time(deltaUT));

    const Matrix rotW   = (rotaryY(Angle(xp)) * rotaryX(Angle(yp))).matrix(); // Polar motion rotation matrix
    const Matrix rotERA = rotaryZ(Angle(S-ERA-sp)).matrix();   // Earth rotation angle rotation matrix
    const Double a      = 0.5 + 0.125*(X*X+Y*Y);
    const Double X2     = X*X;
    const Double Y2     = Y*Y;
    const Matrix rotQ({{1.-a*X2,   -a*X*Y,             X},    // Precession & nutation rotation matrix
                       {  -a*X*Y, 1.-a*Y2,             Y},
                       {      -X,      -Y, 1.-a*(X2+Y2)}});

    const Vector posEarth     = gnss->rotationCrf2Trf(eqn.timeRecv).rotate(eqn.posRecv).vector();
    const Time   timeTemporal = std::max(eqn.timeRecv, gnss->times.at(0));

    // temporal parametrization
    auto designMatrixTemporal = [&](ParametrizationTemporalPtr parametrization, const_MatrixSliceRef B, const GnssParameterIndex &index)
    {
      std::vector<UInt>   idx;
      std::vector<Double> factor;
      parametrization->factors(timeTemporal, idx, factor);
      MatrixSlice Design(A.column(index));
      for(UInt i=0; i<factor.size(); i++)
        axpy(factor.at(i), B, Design.column(B.columns()*idx.at(i), B.columns()));
    };

    // Polar motion
    // ------------
    if(indexParameterPole)
    {
      // partial derivatives of position(CRF) with respect to xp, yp
      Matrix partWxp(3,3), partWyp(3,3);
      partWxp(0,2) = partWyp(2,1) = -1;
      partWxp(2,0) = partWyp(1,2) =  1;

      Matrix B(3,2);
      matMult(DEG2RAD/3600e3, partWxp, posEarth, B.column(0));
      matMult(DEG2RAD/3600e3, partWyp, posEarth, B.column(1));
      B = eqn.A.column(GnssObservationEquation::idxPosRecv,3) * (rotQ * (rotERA * B));
      designMatrixTemporal(parametrizationPole, B, indexParameterPole);
    }

    // Earth rotation angle (UT1)
    // --------------------------
    if(indexParameterUT1)
    {
      // partial derivatives of position(CRF) with respect to UT1
      Matrix partERA({{std::sin(S-ERA-sp), -std::cos(S-ERA-sp),  0},
                      {std::cos(S-ERA-sp),  std::sin(S-ERA-sp),  0},
                      {                 0,                   0,  0}});

      Matrix B = 1e-3 * 2*PI*1.00273781191135448/86400. * (eqn.A.column(GnssObservationEquation::idxPosRecv,3) * (rotQ * (partERA * (rotW * posEarth))));
      designMatrixTemporal(parametrizationUT1, B, indexParameterUT1);
    }

    // Nutation
    // --------
    if(indexParameterNutation)
    {
      // partial derivatives of position(CRF) with respect to X, Y
      Matrix partQX({{-2*X*a-X2*X/4,   -Y*a-Y*X2/4,                  1},
                     {  -Y*a-X2*Y/4,       -X*Y2/4,                  0},
                     {           -1,             0, -2*X*a-X*(X2+Y2)/4}});
      Matrix partQY({{      -X2*Y/4,   -X*a-X*Y2/4,                  0},
                     {  -X*a-X*Y2/4, -2*Y*a-Y*Y2/4,                  1},
                     {            0,            -1, -2*Y*a-Y*(X2+Y2)/4}});

      Matrix B(3,2);
      Matrix p = rotERA * (rotW * posEarth);
      matMult(DEG2RAD/3600e3, partQX, p, B.column(0));
      matMult(DEG2RAD/3600e3, partQY, p, B.column(1));
      B = eqn.A.column(GnssObservationEquation::idxPosRecv,3) * B;
      designMatrixTemporal(parametrizationNutation, B, indexParameterNutation);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Double GnssParametrizationEarthRotation::updateParameter(const GnssNormalEquationInfo &normalEquationInfo, const_MatrixSliceRef x, const_MatrixSliceRef /*Wz*/)
{
  try
  {
    Double maxChange = 0;

    // Polar motion
    // ------------
    if(indexParameterPole)
    {
      Gnss::InfoParameterChange info("mm");
      const Vector dx = x.row(normalEquationInfo.index(indexParameterPole), 2*parametrizationPole->parameterCount());
      xPole += dx;
      std::vector<UInt>   index;
      std::vector<Double> factor;
      for(UInt idEpoch=0; idEpoch<gnss->times.size(); idEpoch++)
      {
        parametrizationPole->factors(gnss->times.at(idEpoch), index, factor);
        Vector p(2);
        for(UInt k=0; k<factor.size(); k++)
          axpy(DEG2RAD/3600e3 * factor.at(k), dx.row(2*index.at(k),2), p);

        // update xp,yp
        gnss->eop(idEpoch, 0) += p(0);
        gnss->eop(idEpoch, 1) += p(1);
        if(info.update(1e3*DEFAULT_R*norm(p)))
          info.info = "earth rotation (pole)";
      } // for(idEpoch)
      info.synchronizeAndPrint(normalEquationInfo.comm, 1e-3, maxChange);
    } // if(ParameterPole)

    // Earth Rotation Angle (UT1)
    // --------------------------
    if(indexParameterUT1)
    {
      Gnss::InfoParameterChange info("mm");
      const Vector dx = x.row(normalEquationInfo.index(indexParameterUT1), parametrizationUT1->parameterCount());
      xUT1 += dx;
      std::vector<UInt>   index;
      std::vector<Double> factor;
      for(UInt idEpoch=0; idEpoch<gnss->times.size(); idEpoch++)
      {
        parametrizationUT1->factors(gnss->times.at(idEpoch), index, factor);
        Double p = 0;
        for(UInt k=0; k<factor.size(); k++)
          p += 1e-3 * factor.at(k) * dx(index.at(k));

        gnss->eop(idEpoch, 3) += p;
        if(info.update(1e3*DEFAULT_R*2*PI*1.00273781191135448/86400*p))
          info.info = "earth rotation (UT1)";
      } // for(idEpoch)
      info.synchronizeAndPrint(normalEquationInfo.comm, 1e-3, maxChange);
    } // if(ParameterUT1)

    // Nutation
    // --------
    if(indexParameterNutation)
    {
      Gnss::InfoParameterChange info("mm");
      const Vector dx = x.row(normalEquationInfo.index(indexParameterNutation), 2*parametrizationNutation->parameterCount());
      xNutation += dx;
      std::vector<UInt>   index;
      std::vector<Double> factor;
      for(UInt idEpoch=0; idEpoch<gnss->times.size(); idEpoch++)
      {
        parametrizationNutation->factors(gnss->times.at(idEpoch), index, factor);
        Vector p(2);
        for(UInt k=0; k<factor.size(); k++)
          axpy(DEG2RAD/3600e3 * factor.at(k), dx.row(2*index.at(k),2), p);

        // update X,Y
        gnss->eop(idEpoch, 5) += p(0);
        gnss->eop(idEpoch, 6) += p(1);
        if(info.update(1e3*DEFAULT_R*norm(p)))
          info.info = "earth rotation (XY)";
      } // for(idEpoch)
      info.synchronizeAndPrint(normalEquationInfo.comm, 1e-3, maxChange);
    } // if(ParameterNutation)

    return maxChange;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssParametrizationEarthRotation::writeResults(const GnssNormalEquationInfo &normalEquationInfo, const std::string &suffix) const
{
  try
  {
    if(!isEnabled(normalEquationInfo, name) || !Parallel::isMaster(normalEquationInfo.comm))
      return;

    if(!fileNameEOP.empty())
    {
      logStatus<<"write EOP to files <"<<fileNameEOP.appendBaseName(suffix)<<">"<<Log::endl;
      std::vector<Time> times;
      Matrix A(normalEquationInfo.idEpochs.size(), 9);
      for(UInt i=0; i<normalEquationInfo.idEpochs.size(); i++)
      {
        times.push_back(gnss->times.at(normalEquationInfo.idEpochs.at(i)));
        copy(gnss->eop.row(i), A.slice(i, 1, 1, 8));
        A(i, 4) += (gnss->times.at(normalEquationInfo.idEpochs.at(i))-timeGPS2UTC(gnss->times.at(normalEquationInfo.idEpochs.at(i)))).seconds(); // UT1-GPS => UT1-UTC
      }
      InstrumentFile::write(fileNameEOP.appendBaseName(suffix), Arc(times, A));
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
