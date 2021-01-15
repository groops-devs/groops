/***********************************************/
/**
* @file gnssParametrizationEarthRotation.cpp
*
* @brief GNSS Earth rotation.
*
* @author Torsten Mayer-Guerr
* @author Sebastian Strasser
* @date 2014-05-25
*
*/
/***********************************************/

#define DOCSTRING_GnssParametrizationEarthRotation

#include "base/import.h"
#include "base/polynomial.h"
#include "base/planets.h"
#include "parser/dataVariables.h"
#include "config/config.h"
#include "config/configRegister.h"
#include "inputOutput/logging.h"
#include "files/fileInstrument.h"
#include "classes/earthRotation/earthRotation.h"
#include "gnss/gnss.h"
#include "gnss/gnssObservation.h"
#include "gnss/gnssDesignMatrix.h"
#include "gnss/gnssParametrizationEarthRotation.h"

/***********************************************/

GROOPS_REGISTER_CLASS_WITHOUT_SUBS(GnssParametrizationEarthRotation, "gnssParametrizationEarthRotationType")
GROOPS_READCONFIG_CLASS(GnssParametrizationEarthRotation, "gnssParametrizationEarthRotationType")

/***********************************************/

GnssParametrizationEarthRotation::GnssParametrizationEarthRotation(Config &config, const std::string &name)
{
  try
  {
    UInt interpolationDegree;

    readConfigSequence(config, name, Config::MUSTSET, "", "");
    readConfig(config, "outputfileEOP",        fileNameTimeSeriesEOP,   Config::OPTIONAL, "",  "EOP time series (mjd, xp, yp, sp, dUT1, LOD, X, Y, S)");
    readConfig(config, "earthRotation",        earthRotationPtr,        Config::MUSTSET,  "",  "apriori earth rotation");
    readConfig(config, "earthRotationModels",  earthRotationModelsPtr,  Config::OPTIONAL, "",  "will be reduced from reported apriori parameters");
    readConfig(config, "estimatePole",         parametrizationPole,     Config::DEFAULT,  "",  "xp, yp [mas]");
    readConfig(config, "estimateUT1",          parametrizationUT1,      Config::DEFAULT,  "",  "rotation angle [ms]");
    readConfig(config, "estimateNutation",     parametrizationNutation, Config::DEFAULT,  "",  "dX, dY [mas]");
    readConfig(config, "interpolationDegree",  interpolationDegree,     Config::DEFAULT,  "7", "for quaternion interpolation");
    endSequence(config);
    if(isCreateSchema(config)) return;

    polynomial.init(interpolationDegree);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssParametrizationEarthRotation::initIntervalEarly(Gnss::AnalysisType /*analysisType*/, const std::vector<Time> &times, const Time &/*timeMargin*/, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    this->times = times;

    parametrizationPole->setInterval    (this->times.front(), this->times.back(), TRUE);
    parametrizationNutation->setInterval(this->times.front(), this->times.back(), TRUE);
    parametrizationUT1->setInterval     (this->times.front(), this->times.back(), TRUE);

    // Matrix eop columns: xp, yp, sp, deltaUT, LOD, X, Y, S
    // -----------------------------------------------------
    eop = Matrix(this->times.size(), 8);
    for(UInt i=0; i<this->times.size(); i++)
      earthRotationPtr->earthOrientationParameter(this->times.at(i), eop(i,0), eop(i,1), eop(i,2), eop(i,3), eop(i,4), eop(i,5), eop(i,6), eop(i,7));
    // UT1-UTC => UT1-GPS (avoid leap seconds jumps for interpolation)
    for(UInt i=0; i<this->times.size(); i++)
      eop(i,3) -= (this->times.at(i)-timeGPS2UTC(this->times.at(i))).seconds();

    // Matrix eop columns: xp, yp, sp, deltaUT, LOD, X, Y, S (only models)
    // -------------------------------------------------------------------
    eopModels = Matrix(this->times.size(), 8);
    if(earthRotationModelsPtr)
      for(UInt i=0; i<this->times.size(); i++)
        earthRotationModelsPtr->earthOrientationParameter(this->times.at(i), eopModels(i,0), eopModels(i,1), eopModels(i,2), eopModels(i,3), eopModels(i,4), eopModels(i,5), eopModels(i,6), eopModels(i,7));
    // UT1-UTC => UT1-GPS (avoid leap seconds jumps for interpolation)
    for(UInt i=0; i<this->times.size(); i++)
      eopModels(i,3) -= (this->times.at(i)-timeGPS2UTC(this->times.at(i))).seconds();

    eopInitial = eop;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssParametrizationEarthRotation::earthOrientationParameter(const Time &time, Double &xp, Double &yp, Double &sp, Double &deltaUT, Double &LOD, Double &X, Double &Y, Double &S) const
{
  try
  {
    const Matrix eopInterpolated = polynomial.interpolate({time}, times, eop); // EOP interpolation
    xp      = eopInterpolated(0,0);
    yp      = eopInterpolated(0,1);
    sp      = eopInterpolated(0,2);
    deltaUT = eopInterpolated(0,3) + (time-timeGPS2UTC(time)).seconds();
    LOD     = eopInterpolated(0,4);
    X       = eopInterpolated(0,5);
    Y       = eopInterpolated(0,6);
    S       = eopInterpolated(0,7);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Rotary3d GnssParametrizationEarthRotation::rotaryMatrix(const Time &time) const
{
  try
  {
    Double xp, yp, sp, deltaUT, LOD, X, Y, S;
    earthOrientationParameter(time, xp, yp, sp, deltaUT, LOD, X, Y, S);

    const Double ERA = Planets::ERA(timeGPS2UTC(time) + seconds2time(deltaUT));
    const Double r2  = X*X + Y*Y;
    const Double E   = (r2!=0.) ? std::atan2(Y, X) : 0.;
    const Double D   = std::atan(std::sqrt(r2/(1-r2)));

    return  rotaryX(Angle(-yp)) * rotaryY(Angle(-xp)) *
            rotaryZ(Angle(sp+ERA-S-E)) *
            rotaryY(Angle(D)) * rotaryZ(Angle(E));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssParametrizationEarthRotation::initParameter(Gnss::NormalEquationInfo &normalEquationInfo)
{
  try
  {
    indexParameterPole     = Gnss::ParameterIndex();
    indexParameterUT1      = Gnss::ParameterIndex();
    indexParameterNutation = Gnss::ParameterIndex();

    // Polar motion
    // ------------
    if(parametrizationPole->parameterCount() && (normalEquationInfo.estimationType & Gnss::NormalEquationInfo::ESTIMATE_EARTHROTATION_POLE))
    {
      std::vector<ParameterName> parameterNames;
      parametrizationPole->parameterName({ParameterName("earth", "polarMotion.xp"), ParameterName("earth", "polarMotion.yp")}, parameterNames);
      indexParameterPole = normalEquationInfo.parameterNamesOther(parameterNames);
    }

    // Earth rotation angle (UT1)
    // -------------------------
    if(parametrizationUT1->parameterCount() && (normalEquationInfo.estimationType & Gnss::NormalEquationInfo::ESTIMATE_EARTHROTATION_UT1))
    {
      std::vector<ParameterName> parameterNames;
      parametrizationUT1->parameterName({ParameterName("earth", "UT1")}, parameterNames);
      indexParameterUT1 = normalEquationInfo.parameterNamesOther(parameterNames);
    }

    // Nutation
    // --------
    if(parametrizationNutation->parameterCount() && (normalEquationInfo.estimationType & Gnss::NormalEquationInfo::ESTIMATE_EARTHROTATION_NUTATION))
    {
      std::vector<ParameterName> parameterNames;
      parametrizationNutation->parameterName({ParameterName("earth", "nutation.X"), ParameterName("earth", "nutation.Y")}, parameterNames);
      indexParameterNutation = normalEquationInfo.parameterNamesOther(parameterNames);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssParametrizationEarthRotation::aprioriParameter(const Gnss::NormalEquationInfo &normalEquationInfo, MatrixSliceRef x0) const
{
  try
  {
    if(!Parallel::isMaster(normalEquationInfo.comm))
      return;

    Matrix eop0 = eop - eopModels;

    // Polar motion
    // ------------
    if(indexParameterPole)
    {
      Vector l(2*times.size());
      Matrix A(2*times.size(), 2*parametrizationPole->parameterCount());
      for(UInt idEpoch=0; idEpoch<times.size(); idEpoch++)
      {
        copy(eop0.slice(idEpoch, 0, 1, 2).trans(), l.row(2*idEpoch, 2));
        parametrizationPole->designMatrix(times.at(idEpoch), identityMatrix(2), A.row(2*idEpoch, 2));
      }
      axpy(3600e3*RAD2DEG, leastSquares(A, l), x0.row(normalEquationInfo.index(indexParameterPole), 2*parametrizationPole->parameterCount()));
    } // if(ParameterPole)

    // Earth rotation angle (UT1)
    // --------------------------
    if(indexParameterUT1)
    {
      // IGS convention seems to be that UT1/LOD parameters refer to UT1-UTC without being distorted by leap second jumps.
      // Therefore we remove the leap second jump by correcting for GPS2UTC-median(GPS2UTC).
      // Since a priori EOP parameters are undefined, all of this is weird and is only done to comply with IGS convention.
      Vector gps2Utc(times.size());
      for(UInt idEpoch=0; idEpoch<times.size(); idEpoch++)
        gps2Utc(idEpoch) = (times.at(idEpoch)-timeGPS2UTC(times.at(idEpoch))).seconds();
      gps2Utc -= median(gps2Utc);

      Vector l(times.size());
      Matrix A(times.size(), parametrizationUT1->parameterCount());
      for(UInt idEpoch=0; idEpoch<times.size(); idEpoch++)
      {
        l(idEpoch) = eop0(idEpoch, 3) - gps2Utc(idEpoch); // apply correction (see comment above)
        copy(parametrizationUT1->factors(times.at(idEpoch)).trans(), A.row(idEpoch));
      } // for(idEpoch)
      axpy(1e3, leastSquares(A, l), x0.row(normalEquationInfo.index(indexParameterUT1), parametrizationUT1->parameterCount()));
    } // if(ParameterUT1)

    // Nutation
    // --------
    if(indexParameterNutation)
    {
      Vector l(2*times.size());
      Matrix A(2*times.size(), 2*parametrizationNutation->parameterCount());
      for(UInt idEpoch=0; idEpoch<times.size(); idEpoch++)
      {
        copy(eop0.slice(idEpoch, 5, 1, 2).trans(), l.row(2*idEpoch, 2));
        parametrizationNutation->designMatrix(times.at(idEpoch), identityMatrix(2), A.row(2*idEpoch, 2));
      }
      axpy(3600e3*RAD2DEG, leastSquares(A, l), x0.row(normalEquationInfo.index(indexParameterNutation), 2*parametrizationNutation->parameterCount()));
    } // if(ParameterNutation)
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Bool GnssParametrizationEarthRotation::isDesignMatrix(const Gnss::NormalEquationInfo &/*normalEquationInfo*/, UInt /*idRecv*/, UInt /*idTrans*/, UInt /*idEpoch*/) const
{
  try
  {
//     if(!gnss.receiver.at(idRecv)->isEarthFixed())
//       return FALSE;
    return (indexParameterPole || indexParameterUT1 || indexParameterNutation);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssParametrizationEarthRotation::designMatrix(const Gnss::NormalEquationInfo &/*normalEquationInfo*/, const Gnss::ObservationEquation &eqn, Gnss::DesignMatrix &A) const
{
  try
  {
//     if(!eqn.receiver->isEarthFixed())
//       return;
    if(!(indexParameterPole || indexParameterUT1 || indexParameterNutation))
      return;

    Double xp, yp, sp, deltaUT, LOD, X, Y, S;
    earthOrientationParameter(eqn.timeRecv, xp, yp, sp, deltaUT, LOD, X, Y, S);
    const Double ERA = Planets::ERA(timeGPS2UTC(eqn.timeRecv) + seconds2time(deltaUT));

    // Polar motion rotation matrix
    const Matrix rotW = (rotaryY(Angle(xp)) * rotaryX(Angle(yp))).matrix();
    // Earth rotation angle rotation matrix
    const Matrix rotERA = rotaryZ(Angle(S-ERA-sp)).matrix();
    // Precession & nutation rotation matrix
    const Double a  = 0.5 + 0.125*(X*X+Y*Y);
    const Double X2 = X*X;
    const Double Y2 = Y*Y;
    const Matrix rotQ({{1.-a*X*X,   -a*X*Y,              X},
                       {  -a*X*Y, 1.-a*Y*Y,              Y},
                       {      -X,       -Y, 1.-a*(X*X+Y*Y)}});

    const Vector posEarth = rotaryMatrix(eqn.timeRecv).rotate(eqn.posRecv).vector();
    const Time timeTemporal = std::max(eqn.timeRecv, times.at(0));

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
      B = eqn.A.column(Gnss::ObservationEquation::idxPosRecv,3) * (rotQ * (rotERA * B));
      designMatrixTemporal(parametrizationPole, timeTemporal, B, indexParameterPole, A);
    }

    // Earth rotation angle (UT1)
    // --------------------------
    if(indexParameterUT1)
    {
      // partial derivatives of position(CRF) with respect to UT1
      Matrix partERA({{std::sin(S-ERA-sp), -std::cos(S-ERA-sp),  0},
                      {std::cos(S-ERA-sp),  std::sin(S-ERA-sp),  0},
                      {                 0,                   0,  0}});

      Matrix B = 1e-3 * 2*PI*1.00273781191135448/86400. * (eqn.A.column(Gnss::ObservationEquation::idxPosRecv,3) * (rotQ * (partERA * (rotW * posEarth))));
      designMatrixTemporal(parametrizationUT1, timeTemporal, B, indexParameterUT1, A);
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
      B = eqn.A.column(Gnss::ObservationEquation::idxPosRecv,3) * B;
      designMatrixTemporal(parametrizationNutation, timeTemporal, B, indexParameterNutation, A);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssParametrizationEarthRotation::designMatrixTemporal(ParametrizationTemporalPtr parametrization, const Time &time, const_MatrixSliceRef B,
                                             const Gnss::ParameterIndex &index, Gnss::DesignMatrix &A) const
{
  try
  {
    std::vector<UInt>   idx;
    std::vector<Double> factor;
    parametrization->factors(time, idx, factor);
    MatrixSlice Design(A.column(index));
    for(UInt i=0; i<factor.size(); i++)
      axpy(factor.at(i), B, Design.column(B.columns()*idx.at(i), B.columns()));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Double GnssParametrizationEarthRotation::updateParameter(const Gnss::NormalEquationInfo &normalEquationInfo, const_MatrixSliceRef x, const_MatrixSliceRef /*Wz*/, Bool printStatistics)
{
  try
  {
    Double      maxChangePole = 0;
    Double      maxChangeUT1  = 0;
    Double      maxChangeXY   = 0;
    std::string infoMaxChangePole, infoMaxChangeUT1, infoMaxChangeXY;

    // Polar motion
    // ------------
    if(indexParameterPole)
    {
      const Vector dx = DEG2RAD/3600e3 * x.row(normalEquationInfo.index(indexParameterPole), 2*parametrizationPole->parameterCount());
      std::vector<UInt>   index;
      std::vector<Double> factor;
      for(UInt idEpoch=0; idEpoch<times.size(); idEpoch++)
      {
        parametrizationPole->factors(times.at(idEpoch), index, factor);
        Vector p(2);
        for(UInt k=0; k<factor.size(); k++)
          axpy(factor.at(k), dx.row(2*index.at(k),2), p);

        // update xp,yp
        eop(idEpoch, 0) += p(0);
        eop(idEpoch, 1) += p(1);

        const Double dr = DEFAULT_R * norm(p);
        if(dr > maxChangePole)
        {
          maxChangePole = dr;
          infoMaxChangePole = "  earth rotation (pole):    posChange      = "+(1e3*dr)%"%6.1f mm"s;
        }
      } // for(idEpoch)
    } // if(ParameterPole)

    // Earth Rotation Angle (UT1)
    // --------------------------
    if(indexParameterUT1)
    {
      const Vector dx = 1e-3 * x.row(normalEquationInfo.index(indexParameterUT1), parametrizationUT1->parameterCount());
      std::vector<UInt>   index;
      std::vector<Double> factor;
      for(UInt idEpoch=0; idEpoch<times.size(); idEpoch++)
      {
        parametrizationUT1->factors(times.at(idEpoch), index, factor);
        Double p = 0;
        for(UInt k=0; k<factor.size(); k++)
          p += factor.at(k) * dx(index.at(k));

        // update UT1
        eop(idEpoch, 3) += p;

        const Double dr = DEFAULT_R * 2*PI*1.00273781191135448/86400 * std::fabs(p);
        if(dr > maxChangeUT1)
        {
          maxChangeUT1 = dr;
          infoMaxChangeUT1 = "  earth rotation (UT1):     posChange      = "+(1e3*dr)%"%6.1f mm"s;
        }
      } // for(idEpoch)
    } // if(ParameterUT1)

    // Nutation
    // --------
    if(indexParameterNutation)
    {
      const Vector dx = DEG2RAD/3600e3 * x.row(normalEquationInfo.index(indexParameterNutation), 2*parametrizationNutation->parameterCount());
      std::vector<UInt>   index;
      std::vector<Double> factor;
      for(UInt idEpoch=0; idEpoch<times.size(); idEpoch++)
      {
        parametrizationNutation->factors(times.at(idEpoch), index, factor);
        Vector p(2);
        for(UInt k=0; k<factor.size(); k++)
          axpy(factor.at(k), dx.row(2*index.at(k),2), p);

        // update X,Y
        eop(idEpoch, 5) += p(0);
        eop(idEpoch, 6) += p(1);

        const Double dr = DEFAULT_R * norm(p);
        if(dr > maxChangeXY)
        {
          maxChangeXY = dr;
          infoMaxChangeXY = "  earth rotation (XY):      posChange      = "+(1e3*dr)%"%6.1f mm"s;
        }
      } // for(idEpoch)
    } // if(ParameterNutation)

    Gnss::checkMaxChange(maxChangePole, infoMaxChangePole, printStatistics, normalEquationInfo.comm);
    Gnss::checkMaxChange(maxChangeUT1,  infoMaxChangeUT1,  printStatistics, normalEquationInfo.comm);
    Gnss::checkMaxChange(maxChangeXY,   infoMaxChangeXY,   printStatistics, normalEquationInfo.comm);
    return std::max(maxChangePole, std::max(maxChangeUT1, maxChangeXY));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void GnssParametrizationEarthRotation::writeResults(const Gnss::NormalEquationInfo &normalEquationInfo, const std::string &suffix)
{
  try
  {
    if(!Parallel::isMaster(normalEquationInfo.comm))
      return;

    if(!fileNameTimeSeriesEOP.empty())
    {
      VariableList fileNameVariableList;
      addTimeVariables(fileNameVariableList);
      evaluateTimeVariables(0, times.at(0), times.back(), fileNameVariableList);
      logStatus<<"write EOP to files <"<<fileNameTimeSeriesEOP(fileNameVariableList).appendBaseName(suffix)<<">"<<Log::endl;
      Matrix A(times.size(), 17);
      copy(eop,            A.column(1, 8));
      copy(eop-eopInitial, A.column(9, 8));
      for(UInt idEpoch = 0; idEpoch < times.size(); idEpoch++)
        A(idEpoch,4) += (times.at(idEpoch)-timeGPS2UTC(times.at(idEpoch))).seconds(); // UT1-GPS => UT1-UTC
      InstrumentFile::write(fileNameTimeSeriesEOP(fileNameVariableList).appendBaseName(suffix), Arc(times, A));
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
