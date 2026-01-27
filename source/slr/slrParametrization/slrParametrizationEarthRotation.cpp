/***********************************************/
/**
* @file slrParametrizationEarthRotation.cpp
*
* @brief EarthRotation.
* @see SlrParametrization
*
* @author Sebastian Strasser
* @author Torsten Mayer-Guerr
* @date 2022-04-28
*
*/
/***********************************************/

#include "base/import.h"
#include "config/config.h"
#include "base/planets.h"
#include "files/fileInstrument.h"
#include "classes/parametrizationTemporal/parametrizationTemporal.h"
#include "slr/slrParametrization/slrParametrizationEarthRotation.h"

/***** CLASS ***********************************/

/** @brief EarthRotation.
* @ingroup slrParametrizationGroup
* @see SlrParametrization */

/***********************************************/

SlrParametrizationEarthRotation::SlrParametrizationEarthRotation(Config &config)
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

void SlrParametrizationEarthRotation::init(Slr *slr, const std::vector<const SlrParametrizationGravityField*> &/*paramGravityField*/)
{
  try
  {
    this->slr = slr;
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

void SlrParametrizationEarthRotation::initParameter(SlrNormalEquationInfo &normalEquationInfo)
{
  try
  {
    indexParameterPole     = SlrParameterIndex();
    indexParameterUT1      = SlrParameterIndex();
    indexParameterNutation = SlrParameterIndex();

    if(!isEnabled(normalEquationInfo, name))
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

void SlrParametrizationEarthRotation::aprioriParameter(const SlrNormalEquationInfo &normalEquationInfo, MatrixSliceRef x0) const
{
  try
  {
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

void SlrParametrizationEarthRotation::designMatrix(const SlrNormalEquationInfo &/*normalEquationInfo*/, const SlrObservationEquation &eqn, SlrDesignMatrix &A) const
{
  try
  {
    if(!(indexParameterPole || indexParameterUT1 || indexParameterNutation))
      return;

    const Matrix eop = slr->polynomialEop.interpolate(eqn.timesStat, slr->eop);
    for(UInt idEpoch=0; idEpoch<eqn.timesStat.size(); idEpoch++)
    {
      const Time   time     = eqn.timesStat.at(idEpoch);
      const Double xp       = eop(idEpoch, 0);
      const Double yp       = eop(idEpoch, 1);
      const Double sp       = eop(idEpoch, 2);
      const Double deltaUT  = eop(idEpoch, 3) + (time-timeGPS2UTC(time)).seconds();
      const Double X        = eop(idEpoch, 5);
      const Double Y        = eop(idEpoch, 6);
      const Double S        = eop(idEpoch, 7);
      const Double ERA      = Planets::ERA(timeGPS2UTC(time) + seconds2time(deltaUT));
      const Matrix rotW     = (rotaryY(Angle(xp)) * rotaryX(Angle(yp))).matrix(); // Polar motion rotation matrix
      const Matrix rotERA   = rotaryZ(Angle(S-ERA-sp)).matrix();   // Earth rotation angle rotation matrix
      const Double X2       = X*X;
      const Double Y2       = Y*Y;
      const Double a        = 0.5 + 0.125*(X2+Y2);
      const Matrix rotQ({{1.-a*X2,   -a*X*Y,             X},    // Precession & nutation rotation matrix
                         {  -a*X*Y, 1.-a*Y2,             Y},
                         {      -X,      -Y, 1.-a*(X2+Y2)}});
      const Vector posEarth = slr->rotationCrf2Trf(time).rotate(eqn.posStat.at(idEpoch)).vector();

      // temporal parametrization
      auto designMatrixTemporal = [&](ParametrizationTemporalPtr parametrization, const_MatrixSliceRef B, const SlrParameterIndex &index)
      {
        std::vector<UInt>   idx;
        std::vector<Double> factor;
        parametrization->factors(time, idx, factor);
        MatrixSlice Design(A.column(index).row(eqn.index.at(idEpoch), eqn.count.at(idEpoch)));
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
        B = eqn.A.slice(eqn.index.at(idEpoch), SlrObservationEquation::idxPosStat, eqn.count.at(idEpoch), 3) * (rotQ * (rotERA * B));
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

        Matrix B = 1e-3 * 2*PI*1.00273781191135448/86400. * (eqn.A.slice(eqn.index.at(idEpoch), SlrObservationEquation::idxPosStat, eqn.count.at(idEpoch), 3) * (rotQ * (partERA * (rotW * posEarth))));
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
        B = eqn.A.slice(eqn.index.at(idEpoch), SlrObservationEquation::idxPosStat, eqn.count.at(idEpoch), 3) * B;
        designMatrixTemporal(parametrizationNutation, B, indexParameterNutation);
      }
    } // for(idEpoch)
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Double SlrParametrizationEarthRotation::updateParameter(const SlrNormalEquationInfo &normalEquationInfo, const_MatrixSliceRef x, const_MatrixSliceRef /*Wz*/)
{
  try
  {
    Double maxChange = 0;

    // Polar motion
    // ------------
    if(indexParameterPole)
    {
      Slr::InfoParameterChange info("mm");
      const Vector dx = DEG2RAD/3600e3 * x.row(normalEquationInfo.index(indexParameterPole), 2*parametrizationPole->parameterCount());
      xPole += dx;
      std::vector<UInt>   index;
      std::vector<Double> factor;
      for(UInt idEpoch=0; idEpoch<slr->times.size(); idEpoch++)
      {
        parametrizationPole->factors(slr->times.at(idEpoch), index, factor);
        Vector p(2);
        for(UInt k=0; k<factor.size(); k++)
          axpy(factor.at(k), dx.row(2*index.at(k),2), p);

        // update xp,yp
        slr->eop(idEpoch, 0) += p(0);
        slr->eop(idEpoch, 1) += p(1);
        if(info.update(1e3*DEFAULT_R*norm(p)))
          info.info = "earth rotation (pole)";
      } // for(idEpoch)
      info.print(1e-3, maxChange);
    } // if(ParameterPole)

    // Earth Rotation Angle (UT1)
    // --------------------------
    if(indexParameterUT1)
    {
      Slr::InfoParameterChange info("mm");
      const Vector dx = 1e-3 * x.row(normalEquationInfo.index(indexParameterUT1), parametrizationUT1->parameterCount());
      xUT1 += dx;
      std::vector<UInt>   index;
      std::vector<Double> factor;
      for(UInt idEpoch=0; idEpoch<slr->times.size(); idEpoch++)
      {
        parametrizationUT1->factors(slr->times.at(idEpoch), index, factor);
        Double p = 0;
        for(UInt k=0; k<factor.size(); k++)
          p += factor.at(k) * dx(index.at(k));

        slr->eop(idEpoch, 3) += p;
        if(info.update(1e3*DEFAULT_R*2*PI*1.00273781191135448/86400*p))
          info.info = "earth rotation (UT1)";
      } // for(idEpoch)
      info.print(1e-3, maxChange);
    } // if(ParameterUT1)

    // Nutation
    // --------
    if(indexParameterNutation)
    {
      Slr::InfoParameterChange info("mm");
      const Vector dx = DEG2RAD/3600e3 * x.row(normalEquationInfo.index(indexParameterNutation), 2*parametrizationNutation->parameterCount());
      xNutation += dx;
      std::vector<UInt>   index;
      std::vector<Double> factor;
      for(UInt idEpoch=0; idEpoch<slr->times.size(); idEpoch++)
      {
        parametrizationNutation->factors(slr->times.at(idEpoch), index, factor);
        Vector p(2);
        for(UInt k=0; k<factor.size(); k++)
          axpy(factor.at(k), dx.row(2*index.at(k),2), p);

        // update X,Y
        slr->eop(idEpoch, 5) += p(0);
        slr->eop(idEpoch, 6) += p(1);
        if(info.update(1e3*DEFAULT_R*norm(p)))
          info.info = "earth rotation (XY)";
      } // for(idEpoch)
      info.print(1e-3, maxChange);
    } // if(ParameterNutation)

    return maxChange;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void SlrParametrizationEarthRotation::writeResults(const SlrNormalEquationInfo &normalEquationInfo, const std::string &suffix) const
{
  try
  {
    if(!isEnabled(normalEquationInfo, name))
      return;

    if(!fileNameEOP.empty())
    {
      logStatus<<"write EOP to files <"<<fileNameEOP.appendBaseName(suffix)<<">"<<Log::endl;
      Matrix A(slr->eop.rows(), 1+slr->eop.columns()); // prepend time column
      copy(slr->eop, A.column(1, slr->eop.columns()));
      InstrumentFile::write(fileNameEOP.appendBaseName(suffix), Arc(slr->times, A));
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
