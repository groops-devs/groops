/***********************************************/
/**
* @file observationGradiometer.cpp
*
* @brief GOCE gradiometer observations.
*
* @author Torsten Mayer-Guerr
* @date 2011-05-14
*
*/
/***********************************************/

#include "base/import.h"
#include "files/fileMatrix.h"
#include "files/fileInstrument.h"
#include "classes/earthRotation/earthRotation.h"
#include "classes/tides/tides.h"
#include "classes/gravityfield/gravityfield.h"
#include "classes/parametrizationGravity/parametrizationGravity.h"
#include "classes/parametrizationTemporal/parametrizationTemporal.h"
#include "misc/observation/observationMisc.h"
#include "classes/observation/observation.h"
#include "classes/observation/observationGradiometer.h"

/***********************************************/

ObservationGradiometer::ObservationGradiometer(Config &config)
{
  try
  {
    FileName  orbitName, starCameraName, covarianceName;
    FileName  sigmaName, covName;

    renameDeprecatedConfig(config, "representation", "parametrizationGravity", date2time(2020, 6, 3));

    readConfig(config, "rightHandSide",          rhs,             Config::MUSTSET,  "",    "input for the observation vector");
    readConfig(config, "inputfileOrbit",         orbitName,       Config::MUSTSET,  "",    "");
    readConfig(config, "inputfileStarCamera",    starCameraName,  Config::MUSTSET,  "",    "");
    readConfig(config, "earthRotation",          earthRotation,   Config::MUSTSET,  "",    "");
    readConfig(config, "ephemerides",            ephemerides,     Config::OPTIONAL, "jpl", "");
    readConfig(config, "parametrizationGravity", parametrization, Config::MUSTSET,  "",    "");
    readConfig(config, "parametrizationBias",    sggBias,         Config::DEFAULT,  "",    "per arc");
    readConfig(config, "useXX",                  useXX,           Config::DEFAULT,  "1",   "");
    readConfig(config, "useYY",                  useYY,           Config::DEFAULT,  "1",   "");
    readConfig(config, "useZZ",                  useZZ,           Config::DEFAULT,  "1",   "");
    readConfig(config, "useXY",                  useXY,           Config::DEFAULT,  "0",   "");
    readConfig(config, "useXZ",                  useXZ,           Config::DEFAULT,  "1",   "");
    readConfig(config, "useYZ",                  useYZ,           Config::DEFAULT,  "0",   "");
    if(readConfigSequence(config, "covarianceSgg", Config::MUSTSET, "", ""))
    {
      readConfig(config, "sigma",                       sigma,     Config::DEFAULT,  "1", "general variance factor");
      readConfig(config, "inputfileSigmasPerArc",       sigmaName, Config::OPTIONAL, "",  "different accuaries for each arc (multplicated with sigma)");
      readConfig(config, "inputfileCovarianceFunction", covName,   Config::OPTIONAL, "",  "covariance function in time");
      endSequence(config);
    }
    if(isCreateSchema(config)) return;

    orbitFile.open(orbitName);
    starCameraFile.open(starCameraName);

    InstrumentFile::checkArcCount({orbitFile, starCameraFile});
    for(UInt rhsNo=0; rhsNo<rhs.size(); rhsNo++)
    {
      InstrumentFile::checkArcCount({orbitFile, *rhs.at(rhsNo)->gradiometerFile});
      for(UInt k=0; k<rhs.at(rhsNo)->referenceFile.size(); k++)
        InstrumentFile::checkArcCount({orbitFile, *rhs.at(rhsNo)->referenceFile.at(k)});
    }

    componentCount = useXX + useXY + useXZ + useYY + useYZ + useZZ;

    if(!sigmaName.empty())
      readFileMatrix(sigmaName, sigmaArc);

    // covariance matrix from covariance function
    // ------------------------------------------
    if(!covName.empty())
    {
      Matrix cov;
      readFileMatrix(covName, cov);

      CovCholesky = Matrix(componentCount*cov.rows(), Matrix::SYMMETRIC, Matrix::UPPER);

      for(UInt i=0; i<cov.rows(); i++)
        for(UInt k=i; k<cov.rows(); k++)
        {
          UInt idx = 0;
          if(useXX) {CovCholesky(componentCount*i+idx, componentCount*k+idx) = cov(k-i, 1+0); idx++;}
          if(useXY) {CovCholesky(componentCount*i+idx, componentCount*k+idx) = cov(k-i, 1+1); idx++;}
          if(useXZ) {CovCholesky(componentCount*i+idx, componentCount*k+idx) = cov(k-i, 1+2); idx++;}
          if(useYY) {CovCholesky(componentCount*i+idx, componentCount*k+idx) = cov(k-i, 1+3); idx++;}
          if(useYZ) {CovCholesky(componentCount*i+idx, componentCount*k+idx) = cov(k-i, 1+4); idx++;}
          if(useZZ) {CovCholesky(componentCount*i+idx, componentCount*k+idx) = cov(k-i, 1+5); idx++;}
        }

      cholesky(CovCholesky);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/************************************************************************/

void ObservationGradiometer::observation(UInt arcNo, Matrix &l, Matrix &A, Matrix &B)
{
  try
  {
    OrbitArc      orbit      = orbitFile.readArc(arcNo);
    StarCameraArc starCamera = starCameraFile.readArc(arcNo);
    const UInt    epochCount = orbit.size();
    const UInt    rhsCount   = rhs.size();
    Arc::checkSynchronized({orbit, starCamera});

    // earth rotation
    // --------------
    std::vector<Rotary3d> rotEarth(epochCount);
    for(UInt i=0; i<epochCount; i++)
      rotEarth.at(i) = earthRotation->rotaryMatrix(orbit.at(i).time);

    // reduced observations
    // ---------------------
    l = Matrix(componentCount*epochCount, rhsCount);

    for(UInt rhsNo=0; rhsNo<rhsCount; rhsNo++)
    {
      GradiometerArc gradiometer = rhs.at(rhsNo)->gradiometerFile->readArc(arcNo);
      Arc::checkSynchronized({orbit, gradiometer});

      std::vector<GradiometerArc> reference(rhs.at(rhsNo)->referenceFile.size());
      for(UInt k=0; k<rhs.at(rhsNo)->referenceFile.size(); k++)
      {
        reference.at(k) = rhs.at(rhsNo)->referenceFile.at(k)->readArc(arcNo);
        Arc::checkSynchronized({orbit, reference.at(k)});
      }

      for(UInt i=0; i<epochCount; i++)
      {
        const Time     time     = orbit.at(i).time;
        const Vector3d posEarth = rotEarth.at(i).rotate(orbit.at(i).position);
        // referencefield + tides
        const Tensor3d tns = rhs.at(rhsNo)->referencefield->gravityGradient(time, posEarth)
                           + rhs.at(rhsNo)->tides->gradient(time, posEarth, rotEarth.at(i), earthRotation, ephemerides);
        // observed minus computed
        Tensor3d gravityGradient = gradiometer.at(i).gravityGradient - starCamera.at(i).rotary.inverseRotate(rotEarth.at(i).inverseRotate(tns));
        // gradients from files
        for(UInt k=0; k<reference.size(); k++)
          gravityGradient -= reference.at(k).at(i).gravityGradient;

        UInt idx = 0;
        if(useXX) l(componentCount*i+idx++, rhsNo) = gravityGradient.xx();
        if(useXY) l(componentCount*i+idx++, rhsNo) = gravityGradient.xy();
        if(useXZ) l(componentCount*i+idx++, rhsNo) = gravityGradient.xz();
        if(useYY) l(componentCount*i+idx++, rhsNo) = gravityGradient.yy();
        if(useYZ) l(componentCount*i+idx++, rhsNo) = gravityGradient.yz();
        if(useZZ) l(componentCount*i+idx++, rhsNo) = gravityGradient.zz();
      }
    }

    // rotary matrix from TRF to satellite system
    // ------------------------------------------
    Matrix rotGRF(componentCount*epochCount, 5*epochCount);
    for(UInt i=0; i<epochCount; i++)
    {
      Matrix rot = inverse(rotEarth.at(i) * starCamera.at(i).rotary).matrix();
      MatrixSlice R(rotGRF.slice(componentCount*i, 5*i, componentCount,5));

      // One row of the rotary matrix for one gradiometer component (e.g. Txy: i=0, k=1)
      auto rotationLine = [&](UInt i, UInt k, UInt row)
      {
        R(row, 0) = rot(i,0)*rot(k,0) - rot(i,2)*rot(k,2);
        R(row, 1) = rot(i,0)*rot(k,1) + rot(i,1)*rot(k,0);
        R(row, 2) = rot(i,0)*rot(k,2) + rot(i,2)*rot(k,0);
        R(row, 3) = rot(i,1)*rot(k,1) - rot(i,2)*rot(k,2);
        R(row, 4) = rot(i,1)*rot(k,2) + rot(i,2)*rot(k,1);
      };

      UInt idx = 0;
      if(useXX) rotationLine(0, 0, idx++);
      if(useXY) rotationLine(0, 1, idx++);
      if(useXZ) rotationLine(0, 2, idx++);
      if(useYY) rotationLine(1, 1, idx++);
      if(useYZ) rotationLine(1, 2, idx++);
      if(useZZ) rotationLine(2, 2, idx++);
    }

    // gradiometer bias for each component
    // -----------------------------------
    B = Matrix();
    const std::vector<Time> times = orbit.times();
    sggBias->setInterval(times.front(), times.back()+medianSampling(times), TRUE);
    if(sggBias->parameterCount())
    {
      B = Matrix(componentCount*epochCount, componentCount*sggBias->parameterCount());
      const Matrix I = identityMatrix(componentCount);
      for(UInt i=0; i<epochCount; i++)
        sggBias->designMatrix(times.at(i), I, B.row(componentCount*i, componentCount));
    }

    // apply sigmas
    // ------------
    Double factor = 1./sigma;
    if(sigmaArc.size())
      factor *= 1./sigmaArc(arcNo);
    if(factor!=1.)
    {
      rotGRF *= factor;
      l      *= factor;
      if(B.size())
        B *= factor;
    }

    // decorrelation
    // -------------
    if(CovCholesky.size())
    {
      if(CovCholesky.rows()<l.rows())
        throw(Exception("covariance matrix to small"));
      const_MatrixSlice W(CovCholesky.slice(0,0,l.rows(),l.rows()).trans());
      triangularSolve(1., W, rotGRF);
      triangularSolve(1., W, l);
      if(B.size())
        triangularSolve(1.,W, B);
    }

    // Design matrix A
    // ---------------
    A = Matrix(componentCount*epochCount, parameterCount());
    Matrix tns(6, parametrization->parameterCount());
    for(UInt i=0; i<epochCount; i++)
    {
      parametrization->gravityGradient(orbit.at(i).time, rotEarth.at(i).rotate(orbit.at(i).position), tns);
      matMult(1., rotGRF.column(5*i,5), tns.row(0,5), A.column(0,tns.columns()));
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
