/***********************************************/
/**
* @file preprocessingGradiometer.cpp
*
* @brief Estimate covariance function / arc weights.
*
* @author Torsten Mayer-Guerr
* @date 2011-05-11
*
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program estimates empirical covariance functions of the gradiometer instrument noise and determine arc wise variances to
downweight arcs with outliers. This program works similar to \program{PreprocessingPod}, see there for details.
Here only the settings explained, which are different.

...
)";

/***********************************************/

#include "programs/program.h"
#include "files/fileMatrix.h"
#include "files/fileInstrument.h"
#include "classes/earthRotation/earthRotation.h"
#include "classes/ephemerides/ephemerides.h"
#include "classes/gravityfield/gravityfield.h"
#include "classes/parametrizationTemporal/parametrizationTemporal.h"
#include "misc/observation/observationMisc.h"
#include "misc/varianceComponentEstimation.h"

/***** CLASS ***********************************/

/** @brief Estimate covariance function / arc weights.
* @ingroup programsGroup */
class PreprocessingGradiometer
{
  SggRightSidePtr            rhs; // right hand sides
  InstrumentFile             orbitFile;
  InstrumentFile             starCameraFile;
  EarthRotationPtr           earthRotation;
  EphemeridesPtr             ephemerides;
  ParametrizationTemporalPtr sggBias;
  Vector                     sigmas, sigmasNew;

  // covariance
  // ----------
  Double sampling;
  UInt   covLength;
  Matrix covFunc;
  Matrix CovCholesky;
  Matrix Psd, ePe, redundancy; // one row for each frequency, one column for each component
  Matrix CosTransform;

  GradiometerArc computeArc(UInt arc);

public:
  void run(Config &config);
};

GROOPS_REGISTER_PROGRAM(PreprocessingGradiometer, PARALLEL, "Estimate covariance function / arc weights.", Preprocessing)

/***********************************************/

void PreprocessingGradiometer::run(Config &config)
{
  try
  {
    FileName outCovName, outResidualsName, outSigmaName;
    FileName orbitName, starCameraName, covName;
    UInt     iterCount;

    readConfig(config, "outputfileCovarianceFunction", outCovName,       Config::OPTIONAL, "", "");
    readConfig(config, "outputfileSigmasPerArc",       outSigmaName,     Config::OPTIONAL, "", "accuracies of each arc");
    readConfig(config, "outputfileSggResiduals",       outResidualsName, Config::OPTIONAL, "", "");
    readConfig(config, "rightHandSide",                rhs,              Config::MUSTSET,  "", "input for the observation vector");
    readConfig(config, "inputfileOrbit",               orbitName,        Config::MUSTSET,  "", "");
    readConfig(config, "inputfileStarCamera",          starCameraName,   Config::MUSTSET,  "", "");
    readConfig(config, "earthRotation",                earthRotation,    Config::MUSTSET,  "", "");
    readConfig(config, "ephemerides",                  ephemerides,      Config::OPTIONAL, "jpl", "");
    readConfig(config, "parametrizationBias",          sggBias,          Config::DEFAULT,  "",    "per arc");
    if(readConfigSequence(config, "covarianceSgg", Config::MUSTSET, "", ""))
    {
      readConfig(config, "inputfileCovarianceFunction",  covName,   Config::OPTIONAL, "",  "approximate covariances in time");
      readConfig(config, "covarianceLength",             covLength, Config::MUSTSET,  "",  "counts observation epochs");
      readConfig(config, "sampling",                     sampling,  Config::MUSTSET,  "5", "[seconds] sampling of the covariance function");
      endSequence(config);
    }
    readConfig(config, "iterationCount", iterCount, Config::DEFAULT,  "3", "for the estimation of calibration parameter and error PSD");
    if(isCreateSchema(config)) return;

    // ======================================================

    // test instrument files
    // ---------------------
    orbitFile.open(orbitName);
    starCameraFile.open(starCameraName);
    InstrumentFile::checkArcCount({orbitFile, starCameraFile, *rhs->gradiometerFile});
    for(UInt k=0; k<rhs->referenceFile.size(); k++)
      InstrumentFile::checkArcCount({orbitFile, *rhs->referenceFile.at(k)});
    const UInt arcCount = orbitFile.arcCount();

    // =======================

    // init covariance function
    // ------------------------
    CosTransform = Vce::cosTransform(covLength);
    covFunc      = Vce::readCovarianceFunction(covName, covLength, 6, sampling);
    Psd          = CosTransform * covFunc.column(1,6);

    // init arc sigmas
    // ---------------
    sigmas = Vector(arcCount);
    for(UInt arcNo=0; arcNo<arcCount; arcNo++)
      sigmas(arcNo) = 1.0;

    // =============================================

    // Iteration
    // ---------
    for(UInt iter=0; iter<iterCount; iter++)
    {
      // new covariance matricies
      // ------------------------
      CovCholesky = Matrix(6*covFunc.rows(), Matrix::SYMMETRIC, Matrix::UPPER);
      for(UInt sggNo=0; sggNo<6; sggNo++)
        for(UInt i=0; i<covFunc.rows(); i++)
          for(UInt k=i; k<covFunc.rows(); k++)
            CovCholesky(6*i+sggNo, 6*k+sggNo) = covFunc(k-i, 1+sggNo);
      cholesky(CovCholesky);

      logStatus<<"set up observation equations"<<Log::endl;
      ePe        = Matrix(covLength, 6); // xx, xy, xz, yy, yz, zz
      redundancy = Matrix(covLength, 6);
      sigmasNew  = Vector(arcCount);
      std::vector<Arc> arcList(arcCount);
      Parallel::forEach(arcList, [this](UInt arcNo) {return computeArc(arcNo);});
      Parallel::reduceSum(ePe);
      Parallel::reduceSum(redundancy);
      Parallel::reduceSum(sigmasNew);

      // sigmas per arc
      // --------------
      // median sigma per arc;
      if(Parallel::isMaster())
      {
        sigmas = sigmasNew;
        Double sigma0 = Vce::meanSigma(sigmas);
        sigmas *= 1./sigma0;

        logInfo<<"  sigma per arc (median): "<<sigma0<<Log::endl;

        if(!outSigmaName.empty())
        {
          logStatus<<"write arc sigma file <"<<outSigmaName<<">"<<Log::endl;
          writeFileMatrix(outSigmaName, sigmas);
        }
      }
      Parallel::broadCast(sigmas);

      // estimate new PSD
      // ----------------
      if(Parallel::isMaster())
      {
        Double maxFactor = 0;
        Vce::estimatePsd(ePe, redundancy, Psd, maxFactor);
        logInfo<<"  max. PSD adjustment factor: "<<maxFactor<<Log::endl;
      } // if(Parallel::isMaster())
      Parallel::broadCast(Psd);
      copy(CosTransform * Psd, covFunc.column(1,Psd.columns())); // compute new covariance function

      if(Parallel::isMaster() && !outCovName.empty())
      {
        logStatus<<"write covariance function file <"<<outCovName<<">"<<Log::endl;
        writeFileMatrix(outCovName, covFunc);
      }

      // gradiometer residuals
      // ---------------------
      if(Parallel::isMaster() && !outResidualsName.empty())
      {
        logStatus<<"write gradiometer residuals file <"<<outResidualsName<<">"<<Log::endl;
        InstrumentFile::write(outResidualsName, arcList);
      }
    } // for(iter)
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

GradiometerArc PreprocessingGradiometer::computeArc(UInt arcNo)
{
  try
  {
    OrbitArc      orbit      = orbitFile.readArc(arcNo);
    StarCameraArc starCamera = starCameraFile.readArc(arcNo);
    UInt          epochCount = orbit.size();

    // =============================================

    // reduced observations
    // ---------------------
    GradiometerArc sggArc = rhs->gradiometerFile->readArc(arcNo);
    Arc::checkSynchronized({sggArc, orbit, starCamera});
    std::vector<GradiometerArc> reference(rhs->referenceFile.size());
    for(UInt k=0; k<rhs->referenceFile.size(); k++)
    {
      reference.at(k) = rhs->referenceFile.at(k)->readArc(arcNo);
      Arc::checkSynchronized({orbit, reference.at(k)});
    }

    for(UInt i=0; i<epochCount; i++)
    {
      const Time     time     = orbit.at(i).time;
      const Rotary3d rotEarth = earthRotation->rotaryMatrix(time);
      const Vector3d posEarth = rotEarth.rotate(orbit.at(i).position);

      // referencefield + tides
      const Tensor3d tns = rhs->referencefield->gravityGradient(time, posEarth)
                          + rhs->tides->gradient(time, posEarth, rotEarth, earthRotation, ephemerides);
      // observed minus computed
      sggArc.at(i).gravityGradient -= starCamera.at(i).rotary.inverseRotate(rotEarth.inverseRotate(tns));
      // gradients from files
      for(UInt k=0; k<reference.size(); k++)
        sggArc.at(i).gravityGradient -= reference.at(k).at(i).gravityGradient;
    }

    // =============================================

    Vector l(6*epochCount);
    for(UInt i=0; i<epochCount; i++)
    {
      l(6*i+0) = sggArc.at(i).gravityGradient.xx();
      l(6*i+1) = sggArc.at(i).gravityGradient.xy();
      l(6*i+2) = sggArc.at(i).gravityGradient.xz();
      l(6*i+3) = sggArc.at(i).gravityGradient.yy();
      l(6*i+4) = sggArc.at(i).gravityGradient.yz();
      l(6*i+5) = sggArc.at(i).gravityGradient.zz();
    }

    Matrix A, B;

    // gradiometer bias for each component
    // -----------------------------------
    const std::vector<Time> times = orbit.times();
    sggBias->setInterval(times.front(), times.back()+medianSampling(times), TRUE);
    if(sggBias->parameterCount())
    {
      B = Matrix(6*epochCount, 6*sggBias->parameterCount());
      const Matrix I = identityMatrix(6);
      for(UInt i=0; i<epochCount; i++)
        sggBias->designMatrix(times.at(i), I, B.row(6*i, 6));
    }

    // =============================================

    // Decorrelation
    // -------------
    Matrix WA = A;
    Matrix WB = B;
    Matrix Wl = l;

    Matrix W  = sigmas(arcNo) * CovCholesky.slice(0,0,6*epochCount,6*epochCount);
    if(WA.size()) triangularSolve(1., W.trans(), WA);
    if(WB.size()) triangularSolve(1., W.trans(), WB);
    if(Wl.size()) triangularSolve(1., W.trans(), Wl);

    // =============================================

    // solution vector & decorrelated residuals
    // ----------------------------------------
    Vector tau = QR_decomposition(WB);
    QTransMult(WB, tau, Wl); // transform observations: l:= Q'l
    if(WA.size())
      QTransMult(WB, tau, WA); // transform design matrix A:=Q'A
    Matrix We = Wl; // transformed residuals
//     if(WA.size())
//       matMult(-1., WA, x, We);
    // solution of arc dependent parameters
    Vector y = We.row(0, WB.columns());
    triangularSolve(1., WB.row(0, tau.rows()), y);
    // residuals
    We.row(0, WB.columns()).setNull(); // remove WB*x
    QMult(WB, tau, We); // back transformation
    generateQ(WB, tau);

    // ============================================

    // Variance component estimation
    // -----------------------------
    std::vector<UInt> index(epochCount);
    for(UInt i=0; i<index.size(); i++)
      index.at(i) = i;
    Double ePeSum=0, redundancySum=0;
    {
      Matrix R;
      Vector WWe;
      Vce::redundancy(W, We, WA, WB, R, WWe);
      Vce::psd(R, WWe, index, sigmas(arcNo), CosTransform, Psd, ePe, redundancy, ePeSum, redundancySum);
    }
    sigmasNew(arcNo) = std::sqrt(ePeSum/redundancySum) * sigmas(arcNo);  // compute new sigma (for this arc)

    // =============================================

    // compute residuals
    // -----------------
    const Vector e = l - B*y;
//     if(A.size())
//       e -= A*x;

    for(UInt i=0; i<epochCount; i++)
    {
      sggArc.at(i).gravityGradient.xx() = e(6*i+0);
      sggArc.at(i).gravityGradient.xy() = e(6*i+1);
      sggArc.at(i).gravityGradient.xz() = e(6*i+2);
      sggArc.at(i).gravityGradient.yy() = e(6*i+3);
      sggArc.at(i).gravityGradient.yz() = e(6*i+4);
      sggArc.at(i).gravityGradient.zz() = e(6*i+5);
    }

    return sggArc;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
