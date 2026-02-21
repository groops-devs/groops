/***********************************************/
/**
* @file preprocessingVariationalEquation.cpp
*
* @brief Integrate Variational Equations.
*
* @author Torsten Mayer-Guerr
* @date 2012-05-30
*
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program integrates an orbit dynamically using the given forces and set up the state transition matrix
for each time step. These are the prerequisites for a least squares adjustment (e.g. gravity field determination) using
the variational equation approach. The variational equations are computed arc-wise as defined by \configFile{inputfileOrbit}{instrument}.
This means for each arc new initial state parameters are set up.

In a first step the \configClass{forces}{forcesType} acting on the satellite are evaluated at the apriori positions given
by \configFile{inputfileOrbit}{instrument}. Non-conservative forces like solar radiation pressure need the orientation of the
satellite (\configFile{inputfileStarCamera}{instrument}) and additionally, a satellite macro model (\config{satelliteModel})
with the surface properties. Furthermore \configFile{inputfileAccelerometer}{instrument} observations are also considered.

In a second step the accelerations are integrated twice to a dynamic orbit using a moving polynomial with the degree
\config{integrationDegree}. The orbit is corrected to be self-consistent. This means the forces should be evaluated
at the new integrated positions instead of the apriori ones. This correction is computed in a linear approximation
using the gradient of the forces with respect to the positions (\config{gradientfield}). As this term is small generally
only the largest force components have to be considered. A low degree spherical harmonic expansion of the static gravity
field (about up to degree 5) is sufficient in almost all cases. In this step also the state transition matrix (the partial
derivatives of the current state, position and velocity) with respect to the initial state is computed.
The integrated orbit together with the state transitions are stored in \configFile{outputfileVariational}{variationalEquation},
the integrated orbit only in \configFile{outputfileOrbit}{instrument}.

To improve the numerical stability a reference ellipse can be reduced beforehand using Enke's method (\config{useEnke}).
Mathematically the result is the same, but as the large central term is removed before and restored
afterwards more digits are available for the computation.

The integrated orbit should be fitted to observations afterwards by the programs
\program{PreprocessingVariationalEquationOrbitFit} and/or \program{PreprocessingVariationalEquationSstFit}.
They apply a least squares adjustment by estimating some satellite parameters (e.g. an accelerometer bias).
If the fitted orbit is too far away from the original \configFile{inputfileOrbit}{instrument} the linearization may not be
accurate enough. In this case \program{PreprocessingVariationalEquation} should be run again with the fitted orbit
as \configFile{inputfileOrbit}{instrument} and introducing the \config{estimatedParameters} as additional forces.
)";

/***********************************************/

#include "programs/program.h"
#include "base/equinoctial.h"
#include "files/fileMatrix.h"
#include "files/fileInstrument.h"
#include "files/fileSatelliteModel.h"
#include "files/fileVariationalEquation.h"
#include "classes/earthRotation/earthRotation.h"
#include "classes/ephemerides/ephemerides.h"
#include "classes/forces/forces.h"
#include "classes/gravityfield/gravityfield.h"
#include "classes/parametrizationAcceleration/parametrizationAcceleration.h"

/***** CLASS ***********************************/

/** @brief Integrate Variational Equations.
* @ingroup programsGroup */
class PreprocessingVariationalEquation
{
  SatelliteModelPtr              satellite;
  InstrumentFile                 orbitFile;
  InstrumentFile                 starCameraFile;
  InstrumentFile                 accelerometerFile;
  ForcesPtr                      forces;
  EarthRotationPtr               earthRotation;
  EphemeridesPtr                 ephemerides;
  GravityfieldPtr                gradientfield;
  ParametrizationAccelerationPtr parameterAcceleration;
  Matrix                         parameter;
  Double                         GM;
  Double                         maxPosDiff;
  UInt                           integrationDegree;
  std::vector<Vector>            coeff_g, coeff_tg;

  VariationalEquationArc computeArc(UInt arcNo);
  Matrix integrate2Position(Double deltaT, const_MatrixSliceRef g) const;
  Matrix integrate2Velocity(Double deltaT, const_MatrixSliceRef g) const;
  Matrix solve(Double deltaT, const std::vector<Tensor3d> &tensor, const_MatrixSliceRef l) const;
  Matrix refine(Double deltaT, const std::vector<Tensor3d> &tensor, const_MatrixSliceRef x) const;
  Matrix approxInverse(Double deltaT, const std::vector<Tensor3d> &tensor, const_MatrixSliceRef x) const;

public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(PreprocessingVariationalEquation, PARALLEL, "Integrate Variational Equations.", Preprocessing, VariationalEquation)

/***********************************************/
/***********************************************/

void PreprocessingVariationalEquation::run(Config &config, Parallel::CommunicatorPtr comm)
{
  try
  {
    FileName fileNameOutVariational, fileNameOutOrbit;
    FileName fileNameSatellite;
    FileName orbitName, starCameraName, accName;
    FileName fileNameSolution;
    GM = 0;

    renameDeprecatedConfig(config, "satelliteModel", "inputfileSatelliteModel", date2time(2020, 8, 19));

    readConfig(config, "outputfileVariational",       fileNameOutVariational, Config::MUSTSET,   "", "approximate position and integrated state matrix");
    readConfig(config, "outputfileOrbit",             fileNameOutOrbit,       Config::OPTIONAL,  "", "integrated orbit");
    readConfig(config, "inputfileSatelliteModel",     fileNameSatellite,      Config::OPTIONAL,  "{groopsDataDir}/satelliteModel/", "satellite macro model");
    readConfig(config, "inputfileOrbit",              orbitName,              Config::MUSTSET,   "", "approximate position, used to evaluate the force");
    readConfig(config, "inputfileStarCamera",         starCameraName,         Config::MUSTSET,   "", "rotation from body frame to CRF");
    readConfig(config, "inputfileAccelerometer",      accName,                Config::OPTIONAL,  "", "non-gravitational forces in satellite reference frame");
    readConfig(config, "forces",                      forces,                 Config::MUSTSET,   "", "");
    if(readConfigSequence(config, "estimatedParameters", Config::OPTIONAL, "", "satellite parameters e.g. from orbit fit"))
    {
      renameDeprecatedConfig(config, "parameter", "parametrizationAcceleration", date2time(2020, 6, 3));

      readConfig(config, "parametrizationAcceleration", parameterAcceleration, Config::MUSTSET,  "", "orbit force parameters");
      readConfig(config, "inputfileParameter",          fileNameSolution,      Config::MUSTSET,  "", "estimated orbit force parameters");
      endSequence(config);
    }
    readConfig(config, "earthRotation",     earthRotation,     Config::MUSTSET,  "",    "");
    readConfig(config, "ephemerides",       ephemerides,       Config::OPTIONAL, "jpl", "");
    readConfig(config, "gradientfield",     gradientfield,     Config::MUSTSET,  "",    "low order field to estimate the change of the gravity by position adjustement");
    readConfig(config, "integrationDegree", integrationDegree, Config::DEFAULT,  "7",   "integration of forces by polynomial approximation of degree n");
    if(readConfigSequence(config, "useEnke", Config::OPTIONAL, "1", "integrate differential forces to an elliptical reference trajectory"))
    {
      readConfig(config, "GM", GM, Config::DEFAULT,  STRING_DEFAULT_GM, "geocentric gravitational constant used for elliptical reference orbit");
      endSequence(config);
    }
    if(isCreateSchema(config)) return;

    if(integrationDegree%2 == 0)
      throw(Exception("polynomial degree for integration must be odd."));

    if(!fileNameSatellite.empty()) readFileSatelliteModel(fileNameSatellite, satellite);

    // open and test instrument files
    // ------------------------------
    orbitFile.open(orbitName);
    starCameraFile.open(starCameraName);
    accelerometerFile.open(accName);
    InstrumentFile::checkArcCount({orbitFile, starCameraFile, accelerometerFile});

    // init calibration parameters
    if(parameterAcceleration)
    {
      readFileMatrix(fileNameSolution, parameter);
      if(parameterAcceleration->parameterCount() > parameter.rows())
        throw(Exception("parameter and inputfileParameter must agree"));
    }

    // init integration
    // ----------------
    coeff_g.resize (integrationDegree);
    coeff_tg.resize(integrationDegree);
    for(UInt idx=0; idx<coeff_g.size(); idx++)
    {
      // polynomial matrix
      Matrix W(integrationDegree+1, integrationDegree+1);
      for(UInt i=0; i<W.columns(); i++)
      {
        W(i,0) = 1.;
        for(UInt n=1; n<W.rows(); n++)
          W(i,n) = (Double(i)-idx) * W(i,n-1);
      }
      inverse(W);

      coeff_g.at(idx)  = Vector(W.rows());
      coeff_tg.at(idx) = Vector(W.rows());
      for(UInt i=0; i<W.rows(); i++)
        for(UInt n=0; n<W.columns(); n++)
        {
          coeff_g.at(idx)(i)  += 1./(n+1.) * W(n,i);
          coeff_tg.at(idx)(i) += 1./(n+2.) * W(n,i);
        }
    }

    // compute arcs
    // ------------
    logStatus<<"integrate arcs"<<Log::endl;
    maxPosDiff = 0;
    std::vector<VariationalEquationArc> arcs(orbitFile.arcCount());
    Parallel::forEach(arcs, [this](UInt arcNo) {return computeArc(arcNo);}, comm);
    Parallel::reduceMax(maxPosDiff, 0, comm);
    logInfo<<"  max diff (orbit-reference orbit) "<<maxPosDiff<<" m"<<Log::endl;

    // write results
    // -------------
    if(Parallel::isMaster(comm))
    {
      logStatus<<"write variational equation to file <"<<fileNameOutVariational<<">"<<Log::endl;
      writeFileVariationalEquation(fileNameOutVariational, satellite, arcs);
    }

    if(Parallel::isMaster(comm) && !fileNameOutOrbit.empty())
    {
      logStatus<<"write orbit to file <"<<fileNameOutOrbit<<">"<<Log::endl;
      std::list<Arc> arcList;
      for(UInt arcNo=0; arcNo<arcs.size(); arcNo++)
        arcList.push_back(arcs.at(arcNo).orbitArc());
      InstrumentFile::write(fileNameOutOrbit, arcList);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

VariationalEquationArc PreprocessingVariationalEquation::computeArc(UInt arcNo)
{
  try
  {
    OrbitArc          orbit         = orbitFile.readArc(arcNo);
    StarCameraArc     starCamera    = starCameraFile.readArc(arcNo);
    AccelerometerArc  accelerometer = accelerometerFile.readArc(arcNo);
    Arc::checkSynchronized({orbit, starCamera, accelerometer});
    std::vector<Time> times         = orbit.times();
    const UInt        epochCount    = orbit.size();
    const Double      T             = (times.back()-times.at(0)).seconds();
    const Double      deltaT        = T/(epochCount-1);

    // computeEarthRotation
    std::vector<Rotary3d> rotEarth(times.size());
    for(UInt k=0; k<times.size(); k++)
      rotEarth.at(k) = earthRotation->rotaryMatrix(times.at(k));

    // computeGravityGradient
    std::vector<Tensor3d> tensor(orbit.size());
    for(UInt i=0; i<orbit.size(); i++)
    {
      Vector3d posEarth = rotEarth.at(i).rotate(orbit.at(i).position);
      tensor.at(i)      = rotEarth.at(i).inverseRotate(gradientfield->gravityGradient(orbit.at(i).time, posEarth));
    }

    // computeForce
    std::vector<Vector3d> force(orbit.size());
    for(UInt k=0; k<orbit.size(); k++)
    {
      // forces are returned in TRF, rotate to CRF
      force.at(k) = rotEarth.at(k).inverseRotate(forces->acceleration(satellite, orbit.at(k).time, orbit.at(k).position, orbit.at(k).velocity,
                                         starCamera.at(k).rotary, rotEarth.at(k), earthRotation, ephemerides));
      // accelerometer
      if(accelerometer.size())
        force.at(k) += starCamera.at(k).rotary.rotate(accelerometer.at(k).acceleration);
    }

    if(parameterAcceleration)
    {
      parameterAcceleration->setIntervalArc(times.at(0), times.back()+seconds2time(deltaT));
      const UInt countA = parameterAcceleration->parameterCount();
      const UInt countB = parameterAcceleration->parameterCountArc();

      for(UInt k=0; k<times.size(); k++)
      {
        Matrix A(3, countA);
        Matrix B(3, countB);
        parameterAcceleration->compute(satellite, times.at(k), orbit.at(k).position, orbit.at(k).velocity,
                                    starCamera.at(k).rotary, rotEarth.at(k), ephemerides, A, B);
        Vector g2(3);
        if(countA) matMult(1., A, parameter.row(0, countA), g2);
        if(countB) matMult(1., B, parameter.row(arcNo*(countB+6)+countA, countB), g2);
        force.at(k) += rotEarth.at(k).inverseRotate(Vector3d(g2(0), g2(1), g2(2)));
      }
    }

    // sort into vector
    Vector g(3*epochCount);
    for(UInt k=0; k<epochCount; k++)
    {
      g(3*k+0) = force.at(k).x();
      g(3*k+1) = force.at(k).y();
      g(3*k+2) = force.at(k).z();
    }

    // =============================================

    // approx position
    // ---------------
    Vector posApprox(3*epochCount);
    for(UInt i=0; i<epochCount; i++)
    {
      posApprox(3*i+0) = orbit.at(i).position.x();
      posApprox(3*i+1) = orbit.at(i).position.y();
      posApprox(3*i+2) = orbit.at(i).position.z();
    }

    // Elliptical reference orbit
    // --------------------------
    Vector posRef(3*epochCount);
    Vector velRef(3*epochCount);
    if(GM>0) // useEnke
    {
      // Collect orbit epochs for mean orbit estimation. Try to use 100 evenly
      // spaced epochs along the arc, or all of the epochs if the arc is shorter.
      UInt stagger = (orbit.size() < 100) ? 1 : orbit.size()/100;

      OrbitArc sparseOrbit;
      for(UInt k=0; k<orbit.size(); k += stagger)
        sparseOrbit.push_back(orbit.at(k));

      // Osculating orbit at first epoch
      Vector3d velocity = sparseOrbit.at(0).velocity;
      if(velocity.r()==0)
        velocity = (1/(sparseOrbit.at(1).time-sparseOrbit.at(0).time).seconds()) * (sparseOrbit.at(1).position-sparseOrbit.at(0).position);

      // We use an equinoctial orbit because Kepler is not accurate enough
      Equinoctial referenceOrbit(sparseOrbit.at(0).time, sparseOrbit.at(0).position, velocity, GM);

      UInt   iter    = 0;
      UInt   maxIter = 100;
      Vector dx;
      do
      {
        // We will fit to positions only
        Vector obs(3*sparseOrbit.size());
        for(UInt i=0; i<sparseOrbit.size(); i++)
        {
          const Vector3d dpos = sparseOrbit.at(i).position - referenceOrbit.position(sparseOrbit.at(i).time);
          obs(3*i+0) = dpos.x();
          obs(3*i+1) = dpos.y();
          obs(3*i+2) = dpos.z();
        }

        Matrix A(3*sparseOrbit.size(), 6);
        for(UInt i=0; i<sparseOrbit.size(); i++)
          copy(referenceOrbit.stateTransitionMatrix(sparseOrbit.at(i).time).row(0,3), A.row(3*i,3));

        dx = leastSquares(A, obs);
        referenceOrbit.a  += dx(0);
        referenceOrbit.h  += dx(1);
        referenceOrbit.k  += dx(2);
        referenceOrbit.p  += dx(3);
        referenceOrbit.q  += dx(4);
        referenceOrbit.l0 += dx(5);
      }
      while((++iter < maxIter) && (maxabs(dx) > 1e-9));

      if(iter == maxIter)
        logWarning<<"no convergence in equinoctial fit"<< Log::endl;

      for(UInt i=0; i<epochCount; i++)
      {
        OrbitEpoch epoch;
        referenceOrbit.orbit(times.at(i), epoch.position, epoch.velocity, epoch.acceleration);

        // acceleration reduced by reference acceleration
        g(3*i+0) -= epoch.acceleration.x();
        g(3*i+1) -= epoch.acceleration.y();
        g(3*i+2) -= epoch.acceleration.z();

        // Observed position and velocity
        posRef(3*i+0) = epoch.position.x();
        posRef(3*i+1) = epoch.position.y();
        posRef(3*i+2) = epoch.position.z();

        velRef(3*i+0) = epoch.velocity.x();
        velRef(3*i+1) = epoch.velocity.y();
        velRef(3*i+2) = epoch.velocity.z();
      }
    } // if(useEnke)

    // =============================================

    // State matrix: position partial derivatives with respect to boundary values
    // --------------------------------------------------------------------------
    Matrix PosState(3*epochCount, 6);
    for(UInt i=0; i<epochCount; i++)
    {
      PosState(3*i+0, 0) = PosState(3*i+1, 1) = PosState(3*i+2, 2) = 1.; // d(pos)/d(pos0)
      PosState(3*i+0, 3) = PosState(3*i+1, 4) = PosState(3*i+2, 5) = i*deltaT/T; //(times.at(i)-times.at(0)).seconds()/T;     // d(pos)/d(vel0)
    }

    // Position
    // --------
    Matrix pos0  = integrate2Position(deltaT, g);
    Vector state = leastSquares(Matrix(PosState), posApprox-posRef-pos0);
    matMult(1, PosState, state, pos0);
    pos0 = solve(deltaT, tensor, pos0 - (posApprox-posRef)) + posApprox;

    // =============================================

    // State with indirect effect
    // --------------------------
    PosState = solve(deltaT, tensor, PosState);

    Matrix AccState(3*epochCount, 6);
    for(UInt i=0; i<epochCount; i++)
      matMult(1., tensor.at(i).matrix(), PosState.row(3*i,3), AccState.row(3*i,3));

    Matrix VelState(3*epochCount, 6);
    for(UInt i=0; i<epochCount; i++)
    {
      VelState(3*i+0, 3) = VelState(3*i+1, 4) = VelState(3*i+2, 5) = 1/T;
    }
    // indirect effect
    axpy(1., integrate2Velocity(deltaT, AccState), VelState);

    // =============================================

    // Acceleration
    // ------------
    Matrix acc0 = g;
    Matrix posDelta = pos0 - posApprox;
    for(UInt i=0; i<epochCount; i++)
      matMult(1., tensor.at(i).matrix(), posDelta.row(3*i,3), acc0.row(3*i,3));

    // Velocity
    // --------
    Matrix vel0 = integrate2Velocity(deltaT, acc0);
    for(UInt i=0; i<epochCount; i++)
      for(UInt j=0; j<vel0.columns(); j++)
        axpy(1/T, state.row(3,3), vel0.slice(3*i,j,3,1));
    vel0 += velRef; // velRef is either 0, or when useEnke, the velocity from equinoctial equations.

    // fit to reference orbit
    // ----------------------
    {
      Vector state = leastSquares(Matrix(PosState), posApprox-pos0);
      matMult(1, PosState, state, pos0);
      matMult(1, VelState, state, vel0);
    }

    // =============================================

    VariationalEquationArc arc;
    arc.times    = times;
    arc.pos0     = pos0;
    arc.vel0     = vel0;
    arc.PosState = PosState;
    arc.VelState = VelState;
    arc.rotEarth = rotEarth; // CRF -> TRF
    arc.rotSat.resize(orbit.size());
    for(UInt k=0; k<orbit.size(); k++)
      arc.rotSat.at(k) = starCamera.at(k).rotary; // Sat -> CRF

    // for statistics
    maxPosDiff = std::max(maxPosDiff, maxabs(arc.pos0-posApprox));

    return arc;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Matrix PreprocessingVariationalEquation::integrate2Position(Double deltaT, const_MatrixSliceRef g) const
{
  try
  {
    const UInt epochCount = g.rows()/3;
    Matrix result(g.rows(), g.columns());  // \int_0^t (t-t') f(t') dt'

    UInt idx       = 0;
    UInt evalPoint = 0;
    Matrix int_g(3, g.columns());   // \int_0^t f(t') dt'
    Matrix int_tg(3, g.columns());  // \int_0^t t'f(t') dt'
    for(UInt i=0; i<epochCount-1; i++)
    {
      Matrix tmp(3, g.columns());
      for(UInt k=0; k<coeff_g.at(idx).rows(); k++)
        axpy(deltaT * coeff_g.at(idx)(k), g.row(3*(i+k-evalPoint),3), tmp);
      for(UInt k=0; k<coeff_tg.at(idx).rows(); k++)
        axpy(deltaT*deltaT * coeff_tg.at(idx)(k), g.row(3*(i+k-evalPoint),3), int_tg);
      axpy(i*deltaT, tmp, int_tg);
      int_g += tmp;

      // select polynomial coefficients
      if(i<integrationDegree/2)
        idx++, evalPoint++;
      if(i>=epochCount-1-integrationDegree/2-1)
        idx++, evalPoint++;

      axpy((i+1)*deltaT, int_g,  result.row(3*(i+1),3)); //  t * \int_0^t f(t') dt'
      axpy(-1.,          int_tg, result.row(3*(i+1),3)); // -1 * \int_0^t t'f(t') dt'
    } // for(i)

    return result;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Matrix PreprocessingVariationalEquation::integrate2Velocity(Double deltaT, const_MatrixSliceRef g) const
{
  try
  {
    const UInt epochCount = g.rows()/3;
    Matrix result(g.rows(), g.columns());

    UInt idx       = 0;
    UInt evalPoint = 0;
    for(UInt i=0; i<epochCount-1; i++)
    {
      copy(result.row(3*i,3), result.row(3*(i+1),3));
      for(UInt k=0; k<coeff_g.at(idx).rows(); k++)
        axpy(deltaT*coeff_g.at(idx)(k), g.row(3*(i+k-evalPoint),3), result.row(3*(i+1),3));

      // select polynomial coefficients
      if(i<integrationDegree/2)
        idx++, evalPoint++;
      if(i>=epochCount-1-integrationDegree/2-1)
        idx++, evalPoint++;
    } // for(i)

    return result;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

// indirect effect (dependency of gravity due to position change)
// x = (I-integrate(T))^-1 l
// with the iterative solver BiCGSTAB
// BiCGSTAB (http://en.wikipedia.org/wiki/Biconjugate_gradient_stabilized_method)
Matrix PreprocessingVariationalEquation::solve(Double deltaT, const std::vector<Tensor3d> &tensor, const_MatrixSliceRef l) const
{
  try
  {
    // Approximate solution
    // --------------------
    Matrix x = approxInverse(deltaT, tensor, l);
    Matrix delta;

    for(UInt iter2=0; iter2<50; iter2++)
    {
      Matrix r  = l - refine(deltaT, tensor, x);
      Matrix r0 = r;
      const UInt rhsCount = x.columns();
      Vector rho(rhsCount), alpha(rhsCount), omega(rhsCount);
      for(UInt i=0; i<rhsCount; i++)
        rho(i) = alpha(i) = omega(i) = 1.;
      Matrix v(x.rows(), rhsCount);
      Matrix p(x.rows(), rhsCount);

      for(UInt iter=0; iter<1; iter++)
      {
        for(UInt i=0; i<rhsCount; i++)
        {
          Double rhoOld = rho(i);
          rho(i) = inner(r0.column(i), r.column(i));
          Double beta = (rho(i)/rhoOld)*(alpha(i)/omega(i));
          axpy(-omega(i), v.column(i), p.column(i));
          p.column(i) *= beta;
          axpy(1., r.column(i), p.column(i));
        }
        Matrix pHat = approxInverse(deltaT, tensor, p);
        v = refine(deltaT, tensor, pHat);
        for(UInt i=0; i<rhsCount; i++)
        {
          alpha(i) = rho(i)/inner(r0.column(i), v.column(i));
          axpy(-alpha(i), v.column(i), r.column(i));
        }
        Matrix sHat = approxInverse(deltaT, tensor, r);
        Matrix t = refine(deltaT, tensor, sHat);
        Matrix xOld = x;
        for(UInt i=0; i<rhsCount; i++)
        {
          omega(i) = inner(t.column(i), r.column(i))/inner(t.column(i), t.column(i));
          axpy( alpha(i), pHat.column(i), x.column(i));
          axpy( omega(i), sHat.column(i), x.column(i));
          axpy(-omega(i), t.column(i), r.column(i));
        }
        delta = x-xOld;
        if(maxabs(delta) < 1e-7)
          break;
      }  // for(iter)

      if(maxabs(delta) < 1e-7)
        break;
    } // for(iter2)

    return x;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

// y = (I-K*T)x = x - integrate(T*x)
Matrix PreprocessingVariationalEquation::refine(Double deltaT, const std::vector<Tensor3d> &tensor, const_MatrixSliceRef x) const
{
  try
  {
    Matrix y  = x;
    Matrix Tx(x.rows(), x.columns());
    for(UInt i=0; i<tensor.size(); i++)
      matMult(1., tensor.at(i).matrix(), x.row(3*i,3), Tx.row(3*i,3));
    y -= integrate2Position(deltaT, Tx);
    return y;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

// y = (I-K*T)^-1 x
Matrix PreprocessingVariationalEquation::approxInverse(Double deltaT, const std::vector<Tensor3d> &tensor, const_MatrixSliceRef x) const
{
  try
  {
    Matrix y   = x;
    Matrix tmp = tensor.at(0).matrix() * y.row(0,3);
    Matrix Ty  = 0.5*deltaT*deltaT * tmp;
    Matrix tTy = -deltaT*deltaT/6  * tmp;

    for(UInt i=1; i<tensor.size(); i++)
    {
      axpy(i,  Ty, y.row(3*i,3));
      axpy(1, tTy, y.row(3*i,3));

      Matrix K(3,3);
      K(0,0) = K(1,1) = K(2,2) = 1.;
      axpy(-1./6.*deltaT*deltaT, tensor.at(i).matrix(), K);
      solveInPlace(K, y.row(3*i,3));

      const Matrix tmp = tensor.at(i).matrix() * y.row(3*i,3);
      axpy( deltaT*deltaT,   tmp, Ty);
      axpy(-deltaT*deltaT*i, tmp, tTy);
    }
    return y;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
