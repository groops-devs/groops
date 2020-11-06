/***********************************************/
/**
* @file observationMiscPodVariational.cpp
*
* @brief Precise Orbit data (variational equations).
*
* @author Torsten Mayer-Guerr
* @date 2015-06-02
*
*/
/***********************************************/

#include "base/import.h"
#include "files/fileInstrument.h"
#include "misc/observation/variationalEquationFromFile.h"
#include "observationMiscPodVariational.h"

/***********************************************/

ObservationMiscPodVariational::ObservationMiscPodVariational(Config &config)
{
  try
  {
    FileName              fileNamePod;
    FileName              fileNameVariational;
    UInt                  integrationDegree;
    UInt                  interpolationDegree;

    renameDeprecatedConfig(config, "representation", "parametrizationGravity",      date2time(2020, 6, 3));
    renameDeprecatedConfig(config, "parameter",      "parametrizationAcceleration", date2time(2020, 6, 3));

    if(readConfigSequence(config, "rightHandSide", Config::MUSTSET, "", "input for observation vectors"))
    {
      readConfig(config, "inputfileOrbit", fileNamePod, Config::MUSTSET, "", "kinematic positions as observations");
      endSequence(config);
    }
    readConfig(config, "inputfileVariational",        fileNameVariational,   Config::MUSTSET,  "",    "approximate position and integrated state matrix");
    readConfig(config, "ephemerides",                 ephemerides,           Config::OPTIONAL, "jpl", "");
    readConfig(config, "parametrizationGravity",      parameterGravity,      Config::DEFAULT,  "",    "gravity field parametrization");
    readConfig(config, "parametrizationAcceleration", parameterAcceleration, Config::DEFAULT,  "",    "orbit force parameters");
    readConfig(config, "integrationDegree",           integrationDegree,     Config::DEFAULT,  "7",   "integration of forces by polynomial approximation of degree n");
    readConfig(config, "interpolationDegree",         interpolationDegree,   Config::DEFAULT,  "7",   "orbit interpolation by polynomial approximation of degree n");
    readConfig(config, "accelerateComputation",       accelerateComputation, Config::DEFAULT,  "0",   "acceleration of computation by transforming the observations");
    if(isCreateSchema(config)) return;

    // =======================

    // init
    // ----
    podFile.open(fileNamePod);
    countArc = podFile.arcCount();
    variationalEquation.open(fileNameVariational, parameterGravity, parameterAcceleration, std::vector<Time>(), ephemerides, integrationDegree);
    polynomial.init(interpolationDegree);

    // =======================

    // count parameters
    // ----------------
    countAParameter = variationalEquation.parameterCount();
    gravityCount    = variationalEquation.parameterCountGravity();
    idxGravity      = 0;
    idxState        = gravityCount;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void ObservationMiscPodVariational::setInterval(const Time &timeStart, const Time &timeEnd)
{
  try
  {
    parameterGravity->setInterval(timeStart, timeEnd);
    parameterAcceleration->setInterval(timeStart, timeEnd);
    variationalEquation.computeIndices();

    // count parameters
    // ----------------
    countAParameter = variationalEquation.parameterCount();
    gravityCount    = variationalEquation.parameterCountGravity();
    idxGravity      = 0;
    idxState        = gravityCount;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void ObservationMiscPodVariational::parameterName(std::vector<ParameterName> &name) const
{
  try
  {
    variationalEquation.parameterName(name);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

ObservationMiscPod::Arc ObservationMiscPodVariational::computeArc(UInt arcNo, CovariancePodPtr covPod)
{
  try
  {
    // read POD observations
    // ---------------------
    OrbitArc pod = podFile.readArc(arcNo);

    const UInt obsCount = 3*pod.size();
    if(obsCount == 0)
      return Arc();

    const std::vector<Time> timePod = pod.times();
    VariationalEquationFromFile::ObservationEquation eqn = variationalEquation.integrateArc(timePod.at(0), timePod.back(), TRUE/*position*/, FALSE/*velocity*/);

    // kinematic orbit observations
    // ----------------------------
    Matrix l(obsCount, 1);
    for(UInt k=0; k<obsCount/3; k++)
    {
      l(3*k+0,0) = pod.at(k).position.x();
      l(3*k+1,0) = pod.at(k).position.y();
      l(3*k+2,0) = pod.at(k).position.z();
    }
    l -= polynomial.interpolate(timePod, eqn.times, eqn.pos0, 3);

    // =============================================

    Matrix A;
    if(!accelerateComputation || (l.rows() <= eqn.PosDesign.rows()))
    {
      // design matrix
      // -------------
      A = polynomial.interpolate(timePod, eqn.times, eqn.PosDesign, 3);

      // decorrelation
      // -------------
      if(covPod)
        covPod->decorrelate(arcNo, pod, {l, A});
    }
    else
    {
      A = eqn.PosDesign;
      Matrix B = polynomial.interpolate(timePod, eqn.times, identityMatrix(A.rows()), 3);
      if(covPod)
        covPod->decorrelate(arcNo, pod, {l, B});
      const Vector tau = QR_decomposition(B);
      QTransMult(B, tau, l);
      triangularMult(1., B.row(0, B.columns()), A);
    }

    // =============================================

    Arc observationArc;
    observationArc.l     = l;
    observationArc.A     = A;
    observationArc.times = timePod;
    observationArc.pod   = pod;
    return observationArc;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
