/***********************************************/
/**
* @file observationMiscSstVariational.cpp
*
* @brief Satellite to satellite tracking (Variational equations).
*
* @author Torsten Mayer-Guerr
* @date 2015-06-02
*
*/
/***********************************************/

#include "base/import.h"
#include "base/polynomial.h"
#include "files/fileInstrument.h"
#include "classes/parametrizationGravity/parametrizationGravity.h"
#include "classes/parametrizationAcceleration/parametrizationAcceleration.h"
#include "classes/parametrizationSatelliteTracking/parametrizationSatelliteTracking.h"
#include "misc/observation/variationalEquation.h"
#include "misc/observation/variationalEquationFromFile.h"
#include "observationMiscSstVariational.h"

/***********************************************/

ObservationMiscSstVariational::ObservationMiscSstVariational(Config &config)
{
  try
  {
    std::vector<FileName> fileNameSst;
    FileName              fileNamePod1, fileNamePod2;
    FileName              fileNameVariational1, fileNameVariational2;
    UInt                  integrationDegree;
    sstType = 0;

    renameDeprecatedConfig(config, "representation", "parametrizationGravity",       date2time(2020, 6, 3));
    renameDeprecatedConfig(config, "parameter1",     "parametrizationAcceleration1", date2time(2020, 6, 3));
    renameDeprecatedConfig(config, "parameter2",     "parametrizationAcceleration2", date2time(2020, 6, 3));
    renameDeprecatedConfig(config, "parameterSst",   "parametrizationSst",           date2time(2020, 6, 3));

    if(readConfigSequence(config, "rightHandSide", Config::MUSTSET, "", "input for observation vectors"))
    {
      readConfig(config, "inputfileSatelliteTracking", fileNameSst,  Config::MUSTSET,  "", "ranging observations and corrections");
      readConfig(config, "inputfileOrbit1",            fileNamePod1, Config::OPTIONAL, "", "kinematic positions of satellite A as observations");
      readConfig(config, "inputfileOrbit2",            fileNamePod2, Config::OPTIONAL, "", "kinematic positions of satellite B as observations");
      endSequence(config);
    }
    std::string choice;
    if(readConfigChoice(config, "sstType", choice, Config::MUSTSET, "", ""))
    {
      if(readConfigChoiceElement(config, "range",     choice, "")) sstType = 0;
      if(readConfigChoiceElement(config, "rangeRate", choice, "")) sstType = 1;
      if(readConfigChoiceElement(config, "none",      choice, "")) sstType = 99;
      endChoice(config);
    }
    readConfig(config, "inputfileVariational1",        fileNameVariational1,   Config::MUSTSET,  "",    "approximate position and integrated state matrix");
    readConfig(config, "inputfileVariational2",        fileNameVariational2,   Config::MUSTSET,  "",    "approximate position and integrated state matrix");
    readConfig(config, "ephemerides",                  ephemerides,            Config::OPTIONAL, "jpl", "");
    readConfig(config, "parametrizationGravity",       parameterGravity,       Config::DEFAULT,  "",    "gravity field parametrization");
    readConfig(config, "parametrizationAcceleration1", parameterAcceleration1, Config::DEFAULT,  "",    "orbit1 force parameters");
    readConfig(config, "parametrizationAcceleration2", parameterAcceleration2, Config::DEFAULT,  "",    "orbit2 force parameters");
    readConfig(config, "parametrizationSst",           parameterSst,           Config::DEFAULT,  "",    "satellite tracking parameter");
    readConfig(config, "integrationDegree",            integrationDegree,      Config::DEFAULT,  "7",   "integration of forces by polynomial approximation of degree n");
    readConfig(config, "interpolationDegree",          interpolationDegree,    Config::DEFAULT,  "7",   "orbit interpolation by polynomial approximation of degree n");
    if(isCreateSchema(config)) return;

    // =======================

    // init
    // ----
    sstFile.resize(fileNameSst.size());
    for(UInt i=0; i<fileNameSst.size(); i++)
      sstFile.at(i) = InstrumentFile::newFile(fileNameSst.at(i));
    pod1File.open(fileNamePod1);
    pod2File.open(fileNamePod2);

    // test instrument files
    // ---------------------
    for(UInt i=1; i<sstFile.size(); i++)
      InstrumentFile::checkArcCount({*sstFile.at(i), *sstFile.at(0)});
    if(sstFile.size())
      InstrumentFile::checkArcCount({pod1File, *sstFile.at(0)});
    InstrumentFile::checkArcCount({pod2File, pod1File});

    countArc = std::max(pod1File.arcCount(), pod2File.arcCount());
    if(sstFile.size()) countArc = std::max(countArc, sstFile.at(0)->arcCount());
    computeVelocity = (sstType==1);

    variationalEquation1.open(fileNameVariational1, parameterGravity, parameterAcceleration1, std::vector<Time>(), ephemerides, integrationDegree);
    variationalEquation2.open(fileNameVariational2, parameterGravity, parameterAcceleration2, std::vector<Time>(), ephemerides, integrationDegree);

    // =======================

    // count parameters
    // ----------------
    gravityCount = variationalEquation1.parameterCountGravity();
    state1Count  = variationalEquation1.parameterCount() - gravityCount;
    state2Count  = variationalEquation2.parameterCount() - gravityCount;

    countAParameter = 0;
    idxGravity  = countAParameter; countAParameter += gravityCount;
    idxState1   = countAParameter; countAParameter += state1Count;
    idxState2   = countAParameter; countAParameter += state2Count;
    idxSstPara  = countAParameter; countAParameter += parameterSst->parameterCount();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Bool ObservationMiscSstVariational::setInterval(const Time &timeStart, const Time &timeEnd)
{
  try
  {
    Bool change = FALSE;
    change = parameterGravity->setInterval(timeStart, timeEnd)       || change;
    change = parameterAcceleration1->setInterval(timeStart, timeEnd) || change;
    change = parameterAcceleration2->setInterval(timeStart, timeEnd) || change;
    change = parameterSst->setInterval(timeStart, timeEnd)           || change;
    if(!change)
      return FALSE;
    variationalEquation1.computeIndices();
    variationalEquation2.computeIndices();

    // count parameters
    // ----------------
    gravityCount = variationalEquation1.parameterCountGravity();
    state1Count  = variationalEquation1.parameterCount() - gravityCount;
    state2Count  = variationalEquation2.parameterCount() - gravityCount;

    countAParameter = 0;
    idxGravity  = countAParameter; countAParameter += gravityCount;
    idxState1   = countAParameter; countAParameter += state1Count;
    idxState2   = countAParameter; countAParameter += state2Count;
    idxSstPara  = countAParameter; countAParameter += parameterSst->parameterCount();

    return change;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void ObservationMiscSstVariational::parameterName(std::vector<ParameterName> &name) const
{
  try
  {
    parameterGravity->parameterName(name);

    std::vector<ParameterName> name1;
    variationalEquation1.parameterName(name1);
    const std::string satelliteName1 = variationalEquation1.satellite() ? variationalEquation1.satellite()->satelliteName : "satellite1";
    for(UInt i=gravityCount; i<name1.size(); i++)
    {
      name1.at(i).object = satelliteName1;
      name.push_back(name1.at(i));
    }

    std::vector<ParameterName> name2;
    variationalEquation2.parameterName(name2);
    const std::string satelliteName2 = variationalEquation2.satellite() ? variationalEquation2.satellite()->satelliteName : "satellite2";
    for(UInt i=gravityCount; i<name2.size(); i++)
    {
      name2.at(i).object = satelliteName2;
      name.push_back(name2.at(i));
    }

    parameterSst->parameterName(name);
    for(UInt i=name.size()-parameterSst->parameterCount(); i<name.size(); i++)
      name.at(i).object = satelliteName1+"."+satelliteName2;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

ObservationMiscSstVariational::Arc ObservationMiscSstVariational::computeArc(UInt arcNo, CovarianceSstPtr covSst, CovariancePodPtr covPod1, CovariancePodPtr covPod2)
{
  try
  {
    // read SST observations
    // ---------------------
    std::vector<SatelliteTrackingArc> sst(sstFile.size());
    for(UInt k=0; k<sstFile.size(); k++)
      sst.at(k) = sstFile.at(k)->readArc(arcNo);
    for(UInt k=1; k<sst.size(); k++)
      ::Arc::checkSynchronized({sst.at(0), sst.at(k)});
    const std::vector<Time> timesSst = (sst.size() && (sstType!=99)) ? sst.at(0).times() : std::vector<Time>();

    // read POD observations
    // ---------------------
    OrbitArc pod1 = pod1File.readArc(arcNo);
    OrbitArc pod2 = pod2File.readArc(arcNo);
    const std::vector<Time> timesPod1 = pod1.times();
    const std::vector<Time> timesPod2 = pod2.times();

    // Init variational equations
    // --------------------------
    Time timeStart(std::numeric_limits<Int>::max(), 0), timeEnd;
    if(timesSst.size())  {timeStart = std::min(timeStart, timesSst.front());  timeEnd = std::max(timeEnd, timesSst.back());}
    if(timesPod1.size()) {timeStart = std::min(timeStart, timesPod1.front()); timeEnd = std::max(timeEnd, timesPod1.back());}
    if(timesPod2.size()) {timeStart = std::min(timeStart, timesPod2.front()); timeEnd = std::max(timeEnd, timesPod2.back());}
    const VariationalEquationArc &arc1 = variationalEquation1.getArc(timeStart);
    const VariationalEquationArc &arc2 = variationalEquation2.getArc(timeStart);

    if((arc1.times != arc2.times) || (timeStart < arc1.times.front()) || (arc1.times.back() < timeEnd))
      throw(Exception("no variational arc found for ["+timeStart.dateTimeStr()+", "+timeEnd.dateTimeStr()+"]"));

    // integration times (and interpolation intervals inbetween)
    std::vector<Time> times;
    const UInt epochStart = std::distance(arc1.times.begin(), std::upper_bound(arc1.times.begin(), arc1.times.end(), timeStart))-1;
    for(UInt idEpoch=epochStart; (idEpoch < arc1.times.size()) && (arc1.times.at(idEpoch) <= timeEnd); idEpoch++)
      times.push_back(arc1.times.at(idEpoch));

    // count observations and calculate index
    // --------------------------------------
    UInt obsCount = 0;
    const UInt idxSst  = obsCount; obsCount += timesSst.size();
    const UInt idxPod1 = obsCount; obsCount += 3*timesPod1.size();
    const UInt idxPod2 = obsCount; obsCount += 3*timesPod2.size();

    // arc related parameters (in matrix B)
    // ------------------------------------
    if(timesSst.size()) parameterSst->setIntervalArc(timesSst.front(), timesSst.back()+medianSampling(timesSst));
    UInt countBParameter = 0;
    const UInt idxSstParaArc = countBParameter; countBParameter += parameterSst->parameterCountArc();
    if(obsCount <= countBParameter)
      return Arc();

    Matrix l(obsCount, 1);
    Matrix A(obsCount, countAParameter);
    Matrix B(obsCount, countBParameter);

    // needed for parameterSst
    Vector sstComputed(timesSst.size());
    Vector sstPos1(3*timesSst.size()), sstPos2(3*timesSst.size());
    Vector sstVel1(3*timesSst.size()), sstVel2(3*timesSst.size());

    // observation equations within each integration interval
    // ------------------------------------------------------
    UInt idEpochSst=0, idEpochPod1=0, idEpochPod2=0;
    std::vector<VariationalEquationFromFile::ObservationEquation> eqn1(interpolationDegree+1);
    std::vector<VariationalEquationFromFile::ObservationEquation> eqn2(interpolationDegree+1);
    for(UInt idInterval=0; idInterval<times.size()-1; idInterval++)
    {
      // integration of eqn1, eqn2  which are cyclicly updated
      // -----------------------------------------------------
      const UInt start = idInterval ? std::max(interpolationDegree+1, idInterval+(interpolationDegree+1)/2) : 0;
      const UInt end   = std::min(std::max(idInterval+(interpolationDegree+1)/2+1, interpolationDegree+1), times.size());
      for(UInt i=start; i<end; i++)
      {
        eqn1.at(i%eqn1.size()) = variationalEquation1.integrateArc(times.at(i), times.at(i), TRUE/*position*/, computeVelocity);
        eqn2.at(i%eqn2.size()) = variationalEquation2.integrateArc(times.at(i), times.at(i), TRUE/*position*/, computeVelocity);
      }

      // satellite to satellite tracking
      // -------------------------------
      for(; (idEpochSst < timesSst.size()) && (timesSst.at(idEpochSst) <= times.at(idInterval+1)); idEpochSst++)
      {
        Matrix pos1, vel1, PosDesign1, VelDesign1;
        Matrix pos2, vel2, PosDesign2, VelDesign2;
        interpolate(timesSst.at(idEpochSst), eqn1, pos1, vel1, PosDesign1, VelDesign1, computeVelocity, interpolationDegree);
        interpolate(timesSst.at(idEpochSst), eqn2, pos2, vel2, PosDesign2, VelDesign2, computeVelocity, interpolationDegree);

        copy(pos1, sstPos1.row(3*idEpochSst, 3));
        copy(pos2, sstPos2.row(3*idEpochSst, 3));
        if(computeVelocity)
        {
          copy(vel1, sstVel1.row(3*idEpochSst, 3));
          copy(vel2, sstVel2.row(3*idEpochSst, 3));
        }

        if(sstType == 0) // range
        {
          Vector3d e12(pos2-pos1);
          const Double rho = e12.normalize();

          sstComputed(idEpochSst) = rho;
          for(UInt k=0; k<sst.size(); k++)
            l(idxSst+idEpochSst, 0) += (sst.at(k).at(idEpochSst).range - sst.at(k).at(0).range);
          l(idxSst+idEpochSst, 0) -= (sstComputed(idEpochSst)-sstComputed(0));

          // dRange/dPos12
          Matrix e = e12.vector().trans();
          if(gravityCount)
            matMult( 1., e, PosDesign2.column(0, gravityCount)-PosDesign1.column(0, gravityCount), A.slice(idxSst+idEpochSst, idxGravity, 1, gravityCount));
          matMult(-1., e, PosDesign1.column(gravityCount, state1Count),  A.slice(idxSst+idEpochSst, idxState1,  1, state1Count));
          matMult( 1., e, PosDesign2.column(gravityCount, state2Count),  A.slice(idxSst+idEpochSst, idxState2,  1, state2Count));
        }
        else if(sstType == 1) // range rate
        {
          Vector3d e12(pos2-pos1);
          const Double rho = e12.normalize();
          Vector3d vel12(vel2-vel1);
          const Double drho = inner(e12, vel12);

          sstComputed(idEpochSst) = drho;
          for(UInt k=0; k<sst.size(); k++)
            l(idxSst+idEpochSst,0) += sst.at(k).at(idEpochSst).rangeRate;
          l(idxSst+idEpochSst,0) -= sstComputed(idEpochSst);

          // dRate/dPos12
          Matrix e = ((1/rho)*vel12 - (drho/rho) * e12).vector().trans();
          if(gravityCount)
            matMult( 1., e, PosDesign2.column(0, gravityCount)-PosDesign1.column(0, gravityCount), A.slice(idxSst+idEpochSst, idxGravity, 1, gravityCount));
          matMult(-1., e, PosDesign1.column(gravityCount, state1Count),  A.slice(idxSst+idEpochSst, idxState1,  1, state1Count));
          matMult( 1., e, PosDesign2.column(gravityCount, state2Count),  A.slice(idxSst+idEpochSst, idxState2,  1, state2Count));

          // dRate/dVel12 = e12
          e = e12.vector().trans();
          if(gravityCount)
            matMult( 1., e, VelDesign2.column(0, gravityCount)-VelDesign1.column(0, gravityCount), A.slice(idxSst+idEpochSst, idxGravity, 1, gravityCount));
          matMult(-1., e, VelDesign1.column(gravityCount, state1Count),  A.slice(idxSst+idEpochSst, idxState1,  1, state1Count));
          matMult( 1., e, VelDesign2.column(gravityCount, state2Count),  A.slice(idxSst+idEpochSst, idxState2,  1, state2Count));
        }
      } // for(idEpochSst)

      // POD 1
      // -----
      for(; (idEpochPod1 < timesPod1.size()) && (timesPod1.at(idEpochPod1) <= times.at(idInterval+1)); idEpochPod1++)
      {
        Matrix pos1, vel1, PosDesign1, VelDesign1;
        interpolate(timesPod1.at(idEpochPod1), eqn1, pos1, vel1, PosDesign1, VelDesign1, FALSE, interpolationDegree);

        l(idxPod1+3*idEpochPod1+0, 0) = pod1.at(idEpochPod1).position.x();
        l(idxPod1+3*idEpochPod1+1, 0) = pod1.at(idEpochPod1).position.y();
        l(idxPod1+3*idEpochPod1+2, 0) = pod1.at(idEpochPod1).position.z();
        l.row(idxPod1+3*idEpochPod1, 3) -= pos1;

        if(gravityCount)
          copy(PosDesign1.column(0, gravityCount),         A.slice(idxPod1+3*idEpochPod1, idxGravity, 3, gravityCount));
        copy(PosDesign1.column(gravityCount, state1Count), A.slice(idxPod1+3*idEpochPod1, idxState1,  3, state1Count));
      } // for(idEpochPod1)

      // POD 2
      // -----
      for(; (idEpochPod2 < timesPod2.size()) && (timesPod2.at(idEpochPod2) <= times.at(idInterval+1)); idEpochPod2++)
      {
        Matrix pos2, vel2, PosDesign2, VelDesign2;
        interpolate(timesPod2.at(idEpochPod2), eqn2, pos2, vel2, PosDesign2, VelDesign2, FALSE, interpolationDegree);

        l(idxPod2+3*idEpochPod2+0, 0) = pod2.at(idEpochPod2).position.x();
        l(idxPod2+3*idEpochPod2+1, 0) = pod2.at(idEpochPod2).position.y();
        l(idxPod2+3*idEpochPod2+2, 0) = pod2.at(idEpochPod2).position.z();
        l.row(idxPod2+3*idEpochPod2, 3) -= pos2;

        if(gravityCount)
          copy(PosDesign2.column(0, gravityCount),         A.slice(idxPod2+3*idEpochPod2, idxGravity, 3, gravityCount));
        copy(PosDesign2.column(gravityCount, state2Count), A.slice(idxPod2+3*idEpochPod2, idxState2,  3, state2Count));
      } // for(idEpochPod2)
    } // for(idInterval)

    // =============================================

    if(parameterSst && timesSst.size())
    {
      parameterSst->compute(sstType, timesSst, sstComputed, sstPos1, sstPos2, sstVel1, sstVel2,
                            interpolateStarCamera(timesSst, arc1.times, arc1.rotSat, 1),  // linear interpolation
                            interpolateStarCamera(timesSst, arc2.times, arc2.rotSat, 1),
                            A.slice(idxSst, idxSstPara,    timesSst.size(), parameterSst->parameterCount()),
                            B.slice(idxSst, idxSstParaArc, timesSst.size(), parameterSst->parameterCountArc()));
    }

    // =============================================

    // decorrelation
    // -------------
    if(covSst  && timesSst.size())  covSst->decorrelate (arcNo, timesSst, {l.row(idxSst,  timesSst.size()),    A.row(idxSst,  timesSst.size()),    B.row(idxSst,  timesSst.size())});
    if(covPod1 && timesPod1.size()) covPod1->decorrelate(arcNo, pod1,     {l.row(idxPod1, 3*timesPod1.size()), A.row(idxPod1, 3*timesPod1.size()), B.row(idxPod1, 3*timesPod1.size())});
    if(covPod2 && timesPod2.size()) covPod2->decorrelate(arcNo, pod2,     {l.row(idxPod2, 3*timesPod2.size()), A.row(idxPod2, 3*timesPod2.size()), B.row(idxPod2, 3*timesPod2.size())});

   // =============================================

    Arc observationArc;
    observationArc.l     = l;
    observationArc.A     = A;
    observationArc.B     = B;
    observationArc.times = {timesSst, timesPod1, timesPod2};
    observationArc.pod1  = pod1;
    observationArc.pod2  = pod2;
    return observationArc;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void ObservationMiscSstVariational::interpolate(const Time &time, const std::vector<VariationalEquationFromFile::ObservationEquation> &eqn,
                                                Matrix &pos0, Matrix &vel0, Matrix &PosDesign, Matrix &VelDesign, Bool computeVelocity, UInt degree)
{
  try
  {
    // time matches? -> no interpolation needed
    for(auto &e : eqn)
      if(std::fabs((time-e.times.front()).seconds()) < 1e-9)
      {
        pos0      = e.pos0;
        PosDesign = e.PosDesign;
        if(computeVelocity) vel0      = e.vel0;
        if(computeVelocity) VelDesign = e.VelDesign;
        return;
      }

    // time does not match, interpolation required
    // -------------------------------------------
    // polynomial interpolation matrix
    Matrix W(degree+1, degree+1);
    for(UInt i=0; i<W.rows(); i++)
    {
      const Double dt = (eqn.at(i).times.front()-time).seconds();
      W(i, 0) = 1.0;
      for(UInt n=1; n<W.columns(); n++)
        W(i, n) = dt * W(i, n-1);
    }
    inverse(W); // interpolation weights are the first line of W^-1, because a_0 is the first column

    pos0      = W(0, 0) * eqn.front().pos0;
    PosDesign = W(0, 0) * eqn.front().PosDesign;
    if(computeVelocity) vel0      = W(0, 0) * eqn.front().vel0;
    if(computeVelocity) VelDesign = W(0, 0) * eqn.front().VelDesign;
    for(UInt k=1; k<eqn.size(); k++)
    {
      axpy(W(0, k), eqn.at(k).pos0,      pos0);
      axpy(W(0, k), eqn.at(k).PosDesign, PosDesign);
      if(computeVelocity) axpy(W(0, k), eqn.at(k).vel0,      vel0);
      if(computeVelocity) axpy(W(0, k), eqn.at(k).VelDesign, VelDesign);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

std::vector<Rotary3d> ObservationMiscSstVariational::interpolateStarCamera(const std::vector<Time> &timesNew, const std::vector<Time> &times, const std::vector<Rotary3d> &rot, UInt degree)
{
  try
  {
    Matrix A(4, rot.size());
    for(UInt i=0; i<rot.size(); i++)
      copy(rot.at(i).quaternion(), A.column(i));

    for(UInt i=1; i<A.columns(); i++)
      if(inner(A.column(i-1), A.column(i)) < 0)
        A.column(i) *= -1.;

    Polynomial polynomial(times, degree);
    A = polynomial.interpolate(timesNew, A.trans()).trans();

    std::vector<Rotary3d> rotNew(A.columns());
    for(UInt i=0; i<rotNew.size(); i++)
      rotNew.at(i) = Rotary3d(A.column(i)/norm(A.column(i)));

    return rotNew;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
