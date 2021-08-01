/***********************************************/
/**
* @file variationalEquation.cpp
*
* @brief Variational equations.
*
* @author Torsten Mayer-Guerr
* @date 2012-05-29
*
*/
/***********************************************/

#include "base/import.h"
#include "files/fileSatelliteModel.h"
#include "files/fileVariationalEquation.h"
#include "classes/parametrizationGravity/parametrizationGravity.h"
#include "classes/parametrizationAcceleration/parametrizationAcceleration.h"
#include "variationalEquation.h"

/***********************************************/

VariationalEquation::VariationalEquation() : parameterCount_(0), gravityCount(0), satCount(0), satArcCount(0)
{
}

/***********************************************/

void VariationalEquation::init(SatelliteModelPtr satellite, ParametrizationGravityPtr parameterGravity, ParametrizationAccelerationPtr parameterAcceleration,
                               const std::vector<Time> &stochasticPulse, EphemeridesPtr ephemerides, UInt integrationDegree)
{
  try
  {
    this->satellite             = satellite;
    this->parameterGravity      = parameterGravity;
    this->parameterAcceleration = parameterAcceleration;
    this->timePulse             = stochasticPulse;
    this->ephemerides           = ephemerides;
    this->integrationDegree     = integrationDegree;

    // count parameters
    // ----------------
    computeIndices();

    // integration in [idx,idx+1]
    // --------------------------
    if(coeffIntegral.size() != integrationDegree)
    {
      coeffIntegral.resize(integrationDegree);
      for(UInt idx=0; idx<coeffIntegral.size(); idx++)
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

        coeffIntegral.at(idx)  = Vector(W.rows());
        for(UInt i=0; i<W.rows(); i++)
          for(UInt n=0; n<W.columns(); n++)
            coeffIntegral.at(idx)(i)  += 1./(n+1.) * W(n,i);
      }
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void VariationalEquation::computeIndices()
{
  try
  {
    parameterCount_  = 0;

    idxGravity       = parameterCount_;
    gravityCount     = (parameterGravity) ? parameterGravity->parameterCount() : 0;
    parameterCount_ += gravityCount;

    idxSat           = parameterCount_;
    satCount         = (parameterAcceleration) ? parameterAcceleration->parameterCount() : 0;
    idxPulse         = parameterCount_ + satCount;
    satCount        += 3*timePulse.size();
    parameterCount_ += satCount;

    idxSatArc        = parameterCount_;
    satArcCount      = (parameterAcceleration) ? parameterAcceleration->parameterCountArc() + 6 : 6; // inclusive state vector
    parameterCount_ += satArcCount;

    idEpochAlpha = NULLINDEX; // current arc is invalid
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void VariationalEquation::setArc(VariationalEquationArc &arc)
{
  try
  {
    this->arc = arc;
    if(arc.times.size() == 0)
      throw(Exception("empty arc"));

    if(parameterAcceleration)
    {
      if(parameterAcceleration->setIntervalArc(arc.times.at(0), arc.times.back()+medianSampling(arc.times)))
        computeIndices();
    }

    initIntegration();

    // init pulses
    // -----------
    stochasticPulse.clear();
    stochasticPulse.resize(timePulse.size());
    for(UInt i=0; i<timePulse.size(); i++)
    {
      if(timePulse.at(i)>=arc.times.back()) break;
      if(timePulse.at(i)<=arc.times.at(0))  continue;
      UInt idEpoch = 0;
      while(timePulse.at(i)>=arc.times.at(idEpoch+1))
        idEpoch++;
      stochasticPulse.at(i) = Matrix(6, 3);
      stochasticPulse.at(i)(3,0) = stochasticPulse.at(i)(4,1) = stochasticPulse.at(i)(5,2) = 1e-6; // [mum/s]
      Matrix Z(6,6);
      copy(arc.PosState.row(3*idEpoch,3), Z.row(0,3));
      copy(arc.VelState.row(3*idEpoch,3), Z.row(3,3));
      solveInPlace(Z, stochasticPulse.at(i));
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void VariationalEquation::parameterName(std::vector<ParameterName> &name) const
{
  parameterNameGravity(name);
  parameterNameSatellite(name);
  parameterNameSatelliteArc(name);
}

/***********************************************/

void VariationalEquation::parameterNameGravity(std::vector<ParameterName> &name) const
{
  if(parameterGravity)
    parameterGravity->parameterName(name);
}

/***********************************************/

void VariationalEquation::parameterNameSatellite(std::vector<ParameterName> &name) const
{
  const std::string satelliteName = satellite ? satellite->satelliteName : "satellite";
  if(parameterAcceleration)
  {
    parameterAcceleration->parameterName(name);
    for(UInt i=name.size(); i-->name.size()-parameterAcceleration->parameterCount();)
      name.at(i).object = satelliteName;
  }

  for(UInt i=0; i<timePulse.size(); i++)
  {
    name.push_back(ParameterName(satelliteName, "stochasticPulse.x", "", timePulse.at(i)));
    name.push_back(ParameterName(satelliteName, "stochasticPulse.y", "", timePulse.at(i)));
    name.push_back(ParameterName(satelliteName, "stochasticPulse.z", "", timePulse.at(i)));
  }
}

/***********************************************/

void VariationalEquation::parameterNameSatelliteArc(std::vector<ParameterName> &name) const
{
  const std::string satelliteName = satellite ? satellite->satelliteName : "satellite";
  if(parameterAcceleration)
  {
    parameterAcceleration->parameterNameArc(name);
    for(UInt i=name.size(); i-->name.size()-parameterAcceleration->parameterCountArc();)
      name.at(i).object = satelliteName;
  }

  name.push_back(ParameterName(satelliteName, "position0.x"));
  name.push_back(ParameterName(satelliteName, "position0.y"));
  name.push_back(ParameterName(satelliteName, "position0.z"));
  name.push_back(ParameterName(satelliteName, "velocity0.x"));
  name.push_back(ParameterName(satelliteName, "velocity0.y"));
  name.push_back(ParameterName(satelliteName, "velocity0.z"));
}

/***********************************************/

void VariationalEquation::position(UInt idEpoch, MatrixSliceRef pos0, MatrixSliceRef PosDesign)
{
  try
  {
    if(arc.times.size() == 0)
      throw(Exception("no variational arc given"));
    if(idEpoch>=arc.times.size())
      throw(Exception("epoch outside arc"));

    computeAlpha(idEpoch);

    // gravity and satellite parameters
    // --------------------------------
    if(gravityCount+satCount)
      matMult(1., arc.PosState.row(3*idEpoch,3), Alpha.column(0, gravityCount+satCount), PosDesign.column(0, gravityCount+satCount));

    // stochastic pulses
    // -----------------
    for(UInt i=0; i<timePulse.size(); i++)
    {
      if(stochasticPulse.at(i).size()==0)
        continue;
      if(timePulse.at(i) >  arc.times.at(idEpoch))
        break;
      matMult(1., arc.PosState.row(3*idEpoch,3), stochasticPulse.at(i), PosDesign.column(3*i+idxPulse, 3));
    }

    // arc related parameters
    // ----------------------
    if(satArcCount-6)
      matMult(1., arc.PosState.row(3*idEpoch,3), Alpha.column(gravityCount+satCount, satArcCount-6), PosDesign.column(idxSatArc, satArcCount-6));

    // satellite state
    // ---------------
    copy(arc.PosState.row(3*idEpoch,3), PosDesign.column(idxSatArc+satArcCount-6, 6));
    copy(arc.pos0.row(3*idEpoch,3), pos0);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void VariationalEquation::velocity(UInt idEpoch, MatrixSliceRef vel0, MatrixSliceRef VelDesign)
{
  try
  {
    if(arc.times.size() == 0)
      throw(Exception("no variational arc given"));
    if(idEpoch>=arc.times.size())
      throw(Exception("epoch outside arc"));

    computeAlpha(idEpoch);

    // gravity and satellite parameters
    // --------------------------------
    if(gravityCount+satCount)
      matMult(1., arc.VelState.row(3*idEpoch,3), Alpha.column(0, gravityCount+satCount), VelDesign.column(0, gravityCount+satCount));

    // stochastic pulses
    // -----------------
    for(UInt i=0; i<timePulse.size(); i++)
    {
      if(stochasticPulse.at(i).size()==0)
        continue;
      if(timePulse.at(i) >  arc.times.at(idEpoch))
        break;
      matMult(1., arc.VelState.row(3*idEpoch,3), stochasticPulse.at(i), VelDesign.column(3*i+idxPulse, 3));
    }

    // arc related parameters
    // ----------------------
    if(satArcCount-6)
      matMult(1., arc.VelState.row(3*idEpoch,3), Alpha.column(gravityCount+satCount, satArcCount-6), VelDesign.column(idxSatArc, satArcCount-6));

    // satellite state
    // ---------------
    copy(arc.VelState.row(3*idEpoch,3), VelDesign.column(idxSatArc+satArcCount-6, 6));
    copy(arc.vel0.row(3*idEpoch,3), vel0);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void VariationalEquation::initIntegration()
{
  try
  {
    Alpha        = Matrix(6, parameterCount_-6); // without satellite state
    idEpochAlpha = 0;
    idCoeff      = 0;
    deltaT       = (arc.times.at(1) - arc.times.at(0)).seconds();
    Integrand.resize(integrationDegree+1);
    for(UInt i=0; i<Integrand.size(); i++)
      Integrand.at(i) = computeIntegrand(i);
    idIntegrand  = Integrand.size();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void VariationalEquation::computeAlpha(UInt idEpoch)
{
  try
  {
    if((idEpoch<idEpochAlpha) || (idEpochAlpha==NULLINDEX))
      initIntegration(); // start integration at beginning

    if(Alpha.size() == 0)
      return;

    while(idEpochAlpha < idEpoch)
    {
      idEpochAlpha++; // next epoch to integrate

      // integration
      for(UInt i=0; i<coeffIntegral.at(idCoeff).rows(); i++)
        axpy(deltaT*coeffIntegral.at(idCoeff)(i), Integrand.at((idIntegrand+i)%Integrand.size()), Alpha);

      if((idEpochAlpha<(integrationDegree+1)/2) || (idEpochAlpha+(integrationDegree+1)/2 >= arc.times.size()))
        idCoeff++; // select next integration interval within polynomial
      else
      {
        // shift integration polynomial: replace oldest integrand
        Integrand.at(idIntegrand % Integrand.size()) = computeIntegrand(idIntegrand);
        idIntegrand++;
      }
    } // while(idEpochAlpha < idEpoch)
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Matrix VariationalEquation::computeIntegrand(UInt idEpoch)
{
  try
  {
    if(Alpha.columns()==0)
      return Matrix(6,0);

    Vector3d pos  = Vector3d(arc.pos0(3*idEpoch+0,0),
                             arc.pos0(3*idEpoch+1,0),
                             arc.pos0(3*idEpoch+2,0));
    Vector3d vel  = Vector3d(arc.vel0(3*idEpoch+0,0),
                             arc.vel0(3*idEpoch+1,0),
                             arc.vel0(3*idEpoch+2,0));

    Matrix F(6, Alpha.columns());
    if(gravityCount)
      parameterGravity->gravity(arc.times.at(idEpoch), arc.rotEarth.at(idEpoch).rotate(pos), F.slice(3, idxGravity, 3, gravityCount));
    if(satCount || satArcCount)
    {
      parameterAcceleration->compute(satellite, arc.times.at(idEpoch), pos, vel, arc.rotSat.at(idEpoch), arc.rotEarth.at(idEpoch), ephemerides,
                                     F.slice(3, idxSat, 3, satCount), F.slice(3, idxSatArc, 3, satArcCount-6));
    }

    Matrix Z(6,6);
    matMult(1., arc.rotEarth.at(idEpoch).matrix(), arc.PosState.row(3*idEpoch,3), Z.row(0,3));
    matMult(1., arc.rotEarth.at(idEpoch).matrix(), arc.VelState.row(3*idEpoch,3), Z.row(3,3));
    solveInPlace(Z, F);

    return F;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/
