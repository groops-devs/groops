/***********************************************/
/**
* @file variationalEquationFromFile.cpp
*
* @brief Variational equations.
*
* @author Torsten Mayer-Guerr
* @date 2014-03-24
*
*/
/***********************************************/

#include "base/import.h"
#include "files/fileVariationalEquation.h"
#include "classes/parametrizationGravity/parametrizationGravity.h"
#include "classes/parametrizationAcceleration/parametrizationAcceleration.h"
#include "misc/observation/variationalEquation.h"
#include "misc/observation/variationalEquationFromFile.h"

/***********************************************/

void VariationalEquationFromFile::open(const FileName &fileName, ParametrizationGravityPtr parameterGravity, ParametrizationAccelerationPtr parameterAcceleration,
                                       const std::vector<Time> &stochasticPulse, EphemeridesPtr ephemerides, UInt integrationDegree)
{
  try
  {
    file.open(fileName);

    variationalEquation.init(file.satellite(), parameterGravity, parameterAcceleration, stochasticPulse, ephemerides, integrationDegree);
    arcNo = 0;
    VariationalEquationArc arc = file.readArc(arcNo);
    variationalEquation.setArc(arc);
    computeIndices();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void VariationalEquationFromFile::close()
{
  try
  {
    file.close();
    _parameterCount = 0;
    gravityCount = 0;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void VariationalEquationFromFile::computeIndices()
{
  try
  {
    variationalEquation.computeIndices();

    _parameterCount  = 0;

    idxGravity       = _parameterCount;
    gravityCount     = variationalEquation.parameterCountGravity();
    _parameterCount += gravityCount;

    idxSat           = _parameterCount;
    satCount         = variationalEquation.parameterCountSatellite();
    _parameterCount += satCount;

    idxSatArc        = _parameterCount;
    satArcCount      = variationalEquation.parameterCountSatelliteArc();
    _parameterCount += file.arcCount() * satArcCount;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}


/***********************************************/

void VariationalEquationFromFile::parameterName(std::vector<ParameterName> &name) const
{
  try
  {
    parameterNameGravity(name);
    parameterNameSatellite(name);
    parameterNameSatelliteArc(name);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void VariationalEquationFromFile::parameterNameGravity(std::vector<ParameterName> &name) const
{
  try
  {
    variationalEquation.parameterNameGravity(name);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void VariationalEquationFromFile::parameterNameSatellite(std::vector<ParameterName> &name) const
{
  try
  {
    variationalEquation.parameterNameSatellite(name);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void VariationalEquationFromFile::parameterNameSatelliteArc(std::vector<ParameterName> &name) const
{
  try
  {
    std::vector<ParameterName> nameArc;
    variationalEquation.parameterNameSatelliteArc(nameArc);
    for(UInt arcNo=0; arcNo<arcCount(); arcNo++)
    {
      const std::string str = "arc"+arcNo%"%i"s+".";
      for(UInt i=0; i<nameArc.size(); i++)
      {
        ParameterName param = nameArc.at(i);
        param.type = str+param.type;
        name.push_back(param);
      }
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

VariationalEquationFromFile::ObservationEquation VariationalEquationFromFile::integrateArc(Time timeStart, Time timeEnd, Bool computePosition, Bool computeVelocity)
{
  try
  {
    const VariationalEquationArc &arc = getArc(timeStart);
    if(timeEnd > arc.times.back())
    {
      timeEnd = arc.times.back();
//       throw(Exception("cannot integrate over different variational arcs ("+timeStart.dateTimeStr()+" - "+timeEnd.dateTimeStr()+") vs. ("+arc.times.at(0).dateTimeStr()+" - "+arc.times.back().dateTimeStr()+")"));
    }

    // find epoch interval
    const UInt epochStart = std::distance(arc.times.begin(), std::upper_bound(arc.times.begin(), arc.times.end(), timeStart))-1;
    const UInt epochEnd   = std::distance(arc.times.begin(), std::lower_bound(arc.times.begin(), arc.times.end(), timeEnd));
    ObservationEquation eqn;
    for(UInt idEpoch=epochStart; idEpoch<=epochEnd; idEpoch++)
      eqn.times.push_back(arc.times.at(idEpoch));

    if(computePosition)
    {
      eqn.pos0      = Vector(3*eqn.times.size());
      eqn.PosDesign = Matrix(3*eqn.times.size(), parameterCount());
    }
    if(computeVelocity)
    {
      eqn.vel0      = Vector(3*eqn.times.size());
      eqn.VelDesign = Matrix(3*eqn.times.size(), parameterCount());
    }
    eqn.rotSat.resize(eqn.times.size());
    eqn.rotEarth.resize(eqn.times.size());

    for(UInt i=0; i<eqn.times.size(); i++)
    {
      const UInt idEpoch = epochStart+i;

      if(computePosition)
      {
        Matrix PosDesign(3, variationalEquation.parameterCount());
        variationalEquation.position(idEpoch, eqn.pos0.row(3*i,3), PosDesign);
        if(gravityCount+satCount)
          copy(PosDesign.column(0, gravityCount+satCount), eqn.PosDesign.slice(3*i,0,3,gravityCount+satCount));
        if(satArcCount)
          copy(PosDesign.column(idxSatArc, satArcCount), eqn.PosDesign.slice(3*i,idxSatArc+arcNo*satArcCount,3,satArcCount));
      }

      if(computeVelocity)
      {
        Matrix VelDesign(3, variationalEquation.parameterCount());
        variationalEquation.velocity(idEpoch, eqn.vel0.row(3*i,3), VelDesign);
        if(gravityCount+satCount)
          copy(VelDesign.column(0, gravityCount+satCount), eqn.VelDesign.slice(3*i,0,3,gravityCount+satCount));
        if(satArcCount)
          copy(VelDesign.column(idxSatArc, satArcCount), eqn.VelDesign.slice(3*i,idxSatArc+arcNo*satArcCount,3,satArcCount));
      }

      eqn.rotSat.at(i)   = arc.rotSat.at(idEpoch);
      eqn.rotEarth.at(i) = arc.rotEarth.at(idEpoch);
    } // for(i=times)

    return eqn;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

const VariationalEquationArc &VariationalEquationFromFile::getArc(const Time &time)
{
  try
  {
    if((variationalEquation.arc().times.front() <= time) && (time <= variationalEquation.arc().times.back()))
      return variationalEquation.arc();

    // restart file?
    VariationalEquationArc arc = variationalEquation.arc();
    if(time < arc.times.front())
      arc = file.readArc(arcNo = 0);

    // find arc
    while((time > arc.times.back()) && (++arcNo < file.arcCount()))
      arc = file.readArc(arcNo);

    if((time > arc.times.back()) || (time < arc.times.front()))
      throw(Exception("no variational arc found for "+time.dateTimeStr()));

    const UInt countOld = variationalEquation.parameterCountSatellite();
    variationalEquation.setArc(arc);
    if(variationalEquation.parameterCountSatellite() != countOld)
      throw(Exception("size of arc related parameters changed from arc to arc."));

    return variationalEquation.arc();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

VariationalEquationArc VariationalEquationFromFile::refineVariationalEquationArc(UInt arcNo, const_MatrixSliceRef x)
{
  try
  {
    const UInt countOld = variationalEquation.parameterCountSatellite();
    VariationalEquationArc arc = file.readArc(arcNo);
    variationalEquation.setArc(arc);
    if(variationalEquation.parameterCountSatellite() != countOld)
      throw(Exception("size of arc related parameters changed from arc to arc."));

    for(UInt idEpoch=0; idEpoch<arc.times.size(); idEpoch++)
    {
      Vector pos0(3);
      Matrix PosDesign(3, variationalEquation.parameterCount());
      variationalEquation.position(idEpoch, pos0, PosDesign);

      Vector vel0(3);
      Matrix VelDesign(3, variationalEquation.parameterCount());
      variationalEquation.velocity(idEpoch, vel0, VelDesign);

      if(gravityCount+satCount)
      {
        matMult(1, PosDesign.column(0, gravityCount+satCount), x.slice(0, 0, gravityCount+satCount, 1), arc.pos0.row(3*idEpoch,3));
        matMult(1, VelDesign.column(0, gravityCount+satCount), x.slice(0, 0, gravityCount+satCount, 1), arc.vel0.row(3*idEpoch,3));
      }

      if(satArcCount)
      {
        matMult(1, PosDesign.column(idxSatArc, satArcCount), x.slice(idxSatArc+arcNo*satArcCount, 0, satArcCount, 1), arc.pos0.row(3*idEpoch,3));
        matMult(1, VelDesign.column(idxSatArc, satArcCount), x.slice(idxSatArc+arcNo*satArcCount, 0, satArcCount, 1), arc.vel0.row(3*idEpoch,3));
      }
    } // for(idEpoch)

    return arc;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/
