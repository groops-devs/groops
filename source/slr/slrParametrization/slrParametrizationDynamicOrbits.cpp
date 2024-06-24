/***********************************************/
/**
* @file slrParametrizationDynamicOrbits.cpp
*
* @brief Orbits by variational equations.
* @see SlrParametrization
*
* @author Torsten Mayer-Guerr
* @date 2022-04-28
*
*/
/***********************************************/

#include "base/import.h"
#include "config/config.h"
#include "files/fileInstrument.h"
#include "files/fileMatrix.h"
#include "files/fileParameterName.h"
#include "classes/ephemerides/ephemerides.h"
#include "classes/parametrizationAcceleration/parametrizationAcceleration.h"
#include "classes/timeSeries/timeSeries.h"
#include "classes/timeSeries/timeSeries.h"
#include "classes/platformSelector/platformSelector.h"
#include "misc/observation/variationalEquationFromFile.h"
#include "slr/slr.h"
#include "slr/slrParametrization/slrParametrizationDynamicOrbits.h"

/***********************************************/

SlrParametrizationDynamicOrbits::SlrParametrizationDynamicOrbits(Config &config)
{
  try
  {
    TimeSeriesPtr stochasticPulse;

    readConfig(config, "name",                        name,                        Config::OPTIONAL, "parameter.dynamicOrbits", "used for parameter selection");
    readConfig(config, "selectSatellites",            selectorSatellites,          Config::MUSTSET,  "",     "");
    readConfig(config, "outputfileOrbit",             fileNameOrbit,               Config::OPTIONAL, "",     "variable {satellite} available");
    readConfig(config, "outputfileParameters",        fileNameParameter,           Config::OPTIONAL, "",     "variable {satellite} available");
    readConfig(config, "inputfileVariational",        fileNameVariational,         Config::MUSTSET,  "variational_{loopTime:%D}.{satellite}.dat", "variable {satellite} available");
    readConfig(config, "stochasticPulse",             stochasticPulse,             Config::DEFAULT,  "",     "[mu/s] parametrization of stochastic pulses");
    readConfig(config, "parametrizationAcceleration", parametrizationAcceleration, Config::DEFAULT,  "",     "orbit force parameters");
    readConfig(config, "ephemerides",                 ephemerides,                 Config::MUSTSET,  "",     "");
    readConfig(config, "integrationDegree",           integrationDegree,           Config::DEFAULT,  "7",    "integration of forces by polynomial approximation of degree n");
    readConfig(config, "interpolationDegree",         interpolationDegree,         Config::DEFAULT,  "7",    "for orbit interpolation and velocity calculation");
    if(isCreateSchema(config)) return;

    pulses = stochasticPulse->times();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

SlrParametrizationDynamicOrbits::~SlrParametrizationDynamicOrbits()
{
  for(Parameter *para : parameters)
    delete para;
}

/***********************************************/

void SlrParametrizationDynamicOrbits::init(Slr *slr, const std::vector<const SlrParametrizationGravityField*> &paramGravityField)
{
  try
  {
    this->slr = slr;
    auto selectedSatellites = slr->selectSatellites(selectorSatellites);
    VariableList varList;
    varList.setVariable("satellite", "****");

    this->paramGravityField = paramGravityField;
    std::vector<ParametrizationGravityPtr> parameterGravity(paramGravityField.size());
    for(UInt i=0; i<paramGravityField.size(); i++)
      parameterGravity.at(i) = paramGravityField.at(i)->parametrization;

    parameters.resize(slr->satellites.size(), nullptr);
    for(UInt idSat=0; idSat<slr->satellites.size(); idSat++)
      if(selectedSatellites.at(idSat) && slr->satellites.at(idSat)->useable())
      {
        auto para = new Parameter();
        parameters.at(idSat) = para;
        para->sat = slr->satellites.at(idSat);

        varList.setVariable("satellite", slr->satellites.at(idSat)->name());
        FileVariationalEquation file(fileNameVariational(varList));
        if(file.arcCount() > 1)
          throw(Exception("<"+fileNameVariational(varList).str()+"> must not contain more than one arc"));
        VariationalEquationArc arc = file.readArc(0/*arcNo*/);

        VariationalEquation variationalEquation;
        variationalEquation.init(file.satellite(), parameterGravity, parametrizationAcceleration, pulses, ephemerides, integrationDegree);
        variationalEquation.setArc(arc);

        // parameter names
        variationalEquation.parameterNameSatellite(para->parameterNames);
        variationalEquation.parameterNameSatelliteArc(para->parameterNames);
        for(auto &name : para->parameterNames)
          name.object = para->sat->name();

        // integrate arc
        logStatus<<"integrate variational equations of "<<slr->satellites.at(idSat)->name()<<Log::endl;
        para->times     = arc.times;
        para->PosDesign = Matrix(3*para->times.size(), variationalEquation.parameterCount());
        para->VelDesign = Matrix(3*para->times.size(), variationalEquation.parameterCount());
        para->pos       = Vector(3*para->times.size());
        para->vel       = Vector(3*para->times.size());
        para->x         = Vector(variationalEquation.parameterCount());
        para->idxParameterSatellite = variationalEquation.parameterCountGravity();
        para->polynomial.init(para->times, interpolationDegree);
        // Log::Timer timer(para->times.size());
        for(UInt idEpoch=0; idEpoch<para->times.size(); idEpoch++)
        {
          // timer.loopStep(idEpoch);
          variationalEquation.position(idEpoch, para->pos.row(3*idEpoch,3), para->PosDesign.row(3*idEpoch,3));
          variationalEquation.velocity(idEpoch, para->vel.row(3*idEpoch,3), para->VelDesign.row(3*idEpoch,3));
        } // for(idEpoch)
        // timer.loopEnd();

        // set satellites position and velocity to arc values
        reshape(para->polynomial.interpolate(para->sat->times, para->pos, 3), para->sat->pos.trans());
        reshape(para->polynomial.interpolate(para->sat->times, para->vel, 3), para->sat->vel.trans());
      } // for(idSat)
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void SlrParametrizationDynamicOrbits::initParameter(SlrNormalEquationInfo &normalEquationInfo)
{
  try
  {
   for(auto para : parameters)
      if(para)
        para->index = SlrParameterIndex();
    if(!isEnabled(normalEquationInfo, name))
      return;

    UInt countPara = 0;
    for(auto para : parameters)
      if(para && para->sat->useable() && normalEquationInfo.estimateSatellite.at(para->sat->idSat()))
      {
        para->index = normalEquationInfo.parameterNamesSatellite(para->sat->idSat(), para->parameterNames);
        countPara += para->parameterNames.size();
      }
    if(countPara)
      logInfo<<countPara%"%9i dynamic orbit parameters"s<<Log::endl;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void SlrParametrizationDynamicOrbits::aprioriParameter(const SlrNormalEquationInfo &normalEquationInfo, MatrixSliceRef x0) const
{
  try
  {
    for(auto para : parameters)
      if(para && para->index)
        copy(para->x.row(para->idxParameterSatellite, normalEquationInfo.count(para->index)),
             x0.row(normalEquationInfo.index(para->index), normalEquationInfo.count(para->index)));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void SlrParametrizationDynamicOrbits::designMatrix(const SlrNormalEquationInfo &/*normalEquationInfo*/, const SlrObservationEquation &eqn, SlrDesignMatrix &A) const
{
  try
  {
    auto para = parameters.at(eqn.satellite->idSat());
    if(!para || !para->index)
      return;

    const Matrix PosDesign = para->polynomial.interpolate(eqn.timesBounce, para->PosDesign, 3);

    // gravity
    UInt idx = 0;
    for(auto paraGravity : paramGravityField)
    {
      if(paraGravity->indexParameter)
      {
        MatrixSlice Design(A.column(paraGravity->indexParameter));
        for(UInt idEpoch=0; idEpoch<eqn.timesBounce.size(); idEpoch++)
          matMult(1., eqn.A.slice(idEpoch, SlrObservationEquation::idxPosSat, 1, 3),
                  PosDesign.slice(3*idEpoch, idx, 3, Design.columns()),
                  Design.row(idEpoch));
      }
      idx += paraGravity->parametrization->parameterCount();
    }

    // satellite related parameters, includes satellite state (6 parameter)
    MatrixSlice Design(A.column(para->index));
    for(UInt idEpoch=0; idEpoch<eqn.timesBounce.size(); idEpoch++)
      matMult(1., eqn.A.slice(idEpoch, SlrObservationEquation::idxPosSat, 1, 3),
              PosDesign.slice(3*idEpoch, para->idxParameterSatellite, 3, Design.columns()),
              Design.row(idEpoch));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Double SlrParametrizationDynamicOrbits::updateParameter(const SlrNormalEquationInfo &normalEquationInfo, const_MatrixSliceRef x, const_MatrixSliceRef /*Wz*/)
{
  try
  {
    Double maxChange = 0;
    Slr::InfoParameterChange info("mm");
    for(auto para : parameters)
      if(para && para->index)
      {
        Vector dx(para->x.rows());

        // gravity
        UInt idx = 0;
        for(auto paraGravity : paramGravityField)
        {
          const UInt count = paraGravity->parametrization->parameterCount();
          if(paraGravity->indexParameter)
            copy(x.row(normalEquationInfo.index(paraGravity->indexParameter), count), dx.row(idx, count));
          idx += count;
        }

        // satellite related parameters, includes satellite state (6 parameter)
        copy(x.row(normalEquationInfo.index(para->index), normalEquationInfo.count(para->index)),
             dx.row(para->idxParameterSatellite, normalEquationInfo.count(para->index)));

        para->x   += dx;
        para->pos += para->PosDesign * dx;
        para->vel += para->VelDesign * dx;

        auto sat = para->sat;
        Matrix posOld = sat->pos;

        reshape(para->polynomial.interpolate(sat->times, para->pos, 3), sat->pos.trans());
        reshape(para->polynomial.interpolate(sat->times, para->vel, 3), sat->vel.trans());

        for(UInt idEpoch=0; idEpoch<sat->times.size(); idEpoch++)
        {
          const Double dr = norm(sat->pos.row(idEpoch)-posOld.row(idEpoch));
          if(info.update(1e3*dr))
            info.info = "position satellite ("+sat->name()+", "+sat->times.at(idEpoch).dateTimeStr()+")";
        }
      }
    info.print(1e-3, maxChange);
    return maxChange;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void SlrParametrizationDynamicOrbits::writeResults(const SlrNormalEquationInfo &normalEquationInfo, const std::string &suffix) const
{
  try
  {
    if(!isEnabled(normalEquationInfo, name))
      return;

    if(!fileNameOrbit.empty())
    {
      VariableList varList;
      varList.setVariable("satellite", "***");
      logStatus<<"write satellite orbits to files <"<<fileNameOrbit(varList).appendBaseName(suffix)<<">"<<Log::endl;
      for(auto para : parameters)
        if(para && para->index)
        {

          OrbitArc arc;
          for(UInt idEpoch=0; idEpoch<para->times.size(); idEpoch++)
          {
            OrbitEpoch epoch;
            epoch.time     = para->times.at(idEpoch);
            epoch.position = Vector3d(para->pos.row(3*idEpoch, 3));
            epoch.velocity = Vector3d(para->vel.row(3*idEpoch, 3));
            arc.push_back(epoch);
          }
          varList.setVariable("satellite", para->sat->name());
          InstrumentFile::write(fileNameOrbit(varList).appendBaseName(suffix), arc);
        }
    }

    if(!fileNameParameter.empty())
    {
      VariableList varList;
      varList.setVariable("satellite", "***");
      logStatus<<"write estimated satellite parameters to files <"<<fileNameParameter(varList).appendBaseName(suffix)<<">"<<Log::endl;
      for(auto para : parameters)
        if(para && para->index)
        {
          varList.setVariable("satellite", para->sat->name());
          writeFileMatrix(fileNameParameter(varList).appendBaseName(suffix), para->x);
        }
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
