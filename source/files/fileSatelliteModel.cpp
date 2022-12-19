/***********************************************/
/**
* @file fileSatelliteModel.cpp
*
* @brief Satellite model.
*
* @author Torsten Mayer-Guerr
* @author Sebastian Strasser
* @date 2015-05-25
*
*/
/***********************************************/

#define DOCSTRING_FILEFORMAT_SatelliteModel

#include "base/import.h"
#include "inputOutput/fileArchive.h"
#include "files/fileFormatRegister.h"
#include "files/fileSatelliteModel.h"

GROOPS_REGISTER_FILEFORMAT(SatelliteModel, FILE_SATELLITEMODEL_TYPE)

/***********************************************/

SatelliteModelModulePtr SatelliteModelModule::create(SatelliteModelModule::Type type)
{
  try
  {
    switch(type)
    {
      case SatelliteModelModule::SOLARPANEL:           return SatelliteModelModulePtr(new SatelliteModelModuleSolarPanel());
      case SatelliteModelModule::ANTENNATHRUST:        return SatelliteModelModulePtr(new SatelliteModelModuleAntennaThrust());
      case SatelliteModelModule::MASSCHANGE:           return SatelliteModelModulePtr(new SatelliteModelModuleMassChange());
      case SatelliteModelModule::SPECIFICHEATCAPACITY: return SatelliteModelModulePtr(new SatelliteModelModuleSetSpecificHeatCapacity());
    }

    std::stringstream ss;
    ss<<"unknown satellite module type ("<<type<<")";
    throw(Exception(ss.str()));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/

void SatelliteModelModuleSolarPanel::changeState(SatelliteModel &satellite, const Time &/*time*/,
                                                 const Vector3d &position, const Vector3d &/*velocity*/,
                                                 const Vector3d &positionSun, const Rotary3d &rotSat, const Rotary3d &/*rotEarth*/)
{
  try
  {
    if(positionSun.r() == 0.)
      throw(Exception("Need position of the sun"));
    const Vector3d posSun = rotSat.inverseRotate(positionSun-position); // in SRF
    const Rotary3d rot = Rotary3d(rotationAxis, posSun)
                       * inverse(Rotary3d(rotationAxis, normal)); // from old to new orientation
    normal = rot.rotate(normal);
    for(UInt i=0; i<indexSurface.size(); i++)
      satellite.surfaces.at(indexSurface.at(i)).normal = rot.rotate(satellite.surfaces.at(indexSurface.at(i)).normal);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void SatelliteModelModuleSolarPanel::save(OutArchive &ar) const
{
  try
  {
    ar<<nameValue("rotationAxis", rotationAxis);
    ar<<nameValue("normal",       normal);
    ar<<nameValue("surface",      indexSurface);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void SatelliteModelModuleSolarPanel::load(InArchive  &ar)
{
  try
  {
    ar>>nameValue("rotationAxis", rotationAxis);
    ar>>nameValue("normal",       normal);
    ar>>nameValue("surface",      indexSurface);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/

void SatelliteModelModuleAntennaThrust::changeState(SatelliteModel &satellite, const Time &/*time*/,
                                                    const Vector3d &/*position*/, const Vector3d &/*velocity*/,
                                                    const Vector3d &/*positionSun*/, const Rotary3d &/*rotSat*/, const Rotary3d &/*rotEarth*/)
{
  if(!applied)
    satellite.thrustPower += thrust;
  applied = TRUE;
}

/***********************************************/

void SatelliteModelModuleAntennaThrust::save(OutArchive &ar) const
{
  ar<<nameValue("antennaThrust", thrust);
}

/***********************************************/

void SatelliteModelModuleAntennaThrust::load(InArchive  &ar)
{
  ar>>nameValue("antennaThrust", thrust);
  applied = FALSE;
}

/***********************************************/
/***********************************************/

void SatelliteModelModuleMassChange::changeState(SatelliteModel &satellite, const Time &time,
                                                 const Vector3d &/*position*/, const Vector3d &/*velocity*/,
                                                 const Vector3d &/*positionSun*/, const Rotary3d &/*rotSat*/, const Rotary3d &/*rotEarth*/)
{
  try
  {
    if((times.size()==0)||(time<times.at(0)))
      return;

    if(time>=times.back())
    {
      satellite.mass = mass.back();
      return;
    }

    // linear interpolation
    UInt idx = 0;
    while((idx+1<times.size()) && (times.at(idx+1)<=time))
      idx++;
    const Double tau = (time-times.at(idx)).mjd()/(times.at(idx+1)-times.at(idx)).mjd();
    satellite.mass = (1-tau)*mass.at(idx) + tau*mass.at(idx+1);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void SatelliteModelModuleMassChange::save(OutArchive &ar) const
{
  try
  {
    ar<<nameValue("epochCount", times.size());
    for(UInt i=0; i<times.size(); i++)
    {
      ar<<beginGroup("epoch");
      ar<<nameValue("time", times.at(i));
      ar<<nameValue("mass", mass.at(i));
      ar<<endGroup("epoch");
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void SatelliteModelModuleMassChange::load(InArchive &ar)
{
  try
  {
    UInt count;
    ar>>nameValue("epochCount", count);
    times.resize(count);
    mass.resize(count);
    for(UInt i=0; i<times.size(); i++)
    {
      ar>>beginGroup("epoch");
      ar>>nameValue("time", times.at(i));
      ar>>nameValue("mass", mass.at(i));
      ar>>endGroup("epoch");
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/

void SatelliteModelModuleSetSpecificHeatCapacity::changeState(SatelliteModel &satellite, const Time &/*time*/,
                                                              const Vector3d &/*position*/, const Vector3d &/*velocity*/,
                                                              const Vector3d &/*positionSun*/, const Rotary3d &/*rotSat*/, const Rotary3d &/*rotEarth*/)
{
  try
  {
    for(UInt i=0; i<indexSurface.size(); i++)
      satellite.surfaces.at(indexSurface.at(i)).specificHeatCapacity = specificHeatCapacity.at(i);
    applied = TRUE;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void SatelliteModelModuleSetSpecificHeatCapacity::save(OutArchive &ar) const
{
  try
  {
    ar<<nameValue("specificHeatCapacity", specificHeatCapacity);
    ar<<nameValue("surface",              indexSurface);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void SatelliteModelModuleSetSpecificHeatCapacity::load(InArchive  &ar)
{
  try
  {
    ar>>nameValue("specificHeatCapacity", specificHeatCapacity);
    ar>>nameValue("surface",              indexSurface);
    applied = FALSE;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/

SatelliteModel::Surface::Surface() :
  type(PLATE), area(0),
  absorptionVisible(0), absorptionInfrared(0),
  diffusionVisible(0),  diffusionInfrared(0),
  reflexionVisible(0),  reflexionInfrared(0),
  specificHeatCapacity(0)
{
}

/***********************************************/

template<> void save(OutArchive &ar, const SatelliteModel::Surface &x)
{
  try
  {
    ar<<nameValue("type",                 static_cast<UInt>(x.type));
    ar<<nameValue("normal",               x.normal);
    ar<<nameValue("area",                 x.area);
    ar<<nameValue("reflexionVisible",     x.reflexionVisible);
    ar<<nameValue("diffusionVisible",     x.diffusionVisible);
    ar<<nameValue("absorptionVisible",    x.absorptionVisible);
    ar<<nameValue("reflexionInfrared",    x.reflexionInfrared);
    ar<<nameValue("diffusionInfrared",    x.diffusionInfrared);
    ar<<nameValue("absorptionInfrared",   x.absorptionInfrared);
    ar<<nameValue("hasThermalReemission", (x.specificHeatCapacity != 0));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

template<> void load(InArchive  &ar, SatelliteModel::Surface &x)
{
  try
  {
    UInt type;
    ar>>nameValue("type",               type);
    x.type = static_cast<SatelliteModel::Surface::Type>(type);
    ar>>nameValue("normal",             x.normal);
    ar>>nameValue("area",               x.area);
    ar>>nameValue("reflexionVisible",   x.reflexionVisible);
    ar>>nameValue("diffusionVisible",   x.diffusionVisible);
    ar>>nameValue("absorptionVisible",  x.absorptionVisible);
    ar>>nameValue("reflexionInfrared",  x.reflexionInfrared);
    ar>>nameValue("diffusionInfrared",  x.diffusionInfrared);
    ar>>nameValue("absorptionInfrared", x.absorptionInfrared);

    x.specificHeatCapacity = 0;
    if(ar.version() >= 20190429)
    {
      Bool hasThermalReemission;
      ar>>nameValue("hasThermalReemission", hasThermalReemission);
      x.specificHeatCapacity = (hasThermalReemission ? -1. : 0.);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/

SatelliteModel::SatelliteModel() : satelliteName("NULL"), mass(0), coefficientDrag(0)
{
}

/***********************************************/

SatelliteModel::~SatelliteModel()
{
}

/***********************************************/

void SatelliteModel::changeState(const Time &time, const Vector3d &position, const Vector3d &velocity,
                                 const Vector3d &positionSun, const Rotary3d &rotSat, const Rotary3d &rotEarth)
{
  for(auto module : modules)
    module ->changeState(*this, time, position, velocity, positionSun, rotSat, rotEarth);
}


/***********************************************/
/***********************************************/

template<> void save(OutArchive &ar, const SatelliteModelPtr &x)
{
  try
  {
    if((x==nullptr) || (x->satelliteName=="NULL"))
    {
      ar<<nameValue("satelliteName", std::string("NULL"));
      return;
    }

    ar<<nameValue("satelliteName",   x->satelliteName);
    ar<<nameValue("mass",            x->mass);
    ar<<nameValue("coefficientDrag", x->coefficientDrag);
    ar<<nameValue("surfaceCount",    x->surfaces.size());
    ar.endLine();
    for(UInt i=0; i<x->surfaces.size(); i++)
    {
      ar<<nameValue("surface", x->surfaces.at(i));
      ar.endLine();
    }
    ar<<nameValue("modulCount", x->modules.size());
    for(UInt i=0; i<x->modules.size(); i++)
    {
      ar<<beginGroup("modul");
      ar<<nameValue("type", static_cast<UInt>(x->modules.at(i)->type()));
      x->modules.at(i)->save(ar);
      ar<<endGroup("modul");
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

template<> void load(InArchive &ar, SatelliteModelPtr &x)
{
  try
  {
    std::string satelliteName;
    ar>>nameValue("satelliteName", satelliteName);
    if(satelliteName=="NULL")
    {
      x = SatelliteModelPtr(nullptr);
      return;
    }

    x = SatelliteModelPtr(new SatelliteModel);
    x->satelliteName = satelliteName;

    UInt count, type;
    ar>>nameValue("mass",            x->mass);
    ar>>nameValue("coefficientDrag", x->coefficientDrag);
    ar>>nameValue("surfaceCount",    count);
    x->surfaces.resize(count);
    for(UInt i=0; i<x->surfaces.size(); i++)
      ar>>nameValue("surface", x->surfaces.at(i));
    ar>>nameValue("modulCount",      count);
    x->modules.resize(count);
    for(UInt i=0; i<x->modules.size(); i++)
    {
      ar>>beginGroup("modul");
      ar>>nameValue("type", type);
      x->modules.at(i) = SatelliteModelModule::create(static_cast<SatelliteModelModule::Type>(type));
      x->modules.at(i)->load(ar);
      ar>>endGroup("modul");
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void writeFileSatelliteModel(const FileName &fileName, const SatelliteModelPtr &x)
{
  try
  {
    OutFileArchive file(fileName, FILE_SATELLITEMODEL_TYPE);
    file<<nameValue("satelliteCount", 1);
    file<<nameValue("satellite", x);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void writeFileSatelliteModel(const FileName &fileName, const std::vector<SatelliteModelPtr> &x)
{
  try
  {
    OutFileArchive file(fileName, FILE_SATELLITEMODEL_TYPE);
    file<<nameValue("satelliteCount", x.size());
    for(UInt i=0; i<x.size(); i++)
      file<<nameValue("satellite", x.at(i));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void readFileSatelliteModel(const FileName &fileName, SatelliteModelPtr &x)
{
  try
  {
    InFileArchive file(fileName, FILE_SATELLITEMODEL_TYPE);
    UInt count;
    file>>nameValue("satelliteCount", count);
    if(count>1)
      logWarning<<fileName<<" contain more than one satellite, only the first is used"<<Log::endl;
    file>>nameValue("satellite", x);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void readFileSatelliteModel(const FileName &fileName, std::vector<SatelliteModelPtr> &x)
{
  try
  {
    InFileArchive file(fileName, FILE_SATELLITEMODEL_TYPE);
    UInt count;
    file>>nameValue("satelliteCount", count);
    x.resize(count);
    for(UInt i=0; i<x.size(); i++)
      file>>nameValue("satellite", x.at(i));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/
