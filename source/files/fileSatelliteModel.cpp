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
      case SatelliteModelModule::SOLARPANEL:    return SatelliteModelModulePtr(new SatelliteModelModuleSolarPanel());
      case SatelliteModelModule::ANTENNATHRUST: return SatelliteModelModulePtr(new SatelliteModelModuleAntennaThrust());
      case SatelliteModelModule::MASSCHANGE:    return SatelliteModelModulePtr(new SatelliteModelModuleMassChange());
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

void SatelliteModelModuleAntennaThrust::accelerationThrust(const SatelliteModel &/*satellite*/, Vector3d &a) const
{
  a += thrust;
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

void SatelliteModelModuleMassChange::load(InArchive  &ar)
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

Double SatelliteModel::Surface::sectionalArea(const Vector3d &direction) const
{
  try
  {
    if(type == PLATE)
    {
      const Double cosPhi = inner(direction, normal);
      return (cosPhi>0) ? (area*cosPhi) : 0;
    }
    else if(type == SPHERE)
      return area;
    else if(type == CYLINDER)
      return area*crossProduct(direction, normal).norm();

    throw(Exception("unknown type"));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}


/***********************************************/

Vector3d SatelliteModel::Surface::accelerationDrag(const Vector3d &direction, Double velocity, Double density, Double temperature) const
{
  try
  {
    constexpr Double R     = 8.3144;    // universal gas constant [J mol-1 K-1]
    constexpr Double M     = 15.5/1000; // average modelcular mass [kg/mol]
    constexpr Double alpha = 0.9;       // accomodation coefficientDrag [-]
    constexpr Double Tw    = 273;       // temperature of the surface [K]

    if(type == PLATE)
    {
      const Double cosPhi = -inner(direction, normal);

      const Vector3d ul   = crossProduct(direction, crossProduct(direction, normal));
      const Double   s    = velocity/std::sqrt(2*R*temperature/M); // molecular speed ratio
      const Double   P    = std::exp(-std::pow(cosPhi*s, 2))/s;
      const Double   G    = 1/(2*s*s);
      const Double   Q    = 1 + G;
      const Double   Z    = 1 + std::erf(cosPhi*s);
      const Double   x    = std::sqrt(1./6.*(1 + alpha * (3*R*Tw/(velocity*velocity)-1))) * (cosPhi*std::sqrt(PI)*Z+P);
      return 0.5*density*area*velocity*velocity * ((P/std::sqrt(PI) + cosPhi*(Q*Z+x)) * direction + (G*Z+x) * ul);
    }
    else if(type == SPHERE)
    {
      throw(Exception("sphere not implemented yet"));
    }
    else if(type == CYLINDER)
    {
      throw(Exception("cylinder not implemented yet"));
    }

    throw(Exception("unknown type"));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void SatelliteModel::Surface::accelerationPressure(const Vector3d &direction, Double visible, Double infrared, Vector3d &acc) const
{
  try
  {
    const Double cosPhi = -inner(direction, normal);
    if(cosPhi<=0)
      return;

    // absorptionX + diffusionX + reflexionX = 1.
    const Double reflexion  = reflexionVisible  * visible + reflexionInfrared  * infrared;
    const Double diffusion  = diffusionVisible  * visible + diffusionInfrared  * infrared;
    const Double absorption = absorptionVisible * visible + absorptionInfrared * infrared;

    if(type == PLATE)
    {
      // Source: Montenbruck et al. (2015) Enhanced solar radiation pressure modeling for Galileo satellites. DOI 10.1007/s00190-014-0774-0
      // ATTENTION: our direction vector is defined with opposite sign compared to source => minuses for normal vector parts
      acc += area*cosPhi*(absorption+diffusion)*direction - area*cosPhi*(2./3.*diffusion+2.*cosPhi*reflexion)*normal;
      if(hasThermalReemission)
        acc -= area*cosPhi*(2./3.*absorption)*normal;
      return;
    }
    else if(type == SPHERE)
    {
      throw(Exception("sphere not implemented yet"));
    }
    else if(type == CYLINDER)
    {
      // Source: Rodriguez Solano (2014) Impact of non-conservative force modeling on GNSS satellite orbits and global solutions. PhD thesis
      // ATTENTION: our direction vector is defined with opposite sign compared to source => minuses for normal vector parts
      acc += area*cosPhi*(absorption+diffusion)*direction - area*cosPhi*(PI/6.*diffusion+4./3.*cosPhi*reflexion)*normal;
      if(hasThermalReemission)
        acc -= area*cosPhi*(PI/6.*absorption)*normal;
      return;
    }

    throw(Exception("unknown type"));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
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
    ar<<nameValue("hasThermalReemission", x.hasThermalReemission);
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
    x.hasThermalReemission = FALSE;
    if(ar.version() >= 20190429)
      ar>>nameValue("hasThermalReemission", x.hasThermalReemission);
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

Vector3d SatelliteModel::accelerationDrag(const Vector3d &velocity, Double density, Double temperature) const
{
  try
  {
    if(mass == 0.)
      throw(Exception("No SatelliteModel given: "+satelliteName));

    Vector3d direction = -velocity;
    const Double v = direction.normalize();
    Vector3d a;
    if(temperature > 0)
    {
      for(auto &surface : surfaces)
        a += surface.accelerationDrag(direction, v, density, temperature);
    }
    else
    {
      Double area = 0;
      for(auto &surface : surfaces)
        area += surface.sectionalArea(direction);
      a = 0.5*coefficientDrag*area*density*v*v*direction;
    }

    for(auto module : modules)
      module->accelerationDrag(*this, velocity, density, temperature, a);

    return (1./mass) * a;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Vector3d SatelliteModel::accelerationPressure(const Vector3d &direction, Double visible, Double infrared) const
{
  try
  {
    if(mass == 0.)
      throw(Exception("No SatelliteModel given: "+satelliteName));

    Vector3d a;
    for(UInt i=0; i<surfaces.size(); i++)
      surfaces.at(i).accelerationPressure(direction, visible, infrared, a);

    for(auto module : modules)
      module->accelerationPressure(*this, direction, visible, infrared, a);

    return (1./mass) * a;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Vector3d SatelliteModel::accelerationThrust() const
{
  try
  {
    if(mass == 0.)
      throw(Exception("No SatelliteModel given: "+satelliteName));

    Vector3d a;
    for(auto module : modules)
      module->accelerationThrust(*this, a);

    return (-1./mass/LIGHT_VELOCITY) * a;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
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
