/***********************************************/
/**
* @file ephemerides.cpp
*
* @brief Ephemerides of Sun, Moon and planets.
*
* @author Torsten Mayer-Guerr
* @date 2020-06-07
*
*/
/***********************************************/

#define DOCSTRING_Ephemerides
#define DOCSTRING_EphemeridesPlanet

#include "base/import.h"
#include "config/configRegister.h"
#include "classes/ephemerides/ephemeridesJpl.h"
#include "classes/ephemerides/ephemerides.h"

/***********************************************/

GROOPS_REGISTER_CLASS(Ephemerides, "ephemeridesType",
                      EphemeridesJpl)

GROOPS_READCONFIG_CLASS(Ephemerides, "ephemeridesType")

/***********************************************/

EphemeridesPtr Ephemerides::create(Config &config, const std::string &name)
{
  try
  {
    EphemeridesPtr ptr;
    std::string    type;

    readConfigChoice(config, name, type, Config::MUSTSET, "", "ephemerides of Sun, Moon and planets");
    if(readConfigChoiceElement(config, "jpl", type, ""))
      ptr = std::make_shared<EphemeridesJpl>(config);
    endChoice(config);

    return ptr;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***** CLASS ***********************************/
/***********************************************/

// Wrapper class
class EphemeridesPlanet
{
public:
  Ephemerides::Planet planet;
  EphemeridesPlanet(Config &config, const std::string &name);
  static EphemeridesPlanet create(Config &config, const std::string &name) {return EphemeridesPlanet(config, name);}
};

GROOPS_REGISTER_CLASS_WITHOUT_SUBS(EphemeridesPlanet, "planetType")

/***********************************************/

template<> Bool readConfig(Config &config, const std::string &name, Ephemerides::Planet &var, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation)
{
  try
  {
    if(isCreateSchema(config))
    {
      config.xselement(name, "planetType", mustSet, Config::ONCE, defaultValue, annotation);
      return FALSE;
    }

    if(!hasName(config, name, mustSet))
      return FALSE;
    EphemeridesPlanet tmp(config, name);
    var = tmp.planet;
    return TRUE;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

EphemeridesPlanet::EphemeridesPlanet(Config &config, const std::string &name)
{
  try
  {
    planet = Ephemerides::EARTH;
    std::string choice;

    readConfigChoice(config, name, choice, Config::MUSTSET, "", "planet");
    if(readConfigChoiceElement(config, "earth",               choice, "")) planet = Ephemerides::EARTH;
    if(readConfigChoiceElement(config, "sun",                 choice, "")) planet = Ephemerides::SUN;
    if(readConfigChoiceElement(config, "moon",                choice, "")) planet = Ephemerides::MOON;
    if(readConfigChoiceElement(config, "mercury",             choice, "")) planet = Ephemerides::MERCURY;
    if(readConfigChoiceElement(config, "venus",               choice, "")) planet = Ephemerides::VENUS;
    if(readConfigChoiceElement(config, "mars",                choice, "")) planet = Ephemerides::MARS;
    if(readConfigChoiceElement(config, "jupiter",             choice, "")) planet = Ephemerides::JUPITER;
    if(readConfigChoiceElement(config, "saturn",              choice, "")) planet = Ephemerides::SATURN;
    if(readConfigChoiceElement(config, "uranus",              choice, "")) planet = Ephemerides::URANUS;
    if(readConfigChoiceElement(config, "neptune",             choice, "")) planet = Ephemerides::NEPTUNE;
    if(readConfigChoiceElement(config, "pluto",               choice, "")) planet = Ephemerides::PLUTO;
    if(readConfigChoiceElement(config, "solarBaryCenter",     choice, "")) planet = Ephemerides::SOLARBARYCENTER;
    if(readConfigChoiceElement(config, "earthMoonBaryCenter", choice, "")) planet = Ephemerides::EARTHMOONBARYCENTER;
    endChoice(config);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
