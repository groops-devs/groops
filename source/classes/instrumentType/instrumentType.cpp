/***********************************************/
/**
* @file instrumentType.cpp
*
* @brief Defines the type of an instrument file.
*
* @author Torsten Mayer-Guerr
* @date 2020-05-29
*
*/
/***********************************************/

#define DOCSTRING_InstrumentType

#include "base/import.h"
#include "config/configRegister.h"
#include "files/fileInstrument.h"
#include "instrumentType.h"

/***** CLASS ***********************************/

// Wrapper class
class InstrumentType
{
public:
  Epoch::Type type;
  InstrumentType(Config &config, const std::string &name);
  static InstrumentType create(Config &config, const std::string &name) {return InstrumentType(config, name);}
};

GROOPS_REGISTER_CLASS_WITHOUT_SUBS(InstrumentType, "instrumentTypeType")

/***********************************************/

template<> Bool readConfig(Config &config, const std::string &name, Epoch::Type &var, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation)
{
  try
  {
    if(isCreateSchema(config))
    {
      config.xselement(name, "instrumentTypeType", mustSet, Config::ONCE, defaultValue, annotation);
      return FALSE;
    }

    if(!hasName(config, name, mustSet))
      return FALSE;
    InstrumentType tmp(config, name);
    var = tmp.type;
    return TRUE;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

InstrumentType::InstrumentType(Config &config, const std::string &name)
{
  try
  {
    type = Epoch::EMPTY;
    std::string choice;

    readConfigChoice(config, name, choice, Config::MUSTSET, "",     "instrument type");
    if(readConfigChoiceElement(config, "INSTRUMENTTIME",    choice, "time without data"))                type = Epoch::INSTRUMENTTIME;
    if(readConfigChoiceElement(config, "MISCVALUE",         choice, "single value"))                     type = Epoch::MISCVALUE;
    if(readConfigChoiceElement(config, "MISCVALUES",        choice, "multiple values"))                  type = Epoch::EMPTY;
    if(readConfigChoiceElement(config, "VECTOR3D",          choice, "x, y, z"))                          type = Epoch::VECTOR3D;
    if(readConfigChoiceElement(config, "COVARIANCE3D",      choice, "xx, yy, zz, xy, xz, yz"))           type = Epoch::COVARIANCE3D;
    if(readConfigChoiceElement(config, "ORBIT",             choice, "position [m], velocity [m/s], acceleration [m/s^2] (each x, y, z)")) type = Epoch::ORBIT;
    if(readConfigChoiceElement(config, "STARCAMERA",        choice, "quaternions (q0, qx, qy, qz)"))     type = Epoch::STARCAMERA;
    if(readConfigChoiceElement(config, "ACCELEROMETER",     choice, "x, y, z [m/s^2]"))                  type = Epoch::ACCELEROMETER;
    if(readConfigChoiceElement(config, "SATELLITETRACKING", choice, "range [m], range rate [m/s], range acceleration [m/s^2]")) type = Epoch::SATELLITETRACKING;
    if(readConfigChoiceElement(config, "GRADIOMETER",       choice, "xx, yy, zz, xy, xz, yz [1/s^2]"))   type = Epoch::GRADIOMETER;
    if(readConfigChoiceElement(config, "GNSSRECEIVER",      choice, "GNSS phase/code observations [m]")) type = Epoch::GNSSRECEIVER;
    if(readConfigChoiceElement(config, "OBSERVATIONSIGMA",  choice, "accuracy"))                         type = Epoch::OBSERVATIONSIGMA;
    if(readConfigChoiceElement(config, "MASS",              choice, ""))                                 type = Epoch::MASS;
    if(readConfigChoiceElement(config, "THRUSTER",          choice, ""))                                 type = Epoch::THRUSTER;
    if(readConfigChoiceElement(config, "MAGNETOMETER",      choice, ""))                                 type = Epoch::MAGNETOMETER;
    if(readConfigChoiceElement(config, "ACCHOUSEKEEPING",   choice, ""))                                 type = Epoch::ACCHOUSEKEEPING;
    endChoice(config);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
