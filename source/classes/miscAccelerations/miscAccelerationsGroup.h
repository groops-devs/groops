/***********************************************/
/**
* @file miscAccelerationsGroup.h
*
* @brief Group.
* @see MiscAccelerations
*
* @author Torsten Mayer-Guerr
* @date 2023-11-13
*
*/
/***********************************************/

#ifndef __GROOPS_MISCACCELERATIONSGROUP__
#define __GROOPS_MISCACCELERATIONSGROUP__

// Latex documentation
#ifdef DOCSTRING_MiscAccelerations
static const char *docstringMiscAccelerationsGroup = R"(
\subsection{Group}\label{miscAccelerationsType:group}
Groups a set of \configClass{miscAccelerations}{miscAccelerationsType} and has no further effect itself.
)";
#endif

/***********************************************/

#include "classes/miscAccelerations/miscAccelerations.h"

/***** CLASS ***********************************/

/** @brief Group.
* @ingroup miscAccelerationsGroup
* @see MiscAccelerations */
class MiscAccelerationsGroup : public MiscAccelerationsBase
{
  MiscAccelerationsPtr miscAccelerations;
  Double               factor;

public:
  MiscAccelerationsGroup(Config &config);

  Vector3d acceleration(SatelliteModelPtr satellite, const Time &time, const Vector3d &position, const Vector3d &velocity,
                        const Rotary3d &rotSat, const Rotary3d &rotEarth, EphemeridesPtr ephemerides) override;
};

/***********************************************/

inline MiscAccelerationsGroup::MiscAccelerationsGroup(Config &config)
{
  try
  {
    readConfig(config, "miscAccelerations", miscAccelerations, Config::DEFAULT, "", "");
    readConfig(config, "factor",            factor,            Config::DEFAULT, "1.0", "the result is multiplied by this factor");
    if(isCreateSchema(config)) return;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline Vector3d MiscAccelerationsGroup::acceleration(SatelliteModelPtr satellite, const Time &time,
                                                             const Vector3d &pos, const Vector3d &vel,
                                                             const Rotary3d &rotSat, const Rotary3d &rotEarth, EphemeridesPtr ephemerides)
{
  try
  {
    return factor * miscAccelerations->acceleration(satellite, time, pos, vel, rotSat, rotEarth, ephemerides);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

#endif
