/***********************************************/
/**
* @file constants.cpp
*
* @brief Define constants.
*
* Constants are only valid after call of main().
* New constants must be set in inputOutput/settings.cpp.
*
* @author Torsten Mayer-Guerr
* @date 2010-03-01
*
*/
/***********************************************/

#include "base/importStd.h"
#include "base/time.h"
#include "base/constants.h"

/***** CONSTANTS *******************************/

Double LIGHT_VELOCITY;

/***** DEFAULTS ********************************/

Double DEFAULT_GRS80_a;
Double DEFAULT_GRS80_f;
Double DEFAULT_GM;
Double DEFAULT_R;
std::string STRING_DEFAULT_GRS80_a;
std::string STRING_DEFAULT_GRS80_f;
std::string STRING_DEFAULT_GM;
std::string STRING_DEFAULT_R;

/***** CONSTANTS for Planets and Tides *********/

Double GRAVITATIONALCONSTANT;
Double R_Earth;
Double R_Moon;
Double GM_Earth;
Double GM_Sun;
Double GM_Moon;
Double GM_MERCURY;
Double GM_VENUS;
Double GM_MARS;
Double GM_JUPITER;
Double GM_SATURN;

/***** TIME SHIFTS *****************************/

Double TIME_EPSILON;
Double DELTA_TAI_GPS;
Double DELTA_TT_GPS;
Double J2000;
std::string STRING_J2000;

std::vector<Int>    MJD_UTC_GPS;
std::vector<Double> DELTA_UTC_GPS;

/***** CLASS ***********************************/

// class for constants initialization,
// will be initialized before main()
class InitConstants
{
public:
  InitConstants();
};

static InitConstants initConstants;

/***********************************************/

InitConstants::InitConstants()
{
  try
  {
    LIGHT_VELOCITY = 299792458.0;

    STRING_DEFAULT_GRS80_a = "6378137.0";
    STRING_DEFAULT_GRS80_f = "298.2572221010";
    STRING_DEFAULT_GM      = "3.986004415e+14";
    STRING_DEFAULT_R       = "6378136.3";

    DEFAULT_GM      = std::atof(STRING_DEFAULT_GM.c_str());
    DEFAULT_R       = std::atof(STRING_DEFAULT_R.c_str());
    DEFAULT_GRS80_a = std::atof(STRING_DEFAULT_GRS80_a.c_str());
    DEFAULT_GRS80_f = std::atof(STRING_DEFAULT_GRS80_f.c_str());

    GRAVITATIONALCONSTANT = 6.673e-11;
    R_Earth    = DEFAULT_R;
    R_Moon     = 1738000;
    GM_Earth   = DEFAULT_GM;
    GM_Sun     = 1.32712442076e20;
    GM_Moon    = 0.49028010560e13; //GM_Earth/81.30056;
    GM_MERCURY = GM_Sun/6023600.0;
    GM_VENUS   = GM_Sun/408523.71;
    GM_MARS    = GM_Sun/3098708.0;
    GM_JUPITER = GM_Sun/1047.3486;
    GM_SATURN  = GM_Sun/3497.898;

    TIME_EPSILON  = 1e-5;
    DELTA_TAI_GPS = 19.0;
    DELTA_TT_GPS  = 51.184;
    STRING_J2000  = "51544.5";
    J2000         = std::atof(STRING_J2000.c_str());

    MJD_UTC_GPS.push_back(date2time(2017, 1, 1).mjdInt());  DELTA_UTC_GPS.push_back(-18);
    MJD_UTC_GPS.push_back(date2time(2015, 7, 1).mjdInt());  DELTA_UTC_GPS.push_back(-17);
    MJD_UTC_GPS.push_back(date2time(2012, 7, 1).mjdInt());  DELTA_UTC_GPS.push_back(-16);
    MJD_UTC_GPS.push_back(date2time(2009, 1, 1).mjdInt());  DELTA_UTC_GPS.push_back(-15);
    MJD_UTC_GPS.push_back(date2time(2006, 1, 1).mjdInt());  DELTA_UTC_GPS.push_back(-14);
    MJD_UTC_GPS.push_back(date2time(1999, 1, 1).mjdInt());  DELTA_UTC_GPS.push_back(-13);
    MJD_UTC_GPS.push_back(date2time(1997, 7, 1).mjdInt());  DELTA_UTC_GPS.push_back(-12);
    MJD_UTC_GPS.push_back(date2time(1996, 1, 1).mjdInt());  DELTA_UTC_GPS.push_back(-11);
    MJD_UTC_GPS.push_back(date2time(1994, 7, 1).mjdInt());  DELTA_UTC_GPS.push_back(-10);
    MJD_UTC_GPS.push_back(date2time(1993, 7, 1).mjdInt());  DELTA_UTC_GPS.push_back(-9);
    MJD_UTC_GPS.push_back(date2time(1992, 7, 1).mjdInt());  DELTA_UTC_GPS.push_back(-8);
    MJD_UTC_GPS.push_back(date2time(1991, 1, 1).mjdInt());  DELTA_UTC_GPS.push_back(-7);
    MJD_UTC_GPS.push_back(date2time(1990, 1, 1).mjdInt());  DELTA_UTC_GPS.push_back(-6);
    MJD_UTC_GPS.push_back(date2time(1988, 1, 1).mjdInt());  DELTA_UTC_GPS.push_back(-5);
    MJD_UTC_GPS.push_back(date2time(1985, 7, 1).mjdInt());  DELTA_UTC_GPS.push_back(-4);
    MJD_UTC_GPS.push_back(date2time(1983, 7, 1).mjdInt());  DELTA_UTC_GPS.push_back(-3);
    MJD_UTC_GPS.push_back(date2time(1982, 7, 1).mjdInt());  DELTA_UTC_GPS.push_back(-2);
    MJD_UTC_GPS.push_back(date2time(1981, 7, 1).mjdInt());  DELTA_UTC_GPS.push_back(-1);
    MJD_UTC_GPS.push_back(date2time(1980, 1, 1).mjdInt());  DELTA_UTC_GPS.push_back(0);
    MJD_UTC_GPS.push_back(date2time(1979, 1, 1).mjdInt());  DELTA_UTC_GPS.push_back(1);
    MJD_UTC_GPS.push_back(date2time(1978, 1, 1).mjdInt());  DELTA_UTC_GPS.push_back(2);
    MJD_UTC_GPS.push_back(date2time(1977, 1, 1).mjdInt());  DELTA_UTC_GPS.push_back(3);
    MJD_UTC_GPS.push_back(date2time(1976, 1, 1).mjdInt());  DELTA_UTC_GPS.push_back(4);
    MJD_UTC_GPS.push_back(date2time(1975, 1, 1).mjdInt());  DELTA_UTC_GPS.push_back(5);
    MJD_UTC_GPS.push_back(date2time(1974, 1, 1).mjdInt());  DELTA_UTC_GPS.push_back(6);
    MJD_UTC_GPS.push_back(date2time(1973, 1, 1).mjdInt());  DELTA_UTC_GPS.push_back(7);
    MJD_UTC_GPS.push_back(date2time(1972, 7, 1).mjdInt());  DELTA_UTC_GPS.push_back(8);
    MJD_UTC_GPS.push_back(date2time(1972, 1, 1).mjdInt());  DELTA_UTC_GPS.push_back(9);
    MJD_UTC_GPS.push_back(0);                               DELTA_UTC_GPS.push_back(10);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
