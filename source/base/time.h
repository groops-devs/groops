/***********************************************/
/**
* @file time.h
*
* @brief time representation.
*
* @author Torsten Mayer-Guerr
* @author Andreas Kvas
* @author Sebastian Strasser
* @date 2002-11-12
*
*/
/***********************************************/

#ifndef __GROOPS_TIME__
#define __GROOPS_TIME__

#include "base/importStd.h"
#include "base/constants.h"

/**
* @defgroup timeGroup Time
* @brief Classes and functions for accurate time representation.
* @ingroup base */
/// @{

/***** CLASS ***********************************/

/** @brief time representation.
* Precise representation: integer days and additional double for fraction of a day,
* Implemented Comparable and LinearSpace. */
class Time
{
  Int    _mjdInt;  // integer days
  Double _mjdMod;  // fraction of a day

public:
  /// Default Constructor.
  Time() : _mjdInt(0), _mjdMod(0.0) {}

  /// Constructor.
  explicit Time(Int mjdInt, Double mjdMod);

  /// Time as MJD (modified julian date)
  Double mjd() const;

  /// Time in seconds.
  Double seconds() const;

  /// Decimal year e.g. 2008.135
  Double decimalYear() const;

  /// Day of year (DOY)
  UInt dayOfYear() const;

  /** @brief The Date of Time.
  * @param[out] year year
  * @param[out] month month
  * @param[out] day day
  * @param[out] hour hour
  * @param[out] minute minutes
  * @param[out] second seconds with frations */
  void date(UInt &year, UInt &month, UInt &day, UInt &hour, UInt &minute, Double &second) const;

  /// Date as string (yyyy-mm-dd).
  inline std::string dateStr() const {return *this%"%y-%m-%d"s;}

  /** @brief Date and time as string (yyyy-mm-dd_hh-mm-ss).
  * Useful for file names */
  inline std::string dateTimeStr() const {return *this%"%y-%m-%d_%H-%M-%S"s;}

  /// Time of day as string (hh-mm-ss)
  inline std::string str() const {return *this%"%H-%M-%S"s;}

  /** @brief Started day in MJD without fractional.
  * Example: gives 2 for mjd=2.25 and -3 for mjd=-2.25. */
  Int mjdInt() const;

  /** @brief Fractional part of started day.
  * Example: gives 0.25 for mjd=2.25 and 0.75 for mjd=-2.25. */
  Double mjdMod() const;

  /** @brief Returns TRUE if calling instance is within given interval [@p intervalStart, @p intervalEnd). */
  Bool isInInterval(const Time &intervalStart, const Time &intervalEnd) const {return (intervalStart<=*this) && (*this<intervalEnd);}

  Time &operator-= (const Time &time);
  Time &operator+= (const Time &time);
  Time &operator*= (Double c);
  Time  operator+  () const {return *this;}
  Time  operator-  () const {return Time(-_mjdInt, -_mjdMod);}

  /** @brief Comparisons.
  * Times considered to be equal within @a TIME_EPSILON seconds margin. */
  Bool operator<  (const Time &time) const;
  Bool operator== (const Time &time) const {return !(*this<time) && !(time<*this);}
  Bool operator!= (const Time &time) const {return !(*this==time);}
  Bool operator>  (const Time &time) const {return !(*this<time) && !(*this==time);}
  Bool operator<= (const Time &time) const {return !(*this>time);}
  Bool operator>= (const Time &time) const {return !(*this<time);}
};

/***** FUNCTIONS *******************************/

/// MJD (modified julian date) to Time representation.
inline Time mjd2time(LongDouble mjd);

/// Seconds to Time representation.
inline Time seconds2time(LongDouble seconds);

/// Date to Time representation.
Time date2time(UInt year, UInt month, UInt day, UInt hour=0, UInt minute=0, Double second=0.0);

/** @brief Conversion from date to Time representation.
* @see Time::date */
inline void time2date(const Time &time, UInt &year, UInt &month, UInt &day, UInt &hour, UInt &minute, Double &second);

/// @brief Median sampling of a time series.
Time medianSampling(const std::vector<Time> &times);

/** @brief return true if all points are equally spaced.
* @param times series to be tested
* @param margin threshold for equality [seconds] */
Bool isRegular(const std::vector<Time> &times, Double margin=1e-5);

/// GPS time to TT time.
inline Time timeGPS2TT(const Time &timeGPS);

/// GPS time to julian centuries (JC).
inline Double timeGPS2JC(const Time &timeGPS);

/// TT to GPS-Time.
inline Time timeTT2GPS(const Time &timeTT);

/** @brief GPS time to Universal Time Coordinated (UTC).
* Need leap seconds from constants.cpp */
inline Time timeGPS2UTC(const Time &timeGPS);

/** @brief Universal Time Coordinated (UTC) to GPS time.
* Need leap seconds from constants.cpp */
inline Time timeUTC2GPS(const Time &timeUTC);

inline Time operator+ (const Time &t1, const Time &t2) {return Time(t1) += t2;}
inline Time operator- (const Time &t1, const Time &t2) {return Time(t1) -= t2;}
inline Time operator* (Double c, const Time &t)        {return Time(t)  *=c;}
inline Time operator* (const Time &t, Double c)        {return Time(t)  *=c;}

/// @} group time

/***********************************************/
/***** INLINES   *******************************/
/***********************************************/

inline Time::Time(Int mjdInt, Double mjdMod) : _mjdInt(mjdInt + static_cast<Int>(floor(mjdMod))), _mjdMod(mjdMod - static_cast<Int>(floor(mjdMod)))
{
}

inline Double Time::mjd()     const {return _mjdInt + _mjdMod;}
inline Int    Time::mjdInt()  const {return _mjdInt;}
inline Double Time::mjdMod()  const {return _mjdMod;}
inline Double Time::seconds() const {return mjd()*(24*60*60);}

inline Time mjd2time(LongDouble mjd)
{
  Int _mjdInt = static_cast<Int>(std::floor(mjd));
  return Time(_mjdInt, static_cast<Double>(mjd-_mjdInt));
}

inline Time seconds2time(LongDouble seconds)
{
  return mjd2time(seconds/(24*60*60));
}

inline void time2date(const Time &time, UInt &year, UInt &month, UInt &day, UInt &hour, UInt &minute, Double &second)
{
  time.date(year,month,day,hour,minute,second);
}

inline Time timeGPS2TT (const Time &timeGPS)
{
  return timeGPS + seconds2time(DELTA_TT_GPS);
}

inline Double timeGPS2JC (const Time &timeGPS)
{
  return (timeGPS2TT(timeGPS)-mjd2time(J2000)).mjd()/36525.;
}

inline Time timeTT2GPS (const Time &timeTT)
{
  return timeTT - seconds2time(DELTA_TT_GPS);
}

inline Time timeGPS2UTC(const Time &timeGPS)
{
  const Double mjd = timeGPS.mjd();
  for(UInt i=0; i<DELTA_UTC_GPS.size(); i++)
    if(mjd >= MJD_UTC_GPS.at(i)-DELTA_UTC_GPS.at(i)/86400.)
      return timeGPS + seconds2time(DELTA_UTC_GPS.at(i));
  throw(Exception("no leap seconds available at "+timeGPS.dateTimeStr()));
}

inline Time timeUTC2GPS(const Time &timeUTC)
{
  const Double mjd = timeUTC.mjd();
  for(UInt i=0; i<DELTA_UTC_GPS.size(); i++)
    if(mjd >= MJD_UTC_GPS.at(i))
      return timeUTC - seconds2time(DELTA_UTC_GPS.at(i));
  throw(Exception("no leap seconds available at "+timeUTC.dateTimeStr()));
}

/***********************************************/

#endif /* __GROOPS_TIME__ */
