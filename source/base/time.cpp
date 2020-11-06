/***********************************************/
/**
* @file time.cpp
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

#include "base/importStd.h"
#include "base/time.h"

/***********************************************/

Time date2time(UInt year, UInt month, UInt day, UInt hour, UInt minute, Double second)
{
  if(second<0)
    throw(Exception("In date2time:\nNegative seconds"));

  if((month<1)||(month>13))
    throw(Exception("In date2time:\nmonth must between 1<=month<=12"));

  Int Y = static_cast<Int>(year);
  Int M = static_cast<Int>(month);
  Int D = static_cast<Int>(day);

  Int _mjd = (1461 * (Y + 4800 + (M - 14)/12))/4 +(367 * (M - 2 - 12 * ((M - 14)/12)))/12 - (3 * ((Y + 4900 + (M - 14)/12)/100))/4 + D - 32075-2400001;
  return Time(_mjd, ((second/60.+minute)/60.+hour)/24.);
}

/***********************************************/

// I copied this code from internet, but I cannot remember where it was.
void Time::date(UInt &year, UInt &month, UInt &day, UInt &hour, UInt &minute, Double &second) const
{
  // jd0 is the Julian number for noon on the day in question
  // for example   mjd      jd     jd0   === mjd0
  //               3.0  ...3.5  ...4.0   === 3.5
  //               3.3  ...3.8  ...4.0   === 3.5
  //               3.7  ...4.2  ...4.0   === 3.5
  //               3.9  ...4.4  ...4.0   === 3.5
  //               4.0  ...4.5  ...5.0   === 4.5
  long int jd0 = _mjdInt + 2400001;

  // next we convert to Julian dates to make the rest of the maths easier.
  // JD1867217 = 1 Mar 400, so $b is the number of complete Gregorian
  // centuries since then.  The constant 36524.25 is the number of days
  // in a Gregorian century.  The 0.25 on the other constant ensures that
  // $b correctly rounds down on the last day of the 400 year cycle.
  // For example $b == 15.9999... on 2000 Feb 29 not 16.00000.
  Int b = static_cast<Int>(floor((jd0-1867216.25)/36524.25));

  // b-int(b/4) is the number of Julian leap days that are not counted in
  // the Gregorian calendar, and 1402 is the number of days from 1 Jan 4713BC
  // back to 1 Mar 4716BC.  $c represents the date in the Julian calendar
  // corrected back to the start of a leap year cycle.
  Int c = jd0+(b-b/4)+1402;

  // d is the whole number of Julian years from 1 Mar 4716BC to the date
  // we are trying to find.
  Int d = static_cast<Int>(floor((c+0.9)/365.25));

  // e is the number of days from 1 Mar 4716BC to 1 Mar this year
  // using the Julian calendar
  Int e = 365*d+d/4;

  // c-e is now the remaining days in this year from 1 Mar to our date
  // and we need to work out the magic number f such that f-1 == month
  Int f = static_cast<Int>(floor((c-e+123)/30.6001));

  // int(f*30.6001) is the day of the start of the month
  // so the day of the month is the difference between that and c-e+123
  day = c-e+123-static_cast<Int>(floor(30.6001*f));

  // month is now f-1, except that Jan and Feb are f-13
  // ie f 4 5 6 7 8 9 10 11 12 13 14 15
  //    m 3 4 5 6 7 8  9 10 11 12  1  2
  month = (f-2)%12+1;

  // year is d - 4716 (adjusted for Jan and Feb again)
  year = d - 4716 + (month<3);

  hour    = static_cast<UInt>(floor(_mjdMod*24.));
  minute  = static_cast<UInt>(floor((_mjdMod*24-hour)*60));
  second  = ((_mjdMod*24-hour)*60- minute)*60;
}

/***********************************************/

Double Time::decimalYear() const
{
  UInt   year, month, day;
  UInt   hour, minute;
  Double second;
  date(year, month, day, hour, minute, second);

  Time time1 = date2time(year,  1,1);
  Time time2 = date2time(year+1,1,1);

  return year + (*this-time1).mjd()/(time2-time1).mjd();
}

/***********************************************/

UInt Time::dayOfYear() const
{
  UInt   year, month, day;
  UInt   hour, minute;
  Double second;
  date(year, month, day, hour, minute, second);

  Time time1 = date2time(year, 1,1);
  return static_cast<UInt>(std::floor((*this-time1).mjd()+1));
}

/***********************************************/

Bool Time::operator< (const Time &time) const
{
  return (_mjdInt-time._mjdInt) + (_mjdMod-time._mjdMod) < -TIME_EPSILON/86400;
}

/***********************************************/

Time &Time::operator+= (const Time &time)
{
  _mjdInt += time._mjdInt;
  _mjdMod += time._mjdMod;
  Int mjd = static_cast<Int>(std::floor(_mjdMod));
  _mjdInt += mjd;
  _mjdMod -= mjd;
  return *this;
}

/***********************************************/

Time &Time::operator-= (const Time &time)
{
  _mjdInt -= time._mjdInt;
  _mjdMod -= time._mjdMod;
  Int mjd = static_cast<Int>(std::floor(_mjdMod));
  _mjdInt += mjd;
  _mjdMod -= mjd;
  return *this;
}

/***********************************************/

Time &Time::operator*= (Double c)
{
  Double _tmpInt = std::floor(c*_mjdInt);
  Double _tmpMod = c*_mjdInt - _tmpInt;
  _tmpMod += c*_mjdMod;

  _mjdMod = _tmpMod - std::floor(_tmpMod);
  _mjdInt = static_cast<Int>(_tmpInt) + static_cast<Int>(std::floor(_tmpMod));
  return *this;
}

/***********************************************/
/***********************************************/

Time medianSampling(const std::vector<Time> &times)
{
  if(times.size()<2)
    return Time();
  if(times.size() == 2)
    return times.back()-times.front();

  std::vector<Time> diff(times.size()); // first element is times[0], so start from 1
  std::adjacent_difference(times.begin(), times.end(), diff.begin());

  std::partial_sort(diff.begin()+1, diff.begin()+1+diff.size()/2+1, diff.end()); // first element is times[0], so start from 1

  return (diff.size()-1) % 2 == 0 ? diff[diff.size()/2]*0.5 + diff[diff.size()/2+1]*0.5 : diff[diff.size()/2+1];
}

/***********************************************/

Bool isRegular(const std::vector<Time> &times, Double margin)
{
  if(times.size()<=2)
    return TRUE;
  std::vector<Time> diff(times.size());
  std::adjacent_difference(times.begin(), times.end(), diff.begin());
  auto minmax = std::minmax_element(diff.begin()+1, diff.end());
  return ((*minmax.second - *minmax.first).seconds() < margin);
}

/***********************************************/
