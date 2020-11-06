/***********************************************/
/**
* @file format.h
*
* @brief Format of numbers.
*
* @author Torsten Mayer-Guerr
* @author Sebastian Strasser
* @date 2016-07-15
*
*/
/***********************************************/

#ifndef __GROOPS_FORMAT__
#define __GROOPS_FORMAT__

#include "base/importStd.h"

using namespace std::literals::string_literals;

/***** FUNCTIONS *******************************/

/** @brief format numbers for output.
* @ingroup base
* The number @a x is replaced by the @a format string. The format contains the text to be written
* to output It can contain embedded format specifiers that are replaced by the number
* and formatted as requested (also multiple times).
* In the following an example output is given in the brackets for the example value of 57493.8:
* - %i: integer [57494]
* - %f: Decimal floating point [57493.800000]
* - %e: Scientific notation [5.749380e+04]
* - %g: Use the shortest representation: %e or %f [57493.8]
* - %c: interpret number as ASCII character.
* - %%: write a single %.
*
* The following specifiers interpret the value as MJD (modified julian date):
* - %y: four digit year [2016]
* - %Y: two digit year [16]
* - %m: month [04]
* - %d: day of month [15]
* - %H: Hour [19]
* - %M: Minute [12]
* - %S: Second [00 or custom floating point format, e.g. %04.1S = 00.0]
* - %D: Date (same as %y-%m-%d) [2016-04-15]
* - %T: Time (same as %H-%M-%S) [19-12-00]
* - %W: GPS week [1892]
* - %w: Day of GPS week (0..6) [5]
* - %O: Day of year (1..366)
* ẂARNING: Using multiple %S and/or %T specifiers with different precisions in the same format string
*          can lead to incorrect rounding (60 sec) since the maximum precision is used for rounding!
*
* The output can be specifed further with %[width][.precision]specifier,
* there [width] is the minimum number of characters to be printed. If the value to be printed
* is shorter than this number, the result is padded with blank spaces (or zeros if [width] starts
* with a zero). The [.precision] defines the number of digits after the period (for %g the number
* of significant digits instead).
*
* Example 51544.5 % "%D_T%T"s -> '2000-01-01_T12-00-00'. */
std::string operator%(LongDouble x, const std::string &format);

inline std::string operator%(Double x, const std::string &format) {return static_cast<LongDouble>(x)%format;}
inline std::string operator%(UInt   x, const std::string &format) {return static_cast<LongDouble>(x)%format;}
inline std::string operator%(Int    x, const std::string &format) {return static_cast<LongDouble>(x)%format;}

class Time;
/** @copydoc operator%(LongDouble, const std::string &) */
std::string operator%(const Time &x, const std::string &format);

/*************************************************/

#endif
