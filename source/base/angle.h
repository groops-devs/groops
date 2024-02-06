/***********************************************/
/**
* @file angle.h
*
* @brief Angle in radians.
*
* @author Torsten Mayer-Guerr
* @date 2001-31-05
*
*/
/***********************************************/

#ifndef __GROOPS_ANGLE__
#define __GROOPS_ANGLE__

#include "base/importStd.h"

/***** CLASS ***********************************/

/** @brief Angle in radians.
* @ingroup vector3dGroup
* Its only used to distinguish between doubles and angles in input/output.
* Angle is converted from/to degrees in input/output. */
class Angle
{
  Double value;

public:
  explicit Angle(Double x=0.) : value(x) {}   //!< Default constructor

  Angle &operator=(Double x) {value=x; return *this;} //!< Assignment.

  operator Double() const {return value;} //!< Cast to Double.

  Angle &operator*= (Double c)       {value *= c;       return *this;}
  Angle &operator+= (const Angle &x) {value += x.value; return *this;}
  Angle &operator-= (const Angle &x) {value -= x.value; return *this;}
};

/***********************************************/

inline Angle operator- (const Angle &t)                   {return Angle(t)  *= -1;}
inline Angle operator+ (const Angle &t1, const Angle &t2) {return Angle(t1) += t2;}
inline Angle operator- (const Angle &t1, const Angle &t2) {return Angle(t1) -= t2;}
inline Angle operator* (Double c, const Angle &t)         {return Angle(t)  *=c;}
inline Angle operator* (const Angle &t, Double c)         {return Angle(t)  *=c;}

/*************************************************/

#endif
