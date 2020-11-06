/***********************************************/
/**
* @file borderGlobal.h
*
* @brief Complete sphere.
* @see Border
*
* @author Torsten Mayer-Guerr
* @date 2004-10-28
*
*/
/***********************************************/

#ifndef __GROOPS_BORDERGLOBAL__
#define __GROOPS_BORDERGLOBAL__

#include "classes/border/border.h"

/***** CLASS ***********************************/

/** @brief Complete sphere: All points are inside.
* @ingroup borderGroup
* @see Border */
class BorderGlobal : public BorderBase
{
public:
  Bool isInnerPoint(Angle, Angle) const {return TRUE;}
  Bool isExclude() const {return FALSE;}
};

/***********************************************/

#endif /* __GROOPS_BORDER__ */
