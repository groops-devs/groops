/***********************************************/
/**
* @file borderRectangle.h
*
* @brief Rectangle: along lines of geographical coordinates.
* @see Border
*
* @author Torsten Mayer-Guerr
* @date 2004-10-28
*
*/
/***********************************************/

#ifndef __GROOPS_BORDERRECTANGLE__
#define __GROOPS_BORDERRECTANGLE__

// Latex documentation
#ifdef DOCSTRING_Border
static const char *docstringBorderRectangle = R"(
\subsection{Rectangle}
The region is restricted along lines of geographical coordinates.
\config{minPhi} and \config{maxPhi} describe the lower and the upper bound of the region.
\config{minLambda} and \config{maxLambda} define the left and right bound.
)";
#endif

/***********************************************/

#include "base/vector3d.h"
#include "config/config.h"
#include "classes/border/border.h"

/***** CLASS ***********************************/

/** @brief Rectangle: along lines of geographical coordinates.
* @ingroup borderGroup
* @see Border */
class BorderRectangle : public BorderBase
{
  Angle minLambda, maxLambda;
  Angle minPhi,    maxPhi;
  Bool  exclude;

public:
  BorderRectangle(Config &config);

  Bool isInnerPoint(Angle lambda, Angle phi) const;
  Bool isExclude() const {return exclude;}
};

/***********************************************/

inline BorderRectangle::BorderRectangle(Config &config)
{
  readConfig(config, "minLambda", minLambda, Config::MUSTSET,  "", "");
  readConfig(config, "maxLambda", maxLambda, Config::MUSTSET,  "", "");
  readConfig(config, "minPhi",    minPhi,    Config::MUSTSET,  "", "");
  readConfig(config, "maxPhi",    maxPhi,    Config::MUSTSET,  "", "");
  readConfig(config, "exclude",   exclude,   Config::DEFAULT,  "0", "dismiss points inside");
  if(isCreateSchema(config)) return;

  // valid values: -180<=lambda<=180, -90<=phi<=90 ??
  if(minLambda<-PI)
    minLambda += Angle(2*PI);

  if(maxLambda>PI)
    maxLambda += Angle(-2*PI);
}

/***********************************************/

inline Bool BorderRectangle::isInnerPoint(Angle lambda, Angle phi) const
{
  if((phi<minPhi)||(phi>maxPhi))
    return FALSE;
  if(minLambda<maxLambda)
    return ((lambda>=minLambda) && (lambda<=maxLambda));
  else
    return ((lambda>=minLambda) || (lambda<=maxLambda));
}

/***********************************************/

#endif /* __GROOPS_BORDER__ */
