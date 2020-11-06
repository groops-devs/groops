/***********************************************/
/**
* @file border.cpp
*
* @brief Borders of an area on sphere/ellipsoid.
*
* @author Torsten Mayer-Guerr
* @date 2004-10-28
*
*/
/***********************************************/

#define DOCSTRING_Border

#include "base/import.h"
#include "config/configRegister.h"
#include "borderRectangle.h"
#include "borderCap.h"
#include "borderPolygon.h"
#include "borderGlobal.h"
#include "border.h"

/***********************************************/

GROOPS_REGISTER_CLASS(Border, "borderType",
                      BorderRectangle,
                      BorderCap,
                      BorderPolygon)

GROOPS_READCONFIG_UNBOUNDED_CLASS(Border, "borderType")

/***********************************************/

Border::Border(Config &config, const std::string &name)
{
  try
  {
    std::string type;
    while(readConfigChoice(config, name, type, Config::OPTIONAL, "", "borders of geographical areas"))
    {
      if(readConfigChoiceElement(config, "rectangle", type, "along lines of geographical coordinates"))
        border.push_back(new BorderRectangle(config));
      if(readConfigChoiceElement(config, "cap", type, "spherical cap"))
        border.push_back(new BorderCap(config));
      if(readConfigChoiceElement(config, "polygon", type, "polygon from file"))
        border.push_back(new BorderPolygon(config));
      endChoice(config);
      if(isCreateSchema(config))
        return;
    };

    if(border.empty())
      border.push_back(new BorderGlobal());
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Border::~Border()
{
  for(UInt i=0; i<border.size(); i++)
    delete border.at(i);
}

/***********************************************/

Bool Border::isInnerPoint(Angle lambda, Angle phi) const
{
  Bool inner = border.at(0)->isExclude();
  for(UInt i=0; i<border.size(); i++)
    if(border.at(i)->isInnerPoint(lambda,phi))
      inner = !border.at(i)->isExclude();
  return inner;
}

/***********************************************/

Bool Border::isInnerPoint(const Vector3d &point, const Ellipsoid &ellipsoid) const
{
  Angle  L,B;
  Double h;
  ellipsoid(point, L, B, h);
  return isInnerPoint(L, B);
}

/***********************************************/
