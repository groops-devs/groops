/***********************************************/
/**
* @file plotMapProjection.cpp
*
* @brief Map projections used for plotting.
*
* @author Andreas Kvas
* @author Torsten Mayer-Guerr
* @date 2015-10-23
*
*/
/***********************************************/

#define DOCSTRING_PlotMapProjection

#include "base/import.h"
#include "config/configRegister.h"
#include "plotMapProjection.h"

/***********************************************/

static const char *docstringPlotMapProjectionRobinson = R"(
\subsection{Robinson}
The Robinson projection, presented by Arthur H. Robinson in 1963,
is a modified cylindrical projection that is neither conformal nor equal-area.
Central meridian and all parallels are straight lines; other meridians are curved.
It uses lookup tables rather than analytic expressions to make the world map look right.
)";

class PlotMapProjectionRobinson : public PlotMapProjection
{
  Angle L0;

public:
  PlotMapProjectionRobinson(Config &config)
  {
    readConfig(config, "centralMeridian", L0, Config::DEFAULT, "0", "central meridian [degree]");
  }

  std::string scriptEntry(Double width, Double /*height*/) const override
  {
    std::stringstream ss;
    ss<<"-JN"<<L0*RAD2DEG<<"/"<<width<<"c";
    return ss.str();
  }

  Double aspectRatio() const override {return (1+1./80.)/2;}
};

/***********************************************/

static const char *docstringPlotMapProjectionMollweide = R"(
\subsection{Mollweide}
This pseudo-cylindrical, equal-area projection was developed by Mollweide in 1805. Parallels are unequally spaced straight
lines with the meridians being equally spaced elliptical arcs. The scale is only true along latitudes 40$^{o}$44' north and south.
The projection is used mainly for global maps showing data distributions.
)";

class PlotMapProjectionMollweide : public PlotMapProjection
{
  Angle L0;

public:
  PlotMapProjectionMollweide(Config &config)
  {
    readConfig(config, "centralMeridian", L0, Config::DEFAULT, "0", "central meridian [degree]");
  }

  std::string scriptEntry(Double width, Double /*height*/) const override
  {
    std::stringstream ss;
    ss<<"-JW"<<L0*RAD2DEG<<"/"<<width<<"c";
    return ss.str();
  }

  Double aspectRatio() const override {return 0.5;}
};

/***********************************************/

static const char *docstringPlotMapProjectionOrthographic = R"(
\subsection{Orthographic}
The orthographic azimuthal projection is a perspective projection from infinite distance.
It is therefore often used to give the appearance of a globe viewed from space.
)";

class PlotMapProjectionOrthographic : public PlotMapProjection
{
  Angle L0, B0;

public:
  PlotMapProjectionOrthographic(Config &config)
  {
    readConfig(config, "lambdaCenter", L0, Config::DEFAULT, "30", "central point [degree]");
    readConfig(config, "phiCenter",    B0, Config::DEFAULT, "30", "central point [degree]");
  }

  std::string scriptEntry(Double width, Double /*height*/) const override
  {
    std::stringstream ss;
    ss<<"-JG"<<L0*RAD2DEG<<"/"<<B0*RAD2DEG<<"/"<<width<<"c";
    return ss.str();
  }
};

/***********************************************/

static const char *docstringPlotMapProjectionPerspectiveSphere = R"(
\subsection{Perspective sphere}
The orthographic azimuthal projection is a perspective projection from infinite distance.
It is therefore often used to give the appearance of a globe viewed from space.
)";

class PlotMapProjectionPerspectiveSphere : public PlotMapProjection
{
  Angle  L0, B0;
  Angle  azimuth, tilt;
  Angle  viewpointTwist, viewpointWidth, viewpointHeight;
  Double altitude;

public:
  PlotMapProjectionPerspectiveSphere(Config &config)
  {
    readConfig(config, "lambdaCenter",    L0,              Config::DEFAULT,  "10",   "longitude of central point in degrees");
    readConfig(config, "phiCenter",       B0,              Config::DEFAULT,  "40",   "latitude of central point in degrees");
    readConfig(config, "altitude",        altitude,        Config::DEFAULT,  "1500", "[km]");
    readConfig(config, "azimuth",         azimuth,         Config::DEFAULT,  "0",    "to the east of north of view [degrees]");
    readConfig(config, "tilt",            tilt,            Config::DEFAULT,  "0",    "upward tilt of the plane of projection, if negative, then the view is centered on the horizon [degrees]");
    readConfig(config, "viewpointTwist",  viewpointTwist,  Config::DEFAULT,  "0",    "clockwise twist of the viewpoint [degrees]");
    readConfig(config, "viewpointWidth",  viewpointWidth,  Config::DEFAULT,  "0",    "width of the viewpoint [degrees]");
    readConfig(config, "viewpointHeight", viewpointHeight, Config::DEFAULT,  "0",    "height of the viewpoint [degrees]");
  }

  std::string scriptEntry(Double width, Double /*height*/) const override
  {
    std::stringstream ss;
    ss<<"-JG"<<L0*RAD2DEG<<"/"<<B0*RAD2DEG<<"/"<<altitude<<"/"<<azimuth*RAD2DEG<<"/"<<tilt*RAD2DEG
             <<"/"<<viewpointTwist*RAD2DEG<<"/"<<viewpointWidth*RAD2DEG<<"/"<<viewpointHeight*RAD2DEG
             <<"/"<<width<<"c";
    return ss.str();
  }
};

/***********************************************/

static const char *docstringPlotMapProjectionPolar = R"(
\subsection{Polar}
Stereographic projection around given central point.
)";

class PlotMapProjectionPolar : public PlotMapProjection
{
  Angle L0, B0;

public:
  PlotMapProjectionPolar(Config &config)
  {
   readConfig(config, "lambdaCenter", L0, Config::DEFAULT, "0",   "longitude of central point in degrees");
   readConfig(config, "phiCenter",    B0, Config::DEFAULT, "-90", "latitude of central point in degrees");
  }

  std::string scriptEntry(Double width, Double /*height*/) const override
  {
    std::stringstream ss;
    ss<<"-JS"<<L0*RAD2DEG<<"/"<<B0*RAD2DEG<<"/"<<width<<"c";
    return ss.str();
  }
};

/***********************************************/

static const char *docstringPlotMapProjectionSkyplot = R"(
\subsection{Skyplot}
Skyplot used to plot azimuth/elevation data as generated by
\program{GnssAntennaDefinition2Skyplot} or \program{GnssResiduals2Skyplot}.
)";

class PlotMapProjectionSkyplot : public PlotMapProjection
{
public:
  PlotMapProjectionSkyplot(Config &/*config*/) {}

  std::string scriptEntry(Double width, Double /*height*/) const override
  {
    std::stringstream ss;
    ss<<"-JPa"<<width<<"cr"; // --MAP_POLAR_CAP=none";
    return ss.str();
  }
};

/***********************************************/

static const char *docstringPlotMapProjectionUtm = R"(
\subsection{UTM}
A particular subset of the transverse Mercator is the Universal Transverse Mercator (UTM)
which was adopted by the US Army for large-scale military maps.
Here, the globe is divided into 60 zones between 84$^{o}$S and 84$^{o}$N, most of which are 6$^{o}$ wide.
Each of these UTM zones have their unique central meridian.
)";

class PlotMapProjectionUtm : public PlotMapProjection
{
  std::string zone;

public:
  PlotMapProjectionUtm(Config &config)
  {
    readConfig(config, "zone", zone, Config::MUSTSET, "", "UTM zone code (e.g. 33N)");
  }

  std::string scriptEntry(Double width, Double /*height*/) const override
  {
    std::stringstream ss;
    ss<<"-JU"<<zone<<"/"<<width<<"c";
    return ss.str();
  }

  virtual Double aspectRatio() const override {return 0;}
};

/***********************************************/

static const char *docstringPlotMapProjectionLambert = R"(
\subsection{Lambert}
This conic projection was designed by Lambert (1772) and has been used extensively for mapping of regions with predominantly east-west orientation.
)";

class PlotMapProjectionLambert : public PlotMapProjection
{
  Angle L0, B0, par1, par2;

public:
  PlotMapProjectionLambert(Config &config)
  {
    readConfig(config, "lambda0", L0,   Config::MUSTSET, "", "longitude of projection center [deg]");
    readConfig(config, "phi0",    B0,   Config::MUSTSET, "", "latitude of projection centert [deg]");
    readConfig(config, "phi1",    par1, Config::MUSTSET, "", "latitude of first standard parallel [deg]");
    readConfig(config, "phi2",    par2, Config::MUSTSET, "", "latitude of first standard parallel [deg]");
  }

  std::string scriptEntry(Double width, Double /*height*/) const override
  {
    std::stringstream ss;
    ss<<"-JL"<<L0*RAD2DEG<<"/"<<B0*RAD2DEG<<"/"<<par1*RAD2DEG<<"/"<<par2*RAD2DEG<<"/"<<width<<"c";
    return ss.str();
  }
};

/***********************************************/

static const char *docstringPlotMapProjectionLinear = R"(
\subsection{Linear}
Linear mapping of longitude/latitude to x/y (Plate Caree).
)";

class PlotMapProjectionLinear : public PlotMapProjection
{
public:
  PlotMapProjectionLinear(Config &/*config*/) {}

  std::string scriptEntry(Double width, Double height) const override
  {
    std::stringstream ss;
    ss<<"-JX"<<width<<"cd/"<<height<<"cd";
    return ss.str();
  }

  virtual Double aspectRatio() const override {return 0;}
};

/***********************************************/
/***********************************************/

GROOPS_REGISTER_CLASS(PlotMapProjection, "plotMapProjectionType",
                      PlotMapProjectionRobinson,
                      PlotMapProjectionOrthographic,
                      PlotMapProjectionPerspectiveSphere,
                      PlotMapProjectionPolar,
                      PlotMapProjectionSkyplot,
                      PlotMapProjectionUtm,
                      PlotMapProjectionLambert,
                      PlotMapProjectionLinear,
                      PlotMapProjectionMollweide)

GROOPS_READCONFIG_CLASS(PlotMapProjection, "plotMapProjectionType")

/***********************************************/

PlotMapProjectionPtr PlotMapProjection::create(Config &config, const std::string &name)
{
  try
  {
    PlotMapProjectionPtr proj;
    std::string  type;

    readConfigChoice(config, name, type, Config::MUSTSET, "", "plot layers");
    if(readConfigChoiceElement(config, "robinson",     type, "robinson projection"))
      proj = PlotMapProjectionPtr(new PlotMapProjectionRobinson(config));
    if(readConfigChoiceElement(config, "orthographic", type, "orthographic projection"))
      proj = PlotMapProjectionPtr(new PlotMapProjectionOrthographic(config));
    if(readConfigChoiceElement(config, "perspective",  type, "perspective sphere"))
      proj = PlotMapProjectionPtr(new PlotMapProjectionPerspectiveSphere(config));
    if(readConfigChoiceElement(config, "polar",        type, "polar stereographic"))
      proj = PlotMapProjectionPtr(new PlotMapProjectionPolar(config));
    if(readConfigChoiceElement(config, "skyplot",      type, "skyplot"))
      proj = PlotMapProjectionPtr(new PlotMapProjectionSkyplot(config));
    if(readConfigChoiceElement(config, "linear",       type, "linear (plate carree) projection"))
      proj = PlotMapProjectionPtr(new PlotMapProjectionLinear(config));
    if(readConfigChoiceElement(config, "Utm",          type, "UTM projection"))
      proj = PlotMapProjectionPtr(new PlotMapProjectionUtm(config));
    if(readConfigChoiceElement(config, "lambert",      type, "lambert conformal conic"))
      proj = PlotMapProjectionPtr(new PlotMapProjectionLambert(config));
    if(readConfigChoiceElement(config, "mollweide",    type, "mollweide projection "))
      proj = PlotMapProjectionPtr(new PlotMapProjectionMollweide(config));
    endChoice(config);

    return proj;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
