/***********************************************/
/**
* @file gnssAntennaDefinition2Skyplot.cpp
*
* @brief Create skyplot from Antenna definition.
*
* @author Torsten Mayer-Guerr
* @date 2012-11-16
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
Produce a \file{skyplot}{griddedData} of antenna center variations
which can be plotted with \program{PlotMap}.

The first antenna from \configFile{inputfileAntennaDefinition}{gnssAntennaDefinition}
matching the wildcard patterns of \config{name}, \config{serial}, \config{radome} is used.

For each antenna pattern (gnssType) a separate data column is computed.
A subset of patterns can be selected with \configClass{types}{gnssType}.

Azimuth and elevation are written as ellipsoidal longitude and latitude in a \file{griddedData file}{griddedData}.
The choosen ellipsoid parameters \config{R} and \config{inverseFlattening} are arbitrary but should be the same
as in \configClass{grid}{gridType} and \program{PlotMap}.

\fig{!hb}{1.0}{fileFormatGnssAntennaDefinition}{fig:gnssAntennaDefinition2Skyplot}{Antenna Center Variations of ASH701945D\_M for two frequencies of GPS and GLONASS}
)";

/***********************************************/

#include "programs/program.h"
#include "base/string.h"
#include "files/fileGriddedData.h"
#include "files/fileGnssStationInfo.h"
#include "classes/grid/grid.h"
#include "misc/miscGriddedData.h"

/***** CLASS ***********************************/

/** @brief Create gridded data from Antenna definition.
* @ingroup programsGroup */
class GnssAntennaDefinition2Skyplot
{
public:
  void run(Config &config);
};

GROOPS_REGISTER_PROGRAM(GnssAntennaDefinition2Skyplot, SINGLEPROCESS, "Create gridded data from Antenna definition.", Gnss, Grid)
GROOPS_RENAMED_PROGRAM(GnssAntennaDefinition2GriddedData, GnssAntennaDefinition2Skyplot, date2time(2019, 9, 9))

/***********************************************/

void GnssAntennaDefinition2Skyplot::run(Config &config)
{
  try
  {
    FileName              fileNameGrid;
    FileName              fileNameAntenna;
    std::string           name, serial, radome;
    GridPtr               grid;
    std::vector<GnssType> types;
    Double                a, f;

    readConfig(config, "outputfileGriddedData",      fileNameGrid,    Config::MUSTSET,   "",  "data column for each gnssType");
    readConfig(config, "inputfileAntennaDefinition", fileNameAntenna, Config::MUSTSET,   "",  "");
    readConfig(config, "grid",                       grid,            Config::MUSTSET,   "",  "");
    readConfig(config, "name",                       name,            Config::OPTIONAL,  "*", "");
    readConfig(config, "serial",                     serial,          Config::OPTIONAL,  "*", "");
    readConfig(config, "radome",                     radome,          Config::OPTIONAL,  "*", "");
    readConfig(config, "types",                      types,           Config::OPTIONAL,  "",  "if not set, all types in the file are used");
    readConfig(config, "R",                          a,               Config::DEFAULT,  STRING_DEFAULT_GRS80_a, "reference radius for ellipsoidal coordinates");
    readConfig(config, "inverseFlattening",          f,               Config::DEFAULT,  STRING_DEFAULT_GRS80_f, "reference flattening for ellipsoidal coordinates");
    if(isCreateSchema(config)) return;

    // ============================

    logStatus<<"read antenna center variations <"<<fileNameAntenna<<">"<<Log::endl;
    std::vector<GnssAntennaDefinitionPtr> antennaList;
    readFileGnssAntennaDefinition(fileNameAntenna, antennaList);

    std::regex pattern = String::wildcard2regex(GnssAntennaDefinition::str(name, serial, radome));
    auto iter = std::find_if(antennaList.begin(), antennaList.end(), [&](const GnssAntennaDefinitionPtr &ant){return std::regex_match(ant->str(), pattern);});
    if(iter == antennaList.end())
      throw(Exception("antenna definition not found"));
    GnssAntennaDefinitionPtr antenna = *iter;

    if(types.size() == 0)
      for(auto &pattern : antenna->pattern)
        types.push_back(pattern.type);

    logStatus<<"create values on grid"<<Log::endl;
    std::vector<Vector3d> points = grid->points();
    std::vector<Double>   areas  = grid->areas();
    std::vector< std::vector<Double> > values(types.size(), std::vector<Double>(points.size()));
    for(UInt i=0; i<points.size(); i++)
      for(UInt idType=0; idType<types.size(); idType++)
      {
        const UInt idPattern = antenna->findAntennaPattern(types.at(idType), GnssAntennaDefinition::THROW_EXCEPTION);
        values.at(idType).at(i) = antenna->pattern.at(idPattern).antennaVariations(points.at(i).lambda(), points.at(i).phi(), FALSE);
      }

    logStatus<<"save values to file <"<<fileNameGrid<<">"<<Log::endl;
    GriddedData griddedData(Ellipsoid(a, f), points, areas, values);
    writeFileGriddedData(fileNameGrid, griddedData);
    MiscGriddedData::printStatistics(griddedData);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
