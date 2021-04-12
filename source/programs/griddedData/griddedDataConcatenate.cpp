/***********************************************/
/**
* @file griddedDataConcatenate.cpp
*
* @brief Concatenate gridded data from several files.
*
* @author Torsten Mayer-Guerr
* @date 2020-02-08
*
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program concatenate grid from several \configFile{inputfileGriddedData}{griddedData}
and write it to a new \configFile{outputfileGriddedData}{griddedData}.
Input files must have the same number of data columns.
If \config{sort} is enabled, the points are sorted by latitudes starting from north/west to south east.
Identical points (within a \config{margin}) can be removed with \config{removeDuplicates}.
)";

/***********************************************/

#include "programs/program.h"
#include "files/fileGriddedData.h"
#include "misc/miscGriddedData.h"
#include "classes/border/border.h"

/***********************************************/

/** @brief Concatenate gridded data from several files.
* @ingroup programsGroup */
class GriddedDataConcatenate
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(GriddedDataConcatenate, SINGLEPROCESS, "Concatenate gridded data from several files", Grid)

/***********************************************/

void GriddedDataConcatenate::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName              fileNameOut, fileNameStatistics;
    std::vector<FileName> fileNamesIn;
    BorderPtr             border;
    Bool                  sortPoints;
    Double                a, f;
    std::string           choice;
    enum Keep {KEEPALL, KEEPFIRST, KEEPLAST};
    Keep                  keep = KEEPALL;
    Double                margin;

    readConfig(config, "outputfileGriddedData", fileNameOut,      Config::MUSTSET,  "",  "");
    readConfig(config, "inputfileGriddedData",  fileNamesIn,      Config::MUSTSET,  "",  "");
    readConfig(config, "border",                border,           Config::OPTIONAL, "",  "");
    readConfig(config, "sortPoints",            sortPoints,       Config::DEFAULT,  "0", "sort from north/west to south east");
    if(readConfigChoice(config, "removeDuplicates", choice, Config::OPTIONAL, "", "remove duplicate points"))
    {
      if(readConfigChoiceElement(config, "keepFirst", choice, "keep first point, remove all other identicals"))
      {
        readConfig(config, "margin", margin, Config::DEFAULT, "1e-5", "margin distance for identical points [m]");
        keep = KEEPFIRST;
      }
      if(readConfigChoiceElement(config, "keepLast",  choice, "keep last point, remove all other identicals"))
      {
        readConfig(config, "margin", margin, Config::DEFAULT, "1e-5", "margin distance for identical points [m]");
        keep = KEEPLAST;
      }
      endChoice(config);
    }
    readConfig(config, "R",                 a, Config::DEFAULT, STRING_DEFAULT_GRS80_a, "reference radius for ellipsoidal coordinates");
    readConfig(config, "inverseFlattening", f, Config::DEFAULT, STRING_DEFAULT_GRS80_f, "reference flattening for ellipsoidal coordinates");
    if(isCreateSchema(config)) return;

    // =====================================================

    // read data
    // ---------
    GriddedData grid;
    for(auto fileName : fileNamesIn)
    {
      GriddedData gridFile;
      try
      {
        logStatus<<"reading grid from file <"<<fileName<<">"<<Log::endl;
        readFileGriddedData(fileName, gridFile);
      }
      catch(std::exception &e)
      {
        logWarning<<e.what()<<" continue..."<<Log::endl;
        continue;
      }

      if(grid.points.size() == 0)
      {
        grid = gridFile;
        grid.ellipsoid = Ellipsoid(a,f);
        continue;
      }

      if((grid.areas.size()==0) != (gridFile.areas.size()==0))
        throw(Exception("all grids must have areas or none of them"));
      if(grid.values.size() != gridFile.values.size())
        throw(Exception("all grids must agree in the number of data columns"));

      grid.points.insert(grid.points.end(), gridFile.points.begin(), gridFile.points.end());
      grid.areas.insert (grid.areas.end(),  gridFile.areas.begin(),  gridFile.areas.end());
      for(UInt k=0; k<grid.values.size(); k++)
        grid.values.at(k).insert(grid.values.at(k).end(), gridFile.values.at(k).begin(), gridFile.values.at(k).end());
    }

    // =====================================================

    auto movePoints = [&](UInt from, UInt to)
    {
      grid.points.at(to) = grid.points.at(from);
      if(grid.areas.size())
        grid.areas.at(to)  = grid.areas.at(from);
      for(UInt k=0; k<grid.values.size(); k++)
        grid.values.at(k).at(to) = grid.values.at(k).at(from);
    };

    // remove points outside border
    // ----------------------------
    if(border)
    {
      logStatus<<"remove points outside border"<<Log::endl;
      UInt count = std::distance(grid.points.begin(), std::find_if(grid.points.begin(), grid.points.end(), [&](const Vector3d &p) {return !border->isInnerPoint(p, grid.ellipsoid);}));
      for(UInt i=count+1; i<grid.points.size(); i++)
        if(border->isInnerPoint(grid.points.at(i), grid.ellipsoid))
          movePoints(i, count++);
      logInfo<<" "<<grid.points.size()-count<<" points removed!"<<Log::endl;
      grid.points.resize(count);
      if(grid.areas.size())
        grid.areas.resize(count);
      for(UInt k=0; k<grid.values.size(); k++)
        grid.values.at(k).resize(count);
    }

    // =====================================================

    if(sortPoints || (keep != KEEPALL))
    {
      logStatus<<"sort points"<<Log::endl;
      grid.sort();
    }

    // =====================================================

    // eliminate duplicates
    // --------------------
    if(keep != KEEPALL)
    {
      logStatus<<"eliminate duplicates"<<Log::endl;
      UInt count = std::distance(grid.points.begin(), std::adjacent_find(grid.points.begin(), grid.points.end(),
                                 [&](const Vector3d &x, const Vector3d &y) {return (x-y).r() <= margin;}));
      if(count < grid.points.size())
        for(UInt i=++count; i<grid.points.size(); i++)
        {
          if((grid.points.at(i)-grid.points.at(count-1)).r() > margin)
            movePoints(i, count++);
          else if(keep == KEEPLAST)
            movePoints(i, count-1);
        }
      logInfo<<" "<<grid.points.size()-count<<" duplicates removed!"<<Log::endl;
      grid.points.resize(count);
      if(grid.areas.size())
        grid.areas.resize(count);
      for(UInt k=0; k<grid.values.size(); k++)
        grid.values.at(k).resize(count);
    }

    // =====================================================

    // save grid
    // ---------
    logStatus<<"save grid <"<<fileNameOut<<">"<<Log::endl;
    writeFileGriddedData(fileNameOut, grid);
    MiscGriddedData::printStatistics(grid);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
