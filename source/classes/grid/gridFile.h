/***********************************************/
/**
* @file gridFile.h
*
* @brief Read grid from file.
* @see Grid
*
* @author Torsten Mayer-Guerr
* @date 2004-01-11
*
*/
/***********************************************/

#ifndef __GROOPS_GRIDFILE__
#define __GROOPS_GRIDFILE__

// Latex documentation
#ifdef DOCSTRING_Grid
static const char *docstringGridFile = R"(
\subsection{File}\label{gridType:file}
In this class grid is read from a file, which is given by \configFile{inputfileGrid}{griddedData}.
A corresponding file can be generated with \program{GriddedDataCreate} or with \program{Matrix2GriddedData}.
)";
#endif

/***********************************************/

#include "base/import.h"
#include "config/config.h"
#include "files/fileGriddedData.h"
#include "classes/border/border.h"
#include "classes/grid/grid.h"

/***** CLASS ***********************************/

/** @brief Read grid from file.
* @ingroup gridGroup
* @see Grid */
class GridFile : public GridBase
{
public:
  GridFile(Config &config);
};

/***********************************************/

inline GridFile::GridFile(Config &config)
{
  try
  {
    FileName  fileNameGrid;
    BorderPtr border;

    readConfig(config, "inputfileGrid",  fileNameGrid, Config::MUSTSET,  "", "");
    readConfig(config, "border",         border,       Config::DEFAULT,  "", "");
    if(isCreateSchema(config)) return;

    GriddedData grid;
    readFileGriddedData(fileNameGrid, grid);

    for(UInt i=0; i<grid.points.size(); i++)
      if(border->isInnerPoint(grid.points.at(i), grid.ellipsoid))
      {
        points.push_back(grid.points.at(i));
        if(grid.areas.size())
          areas.push_back(grid.areas.at(i));
      }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

#endif /* __GROOPS__ */
