/***********************************************/
/**
* @file matrixRectangular2GriddedData.cpp
*
* @brief Read gridded data (matrix).
*
* @author Torsten Mayer-Guerr
* @date 2013-11-16
*
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
Read gridded data (matrix).
)";

/***********************************************/

#include "programs/program.h"
#include "files/fileMatrix.h"
#include "files/fileGriddedData.h"
#include "misc/miscGriddedData.h"

/***** CLASS ***********************************/

/** @brief Read gridded data (matrix).
* @ingroup programsGroup */
class MatrixRectangular2GriddedData
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(MatrixRectangular2GriddedData, SINGLEPROCESS, "read gridded data (matrix)", Grid, Matrix)
GROOPS_RENAMED_PROGRAM(AsciiMatrix2GridRectangular, MatrixRectangular2GriddedData, date2time(2020, 2, 2))

/***********************************************/

void MatrixRectangular2GriddedData::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName fileNameIn, fileNameOut;
    Angle    startL, startB, deltaL, deltaB;
    Double   a, f;
    Bool     rowMajor;

    readConfig(config, "outputfileGriddedData", fileNameOut,   Config::MUSTSET,  "",    "");
    readConfig(config, "inputfileMatrix",       fileNameIn,    Config::MUSTSET,  "",    "");
    readConfig(config, "rowMajor",              rowMajor,      Config::DEFAULT,  "1",   "true: data is ordered row by row, false: columnwise");
    readConfig(config, "startLongitude",        startL,        Config::MUSTSET,  "",    "longitude of upper left corner of the grid");
    readConfig(config, "startLatitude",         startB,        Config::MUSTSET,  "",    "latitude of upper left corner of the grid");
    readConfig(config, "deltaLongitude",        deltaL,        Config::MUSTSET,  "",    "sampling, negative for east to west data");
    readConfig(config, "deltaLatitude",         deltaB,        Config::MUSTSET,  "",    "sampling, negative for south to north data");
    readConfig(config, "R",                     a,             Config::DEFAULT,  STRING_DEFAULT_GRS80_a, "reference radius for ellipsoidal coordinates");
    readConfig(config, "inverseFlattening",     f,             Config::DEFAULT,  STRING_DEFAULT_GRS80_f, "reference flattening for ellipsoidal coordinates");
    if(isCreateSchema(config)) return;

    logStatus<<"Read matrix file <"<<fileNameIn<<">"<<Log::endl;
    Matrix A;
    readFileMatrix(fileNameIn, A);
    logInfo<<"  matrix: "<<A.rows()<<" x "<<A.columns()<<Log::endl;

    // order of data
    // -------------
    const Bool top2Down   = (deltaB>0);
    const Bool left2Right = (deltaL>0);

    // Create rectangular grid
    // -----------------------
    GriddedDataRectangular grid;
    grid.ellipsoid = Ellipsoid(a,f);
    grid.longitudes.resize(A.columns());
    grid.latitudes.resize(A.rows());
    grid.heights.resize(A.rows(), 0.);
    grid.values.resize(1, Matrix(A.rows(), A.columns()));
    for(UInt s=0; s<A.columns(); s++)
      grid.longitudes.at((left2Right) ? s : (A.columns()-s-1)) = Angle(startL+s*Double(deltaL));
    for(UInt z=0; z<A.rows(); z++)
      grid.latitudes.at(((top2Down) ? z : (A.rows()-z-1))) = Angle(startB-z*Double(deltaB));

    // assign data
    // -----------
    for(UInt z=0; z<A.rows(); z++)
      for(UInt s=0; s<A.columns(); s++)
      {
        // determine index
        UInt z1 = ((top2Down)   ? z : (A.rows()-z-1));
        UInt s1 = ((left2Right) ? s : (A.columns()-s-1));
        if(!rowMajor)
          std::swap(z1, s1);
        grid.values.at(0)(z1,s1) = A(z,s);
      }

    MiscGriddedData::printStatistics(grid);

    logStatus<<"write gridded data to file <"<<fileNameOut<<">"<<Log::endl;
    writeFileGriddedData(fileNameOut, grid);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
