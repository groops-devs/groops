/***********************************************/
/**
* @file griddedData2GriddedDataStatistics.cpp
*
* @brief Assign gridded data to grid points.
*
* @author Torsten Mayer-Guerr
* @date 2018-07-03
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program assigns values \configFile{inputfileGriddedData}{griddedData} to the nearest points
of a new \configClass{grid}{gridType}. If some of the new points are not filled in with data
\config{emptyValue} is used instead. If multiple points of the input fall on the same node
the result can be selected with \config{statistics} (e.g. mean, root mean square, min, max, \ldots).
It also is possible to simply count the number of data points that were assigned to each point.

Be aware in case borders are given within \configClass{grid}{gridType}, the \configFile{outputfileGriddedData}{griddedData} will have points excluded before the assignement of old points to the new points.
The data from \configFile{inputfileGriddedData}{griddedData} will not be limited by the given borders! See \reference{GriddedDataConcatenate}{GriddedDataConcatenate} to limit the
\configFile{inputfileGriddedData}{griddedData} to given borders.

\fig{!hb}{0.8}{griddedData2GriddedDataStatistics}{fig:griddedData2GriddedDataStatistics}{Assignement of irregular distributed data to grid.}
)";

/***********************************************/

#include "programs/program.h"
#include "files/fileGriddedData.h"
#include "classes/grid/grid.h"
#include "misc/miscGriddedData.h"

/***** CLASS ***********************************/

/** @brief Assign gridded data to grid points.
* @ingroup programsGroup */
class GriddedData2GriddedDataStatistics
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(GriddedData2GriddedDataStatistics, SINGLEPROCESS, "Assign gridded data to grid points", Grid)

/***********************************************/

void GriddedData2GriddedDataStatistics::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    enum Type {MEAN, WMEAN, RMS, WRMS, STD, WSTD, SUM, MIN, MAX, COUNT, FIRST, LAST};
    Type        type = MEAN;
    FileName    fileNameOutGrid, fileNameInGrid;
    GridPtr     gridPtr;
    Double      emptyValue = NAN_EXPR;
    Double      a, f;
    std::string choice;

    readConfig(config, "outputfileGriddedData", fileNameOutGrid, Config::MUSTSET,  "",    "");
    readConfig(config, "inputfileGriddedData",  fileNameInGrid,  Config::MUSTSET,  "",    "");
    readConfig(config, "grid",                  gridPtr,         Config::MUSTSET,  "",    "");
    if(readConfigChoice(config, "statistic", choice, Config::MUSTSET, "mean", "statistic used if multiple values fall on the same cell"))
    {
      if(readConfigChoiceElement(config, "mean",  choice, "mean"))                             type = MEAN;
      if(readConfigChoiceElement(config, "wmean", choice, "area weighted mean"))               type = WMEAN;
      if(readConfigChoiceElement(config, "rms",   choice, "root mean square"))                 type = RMS;
      if(readConfigChoiceElement(config, "wrms",  choice, "area weighted root mean square"))   type = WRMS;
      if(readConfigChoiceElement(config, "std",   choice, "standard deviation"))               type = STD;
      if(readConfigChoiceElement(config, "wstd",  choice, "area weighted standard deviation")) type = WSTD;
      if(readConfigChoiceElement(config, "sum",   choice, "sum"))                              type = SUM;
      if(readConfigChoiceElement(config, "min",   choice, "minimum value"))                    type = MIN;
      if(readConfigChoiceElement(config, "max",   choice, "maximum value"))                    type = MAX;
      if(readConfigChoiceElement(config, "count", choice, "number of values"))                 type = COUNT;
      if(readConfigChoiceElement(config, "first", choice, "first value"))                      type = FIRST;
      if(readConfigChoiceElement(config, "last",  choice, "last value"))                       type = LAST;
    }
    endChoice(config);
    readConfig(config, "emptyValue",        emptyValue, Config::DEFAULT, "nan()", "value for nodes without data");
    readConfig(config, "R",                 a,          Config::DEFAULT, STRING_DEFAULT_GRS80_a, "reference radius for ellipsoidal coordinates");
    readConfig(config, "inverseFlattening", f,          Config::DEFAULT, STRING_DEFAULT_GRS80_f, "reference flattening for ellipsoidal coordinates");
    if(isCreateSchema(config)) return;

    // read grid
    // ---------
    logStatus<<"read grid from file <"<<fileNameInGrid<<">"<<Log::endl;
    GriddedData grid;
    readFileGriddedData(fileNameInGrid, grid);
    MiscGriddedData::printStatistics(grid);
    if((grid.areas.size() == 0) && ((type == WMEAN) || (type == WRMS) || (type == WSTD)))
      throw(Exception("no area elements"));

    // Create grid
    // -----------
    logStatus<<"create new grid"<<Log::endl;
    GriddedData gridNew(Ellipsoid(a, f), gridPtr->points(), gridPtr->areas(), {});
    Double initialValue = 0;
    if(type == MIN) initialValue =  1e99;
    if(type == MAX) initialValue = -1e99;
    gridNew.values.resize(grid.values.size(), std::vector<Double>(gridNew.points.size(), initialValue));

    std::vector<Angle>  lambda, phi;
    std::vector<Double> radius;
    const Bool isRectangle = gridNew.isRectangle(lambda, phi, radius);

    // additional variables
    std::vector<std::vector<Double>> count, wmean, weight;
    count.resize(gridNew.values.size(), std::vector<Double>(gridNew.points.size(), 0));
    if((type == STD) || (type == WSTD))
      wmean.resize(gridNew.values.size(), std::vector<Double>(gridNew.points.size(), 0));
    if((type == WMEAN) || (type == WRMS) || (type == WSTD))
      weight.resize(gridNew.values.size(), std::vector<Double>(gridNew.points.size(), 0));

    // Assign grid
    // -----------
    logStatus<<"assign grid"<<Log::endl;
    Single::forEach(grid.points.size(), [&](UInt i)
    {
      // find nearest neighbor
      UInt idx = 0;
      if(isRectangle)
      {
        const Angle lon = grid.points.at(i).lambda();
        const Angle lat = grid.points.at(i).phi();
        const UInt  col = std::distance(lambda.begin(), std::min_element(lambda.begin(), lambda.end(),
                                        [&](Angle lon1, Angle lon2) {return std::fabs(std::remainder(lon-lon1, 2*PI)) < std::fabs(std::remainder(lon-lon2, 2*PI));}));
        const UInt  row = std::distance(phi.begin(), std::min_element(phi.begin(), phi.end(),
                                        [&](Angle lat1, Angle lat2) {return std::fabs(lat-lat1) < std::fabs(lat-lat2);}));
        idx = row * lambda.size() + col;
      }
      else
      {
        Double minDistance = std::numeric_limits<Double>::max();
        for(UInt k=0; k<gridNew.points.size(); k++)
        {
          const Double distance = (gridNew.points[k] - grid.points[i]).r();
          if(distance < minDistance)
          {
            minDistance = distance;
            idx = k;
          }
        }
      }
      Double w = 1;
      if((type == WMEAN) || (type == WRMS) || (type == WSTD))
        w = grid.areas.at(i);

      for(UInt k=0; k<grid.values.size(); k++)
      {
        const Double v = grid.values.at(k).at(i);
        if(std::isnan(v))
          continue;

        if(count.size())  count.at(k).at(idx)++;
        if(weight.size()) weight.at(k).at(idx) += w;
        if(wmean.size())  wmean.at(k).at(idx)  += w*v;

        switch(type)
        {
          case MEAN:
          case WMEAN:
          case SUM:   gridNew.values.at(k).at(idx) += w*v;   break;
          case RMS:
          case WRMS:
          case STD:
          case WSTD:  gridNew.values.at(k).at(idx) += w*v*v; break;
          case MIN:   gridNew.values.at(k).at(idx)  = std::min(v, gridNew.values.at(k).at(idx)); break;
          case MAX:   gridNew.values.at(k).at(idx)  = std::max(v, gridNew.values.at(k).at(idx)); break;
          case LAST:  gridNew.values.at(k).at(idx)  = v  ; break;
          case FIRST: if(count.at(k).at(idx) == 1) gridNew.values.at(k).at(idx) = v; break;
          default: ;
        }
      }
    });

    // post computation
    // ----------------
    for(UInt k=0; k<gridNew.values.size(); k++)
      for(UInt idx=0; idx<gridNew.values.at(k).size(); idx++)
        if(count.at(k).at(idx))
        {
          switch(type)
          {
            case MEAN:  gridNew.values.at(k).at(idx) /= count.at(k).at(idx);  break;
            case WMEAN: gridNew.values.at(k).at(idx) /= weight.at(k).at(idx); break;
            case RMS:   gridNew.values.at(k).at(idx) = std::sqrt(gridNew.values.at(k).at(idx)/count.at(k).at(idx));  break;
            case WRMS:  gridNew.values.at(k).at(idx) = std::sqrt(gridNew.values.at(k).at(idx)/weight.at(k).at(idx)); break;
            case STD:   wmean.at(k).at(idx) /= count.at(k).at(idx);
                        gridNew.values.at(k).at(idx) = std::sqrt(count.at(k).at(idx)/(count.at(k).at(idx)-1.)*(gridNew.values.at(k).at(idx)/count.at(k).at(idx)-std::pow(wmean.at(k).at(idx), 2))); break;
            case WSTD:  wmean.at(k).at(idx) /= weight.at(k).at(idx);
                        gridNew.values.at(k).at(idx) = std::sqrt(count.at(k).at(idx)/(count.at(k).at(idx)-1.)*(gridNew.values.at(k).at(idx)/weight.at(k).at(idx)-std::pow(wmean.at(k).at(idx), 2))); break;
            case COUNT: gridNew.values.at(k).at(idx) = count.at(k).at(idx); break;
            default: ;
          }
        }
        else
          gridNew.values.at(k).at(idx) = emptyValue;

    // Output
    // ------
    logStatus<<"write grid to file <"<<fileNameOutGrid<<">"<<Log::endl;
    writeFileGriddedData(fileNameOutGrid, gridNew);
    MiscGriddedData::printStatistics(gridNew);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
