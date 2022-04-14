/***********************************************/
/**
* @file dataVariables.h
*
* @brief Create variables to evaluate a matrix/grid/timeSeries.
*
* @author Torsten Mayer-Guerr
* @date 2018-06-18
*
*/
/***********************************************/

#ifndef __GROOPS_DATAVARIABLES__
#define __GROOPS_DATAVARIABLES__

// Latex documentation
#ifdef DOCSTRING_Parser
static const char *docstringParserDataVariables = R"(
\subsection{Variables for data}\label{general.parser:dataVariables}
Some programs (e.g. \program{FunctionsCalculate}, \program{InstrumentArcCalculate},
\program{GriddedDataCalculate}, or the plot programs)
read data (\file{matrix}{matrix}) or \file{gridded data}{griddedData}
and evaluate input/output expressions for each data row.
For these kind of expressions additional variables are automatically defined for each data column
(\verb|X| stands for the data column number: $0\ldots n$):
\begin{itemize}
\item \verb|index|: the row number, starting with zero
\item \verb|dataX|: the value itself
\item \verb|dataXcount|: number of rows
\item \verb|dataXmin|
\item \verb|dataXmax|
\item \verb|dataXsum|
\item \verb|dataXmean|
\item \verb|dataXrms|: root mean square
\item \verb|dataXstd|: standard deviation
\item \verb|dataXmedian|
\item \verb|dataXmad|: median absolute deviation
\item \verb|dataXstep|: the minimal difference between two neighboring data points in the column
\end{itemize}
For \file{gridded data}{griddedData} input the following variables are additionally defined for each data point:
\begin{itemize}
\item \verb|longitude| in degrees
\item \verb|latitude| in degrees
\item \verb|height| in meters
\item \verb|cartesianX| coordinate in meters
\item \verb|cartesianY| coordinate in meters
\item \verb|cartesianZ| coordinate in meters
\item \verb|area| of the unit sphere
\item \verb|dataXwmean|: area-weighted mean
\item \verb|dataXwrms|: area-weighted root mean square
\item \verb|dataXwstd|: area-weighted standard deviation
\end{itemize}
)";
#endif

/***********************************************/

#include "base/importStd.h"
#include "base/matrix.h"
#include "expressionParser.h"

/** @addtogroup parserGroup */
/// @{

/***********************************************/

/** @brief Create variables to loop over intervals.
* @ingroup parserGroup
* The following variable are created:
* index, loopTime, loopTimeStart, loopTimeEnd. */
void addTimeVariables(VariableList &varList);

/** @brief Compute the values of the variables of a temporal loop.
* @ingroup parserGroup */
void evaluateTimeVariables(UInt index, const Time &timeStart, const Time &timeEnd, VariableList &varList);

/***********************************************/

/** @brief Create variables to evaluate a vector.
* @ingroup parserGroup
* The following variable are created if they occur in @a usedName
* prefix, prefixmean, prefixmedian, prefixrms, prefixstd, prefixmad, prefixmin, prefixmax, prefixstep
* @param prefix to build parameter names.
* @param data The vector for which the variables are created.
* @param varList The @a VariableList to which the variables are added.
* @param usedName Only variables are created which are in the list */
void addDataVariables(const std::string &prefix, const_MatrixSliceRef data, VariableList &varList, const std::set<std::string> &usedName);

class Time;
void addDataVariables(const std::string &prefix, const std::vector<Time> &times, VariableList &varList, const std::set<std::string> &usedName);

void addDataVariables(const std::string &prefix, const_MatrixSliceRef data, const_MatrixSliceRef weight, VariableList &varList, const std::set<std::string> &usedName);

/** @brief Create variables to evaluate a matrix.
* @ingroup parserGroup
* The following variable are created if they occur in @a usedName
* index,
* data, data0, data1, ...
* data0mean, data0median, data0rms, data0std, data0mad, data0min, data0max, data0step
* data1mean, data1median, data1rms, data1std, data1mad, data1min, data1max, data1step
* ...
* @param data The matrix for which the variables are created.
* @param varList The @a VariableList to which the variables are added.
* @param usedName Only variables are created which are in the list */
void addDataVariables(const_MatrixSliceRef data, VariableList &varList, const std::set<std::string> &usedName);

/** @brief Compute the values of the variables of a matrix.
* @ingroup parserGroup */
void evaluateDataVariables(const_MatrixSliceRef data, UInt row, VariableList &varList);

/** @brief Set status to undefined of the variables of a matrix.
* @ingroup parserGroup */
void undefineDataVariables(const_MatrixSliceRef data, VariableList &varList);

/***********************************************/

class GriddedData;
class GriddedDataRectangular;
/** @brief Create variables to evaluate a grid.
* @ingroup parserGroup
* The following variable are created
* longitude, latitude, height, area,
* index,
* data, data0, data1, ...
* data0mean, data0rms, data0std, data0min, data0max,
* data1mean, data1rms, data1std, data1min, data1max, ... */
void addDataVariables(const GriddedData            &grid, VariableList &varList, const std::set<std::string> &usedName);
void addDataVariables(const GriddedDataRectangular &grid, VariableList &varList, const std::set<std::string> &usedName);

/** @brief Compute the values of the variables of a grid.
* @ingroup parserGroup */
void evaluateDataVariables(const GriddedData            &grid, UInt row, VariableList &varList);
void evaluateDataVariables(const GriddedDataRectangular &grid, UInt row, UInt col, VariableList &varList);

/** @brief Set status to undefined of the variables of a grid.
* @ingroup parserGroup */
void undefineDataVariables(const GriddedData            &grid, VariableList &varList);
void undefineDataVariables(const GriddedDataRectangular &grid, VariableList &varList);

/***********************************************/

/// @}

#endif /* __GROOPS__ */

