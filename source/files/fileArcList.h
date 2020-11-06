/***********************************************/
/**
* @file fileArcList.h
*
* @brief Read/write Arc list.
*
* @author Torsten Mayer-Guerr
* @date 2013-02-07
*
*/
/***********************************************/

#ifndef __GROOPS_FILEARCLIST__
#define __GROOPS_FILEARCLIST__

// Latex documentation
#ifdef DOCSTRING_FILEFORMAT_ArcList
static const char *docstringArcList = R"(
With the \program{InstrumentSynchronize} an \file{instrument file}{instrument} can
be divided into time intervals and within the intervals into arcs.
This file provides the information about the mapping of arcs to time intervals.

This file can be used for the variational equation approach or \program{KalmanBuildNormals}.

\begin{verbatim}
groops arclist version=20200123
         32  # number of times
# time [MJD]               first arc
# ==================================
 58909.000000000000000000          0
 58910.000000000000000000          8
 58911.000000000000000000         17
 58912.000000000000000000         25
 58913.000000000000000000         29
 58914.000000000000000000         37
 58915.000000000000000000         45
 58916.000000000000000000         53
 58917.000000000000000000         61
 58918.000000000000000000         69
 58919.000000000000000000         78
 58920.000000000000000000         86
 58921.000000000000000000         95
 58922.000000000000000000        103
 58923.000000000000000000        112
 58924.000000000000000000        120
 58925.000000000000000000        128
 58926.000000000000000000        136
 58927.000000000000000000        144
 58928.000000000000000000        153
 58929.000000000000000000        161
 58930.000000000000000000        169
 58931.000000000000000000        177
 58932.000000000000000000        185
 58933.000000000000000000        193
 58934.000000000000000000        201
 58935.000000000000000000        210
 58936.000000000000000000        218
 58937.000000000000000000        226
 58938.000000000000000000        234
 58939.000000000000000000        242
 58940.000000000000000000        250
\end{verbatim}
)";
#endif

/***********************************************/

#include "base/exception.h"
#include "base/time.h"
#include "inputOutput/fileArchive.h"

/** @addtogroup filesGroup */
/// @{

/***** CONSTANTS ********************************/

const char *const FILE_ARCLIST_TYPE = "arclist";

/***** FUNCTIONS *******************************/

/** @brief Write into a ArcList file. */
void writeFileArcList(const FileName &fileName, const std::vector<UInt> &arcsInterval, const std::vector<Time> &timesInterval);

/** @brief Read from a ArcList file. */
void readFileArcList(const FileName &fileName, std::vector<UInt> &arcsInterval, std::vector<Time> &timesInterval);

/// @}

/***********************************************/

#endif /* __GROOPS__ */
