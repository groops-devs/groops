/***********************************************/
/**
* @file fileTideGeneratingPotential.h
*
* @brief Read/write tide genearting potential.
*
* @author Torsten Mayer-Guerr
* @date 2005-01-14
*
*/
/***********************************************/

#ifndef __GROOPS_FILETIDEGENERATINGPOTENTIAL__
#define __GROOPS_FILETIDEGENERATINGPOTENTIAL__

// Latex documentation
#ifdef DOCSTRING_FILEFORMAT_TideGeneratingPotential
static const char *docstringTideGeneratingPotential = R"(

\begin{verbatim}
groops tideGeneratingPotential version=20200123
       7160
# Degree    Dood.    cos                       sin                      name
# ==========================================================================
          3 055.556  0.000000000000000000e+00 -6.569000000000000125e-07 ""
          3 055.635  0.000000000000000000e+00 -2.842000000000000170e-07 ""
          3 055.561  0.000000000000000000e+00 -6.360000000000000040e-08 ""
          2 055.563 -3.122600001000726621e-06  0.000000000000000000e+00 ""
          2 055.565  7.719644799947265879e-02  0.000000000000000000e+00 om1
          3 055.645  0.000000000000000000e+00  2.921429999971616515e-05 ""
          2 055.573  1.975999999999999959e-07  0.000000000000000000e+00 ""
          2 055.575 -7.535264999889109729e-04  0.000000000000000000e+00 om2
\end{verbatim}

)";
#endif

/***********************************************/

#include "base/exception.h"
#include "base/tideGeneratingPotential.h"
#include "inputOutput/fileArchive.h"

/** @addtogroup filesGroup */
/// @{

/***** CONSTANTS ********************************/

const char *const FILE_TIDEGENERATINGPOTENTIAL_TYPE    = "tideGeneratingPotential";
constexpr UInt    FILE_TIDEGENERATINGPOTENTIAL_VERSION = std::max(UInt(20241108), FILE_BASE_VERSION);

/***** FUNCTIONS *******************************/

template<> void save(OutArchive &ar, const TideGeneratingPotential &x);
template<> void load(InArchive  &ar, TideGeneratingPotential &x);

/** @brief Write into a TideGeneratingPotential file. */
void writeFileTideGeneratingPotential(const FileName &fileName, const TideGeneratingPotential &x);

/** @brief Read from a TideGeneratingPotential file. */
void readFileTideGeneratingPotential(const FileName &fileName, TideGeneratingPotential &x);

/// @}

/***********************************************/

#endif /* __GROOPS_FILEGENERATINGTIDE__ */
