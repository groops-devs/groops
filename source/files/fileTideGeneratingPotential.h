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
       3606
# Dood.   cos                       sin                      name
# ===============================================================
 055.563 -3.122600001000726621e-06  0.000000000000000000e+00 ""
 055.565  7.719644799947265879e-02  0.000000000000000000e+00 om1
 055.573  1.975999999999999959e-07  0.000000000000000000e+00 ""
 055.575 -7.535264999889109729e-04  0.000000000000000000e+00 om2
 055.654 -4.037500000326006771e-06  0.000000000000000000e+00 ""
 055.666  6.590000000000000115e-08  0.000000000000000000e+00 ""
 055.753  6.023000000000000398e-07  0.000000000000000000e+00 ""
 056.475  7.920000000000000775e-08  0.000000000000000000e+00 ""
 056.554 -1.360322439950949897e-02  0.000000000000000000e+00 sa
 056.556  7.156881999672112154e-04  0.000000000000000000e+00 ""
 056.564  8.641769999583327993e-05  0.000000000000000000e+00 ""
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

const char *const FILE_TIDEGENERATINGPOTENTIAL_TYPE = "tideGeneratingPotential";

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
