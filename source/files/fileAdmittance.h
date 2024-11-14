/***********************************************/
/**
* @file fileAdmittance.h
*
* @brief Read/write admittance.
*
* @author Torsten Mayer-Guerr
* @date 2013-02-08
*
*/
/***********************************************/

#ifndef __GROOPS_FILEADMITTANCE__
#define __GROOPS_FILEADMITTANCE__

// Latex documentation
#ifdef DOCSTRING_FILEFORMAT_Admittance
static const char *docstringAdmittance = R"(
Interpolation matrix to create ocean minor tides from modeled major tides.
The file can be created with \program{DoodsonHarmonicsCalculateAdmittance} and used e.g. in
\configClass{doodsonHarmonicTide}{tidesType:doodsonHarmonicTide}.

See \program{DoodsonHarmonicsCalculateAdmittance}.
)";
#endif

/***********************************************/

#include "base/import.h"
#include "base/doodson.h"
#include "inputOutput/fileArchive.h"

/** @addtogroup filesGroup */
/// @{

/***** CONSTANTS ********************************/

const char *const FILE_ADMITTANCE_TYPE    = "admittance";
constexpr UInt    FILE_ADMITTANCE_VERSION = std::max(UInt(20200123), FILE_BASE_VERSION);

/***** CLASS ***********************************/

/** @brief Admittance to interpolate minor tides from major tides. */
class Admittance
{
public:
  /// List of major tides.
  std::vector<Doodson> doodsonMajor;

  /// List of minor tides (inclusive major tides).
  std::vector<Doodson> doodsonMinor;

  /** @brief Interpolation matrix.
  Dimension: major count times minor count. */
  Matrix admittance;
};

/***** FUNCTIONS *******************************/

/** @brief Write into a Admittance file. */
void writeFileAdmittance(const FileName &fileName, const Admittance &admittance);

/** @brief Read from a Admittance file. */
void readFileAdmittance(const FileName &fileName, Admittance &admittance);

/// @}

/***********************************************/

#endif /* __GROOPS__ */
