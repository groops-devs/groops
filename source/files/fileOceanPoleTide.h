/***********************************************/
/**
* @file fileOceanPoleTide.h
*
* @brief Read/write OceanPoleTide.
*
* @author Torsten Mayer-Guerr
* @date 2017-12-09
*
*/
/***********************************************/

#ifndef __GROOPS_FILEOCEANPOLETIDE__
#define __GROOPS_FILEOCEANPOLETIDE__

// Latex documentation
#ifdef DOCSTRING_FILEFORMAT_OceanPoleTide
static const char *docstringOceanPoleTide = R"(
Describes the reaction of the ocean mass to the change
of the centrifugal potential (polar wobble) in terms spherical harmonics.

See also \program{Iers2OceanPoleTide}.
)";
#endif

/***********************************************/

#include "base/sphericalHarmonics.h"
#include "inputOutput/fileArchive.h"

/** @addtogroup filesGroup */
/// @{

/***** CONSTANTS ********************************/

const char *const FILE_OCEANPOLETIDE_TYPE = "oceanPoleTide";

/***** FUNCTIONS *******************************/

/** @brief Write into a OceanPoleTide file. */
void writeFileOceanPoleTide(const FileName &fileName, const SphericalHarmonics &harmReal, const SphericalHarmonics &harmImag);


/** @brief Read from a OceanPoleTide file. */
void readFileOceanPoleTide(const FileName &fileName, SphericalHarmonics &harmReal, SphericalHarmonics &harmImag);

/// @}

/***********************************************/

#endif /* __GROOPS__ */
