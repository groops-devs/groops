/***********************************************/
/**
* @file settings.h
*
* @brief Read/write constants.
*
* @author Torsten Mayer-Guerr
* @date 2010-03-01
*
*/
/***********************************************/

#ifndef __GROOPS_SETTINGS__
#define __GROOPS_SETTINGS__

#include "inputOutput/fileName.h"

/***** FUNCTIONS *******************************/

/// @{

/** @brief Read constants to file.
* @ingroup inputOutputGroup
* The default values are overwritten. */
void readFileSettings(const FileName &fileName);

/** @brief Write the currently used constants to file.
* @ingroup inputOutputGroup
@verbatim
groops --write-settings <fileName.xml>
@endverbatim */
void writeFileSettings(const FileName &fileName);

/// @}

/***********************************************/

#endif /* __GROOPS__ */
