/***********************************************/
/**
* @file instrumentType.h
*
* @brief Defines the type of an instrument file.
*
* @author Torsten Mayer-Guerr
* @date 2020-05-29
*
*/
/***********************************************/

#ifndef __GROOPS_INSTRUMENTTYPE__
#define __GROOPS_INSTRUMENTTYPE__

// Latex documentation
#ifdef DOCSTRING_InstrumentType
static const char *docstringInstrumentType = R"(
\section{InstrumentType}\label{instrumentTypeType}
Defines the type of an \file{instrument file}{instrument}.
)";
#endif

/***********************************************/

#include "base/import.h"
#include "config/config.h"
#include "files/fileInstrument.h"

/**
* @defgroup instrumentTypeGroup InstrumentType
* @brief Defines the type of instrument files.
* @ingroup classesGroup
* The interface is given by @ref Epoch::Type.
* An Instance can be created by @ref readConfig. */
/// @{

/***** FUNCTIONS *******************************/

/** @brief Reads an instrument type from config.
* Search for a node with @a name in the Config node.
* if @a name is not found the function returns FALSE and @a var is untouched.
* @param config The config node which includes the node with the options for @a var
* @param name Tag name in the config.
* @param[out] var read from config.
* @param mustSet If is MUSTSET and @a name is not found, this function throws an exception instead of returning with FALSE.
* @param defaultValue Ignored at the moment.
* @param annotation Short description. */
template<> Bool readConfig(Config &config, const std::string &name, Epoch::Type &var, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation);

/***********************************************/

/// @}

#endif /* __GROOPS__ */
