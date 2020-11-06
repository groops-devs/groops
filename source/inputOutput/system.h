/***********************************************/
/**
* @file system.h
*
* @brief Operating system related functions.
*
* @author Andreas Kvas
* @date 2017-03-21
*
*/
/***********************************************/

#ifndef __GROOPS_SYSTEM__
#define __GROOPS_SYSTEM__

#include "base/importStd.h"
#include "inputOutput/fileName.h"

/** @brief System operations.
* @ingroup inputOutputGroupGroup */
namespace System
{
  /** @brief Execute a command using the system shell.
  * @param command command to be executed as string.
  * @param[out] output command output is appended line by line (without newline).
  * @return success if executed command returns without error.  */
  Bool exec(const std::string &command, std::vector<std::string> &output);

  /** @brief Execute a command using the system shell.
  * @param command command to be executed as string
  * @return success if executed command returns without error.  */
  Bool exec(const std::string &command);

  /** @brief Creates the directory and parent directories as needed.
  * @param fileName file to be created.
  * @return TRUE if a directory was created or exists already. */
  Bool createDirectories(const FileName &fileName);

  /** @brief Remove a file or directory.
  * Deletes also the content recursivley if @a fileName is a directory.
  * @param fileName file to be removed.
  * @return TRUE if the file was deleted. */
  Bool remove(const FileName &fileName);

  /** @brief Checks if the given file or path corresponds to an existing file or directory. */
  Bool exists(const FileName &fileName);

  /** @brief Check wether fileName is an existing directory */
  Bool isDirectory(const FileName &fileName);

  /** @brief Current working directory as FileName. */
  FileName currentWorkingDirectory();

  /** @brief Current time as used by the file system. */
  Time now();
}

/***********************************************/

#endif

