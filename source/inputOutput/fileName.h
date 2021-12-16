/***********************************************/
/**
* @file fileName.h
*
* @brief File names
*
* @author Torsten Mayer-Guerr
* @date 2008-07-28
*
*/
/***********************************************/

#ifndef __GROOPS_FILENAME__
#define __GROOPS_FILENAME__

#include "base/importStd.h"

class VariableList;

/***** CLASS ***********************************/

/** @brief File names.
* @ingroup inputOutputGroup */
class FileName
{
  std::string name_;

public:
  FileName() {}                                      //!< Constructor.
  FileName(const FileName &) = default;              //!< Copy constructor.
  FileName(const std::string &name) : name_(name) {} //!< Constructor.
  FileName(const char *name) : name_(name) {}        //!< Constructor.
  FileName &operator=(const FileName &) = default;   //!< Assignement.

  operator std::string()  const {return name_;}         //!< Cast to  string.
  operator const char *() const {return name_.c_str();} //!< Cast to old C type string.
  const char *c_str()     const {return name_.c_str();} //!< Old C type string.
  std::string str()       const {return name_;}         //!< String.

  /** @brief Is FileName empty? */
  Bool empty() const {return name_.empty();}

  /** @brief Append a fileName to a directory.
  * Directory and file will be separated by the directory separator (e.g. '/') if needed. */
  FileName append(const FileName& fileName) const;

  /** @brief Extension.
  * Returns the string of all characters in the FileName
  * after (but not including) the last '.' character plus an additional ".gz" or ".z".
  * Example: FileName("name.txt.gz").typeExtension() returns "txt.gz". */
  FileName fullExtension() const;

  /** @brief Extension.
  * Returns the string of all characters in the FileName
  * after (but not including) the last '.' character without an additional ".gz" or ".z".
  * Example: FileName("name.txt.gz").typeExtension() returns "txt". */
  FileName typeExtension() const;

  /** @brief Extension.
  * Example: FileName("name.txt.gz").typeExtension() returns "gz". */
  FileName packExtension() const;

  /** @brief File name without extension.
  * Returns the string of all characters in the FileName
  * before (but not including) the last '.' character. */
  FileName stripFullExtension() const;

  /** @brief Replaces the extension with text.
  * Firstly, if this path has an extension(), it is removed.
  * Then, a dot character is appended if text is not empty or does not begin with a dot character.
  * Then text is appended to the path. */
  FileName replaceFullExtension(const std::string &text) const;

  /** @brief Directory name of file.
  * Returns the string of all characters in the FileName
  * name before the last '/' character. */
  FileName directory() const;

  /** @brief File name without directory.
  * Returns the string of all characters in the FileName
  * name after (but not including) the last '/' character. */
  FileName stripDirectory() const;

  /** @brief File name without directory and extension.
  * Returns (*this).stripDirectory().stripFullExtension(). */
  FileName baseName() const {return stripDirectory().stripFullExtension();}

  /** @brief Extent the baseName.
  * Returns a FileName with appending text before '.'+extension(). */
  FileName appendBaseName(const std::string &text) const;

  /** Returns a FileName with replaced variable {name} names. */
  FileName operator()(const VariableList &varList) const;

  /** @brief 'Less than' operator, required for sorting. */
  Bool operator< (const FileName &fileName) const {return name_ < fileName.name_;}
};

/*************************************************/

#endif /* __GROOPS__ */


