/***********************************************/
/**
* @file generateDocumentation.h
*
* @brief Generates a latex user documentation file.
*
* @author Torsten Mayer-Guerr
* @author Andreas Kvas
* @date 2008-10-14
*/
/***********************************************/

#ifndef __GROOPS_GENERATEDOCUMENTATION__
#define __GROOPS_GENERATEDOCUMENTATION__

#include "inputOutput/fileName.h"
#include "config/config.h"

/***** CLASS ***********************************/

/** @brief Generates a user documentation.
* @ingroup configGroup */
class Documentation
{
public:
  virtual ~Documentation() {}

  virtual void writeText(const std::string &text) = 0;
  virtual void writeConfigTable(Config &config) = 0;

  /** @brief Generates a hmtl user documentation files.
  * @ingroup configGroup */
  static void writeHtml(const FileName &directoryName);

  /** @brief Generates a latex user documentation files.
  * @ingroup configGroup */
  static void writeLatex(const FileName &directoryName);

  /** @brief Generates documentation files.
  * @ingroup configGroup */
  static void write(const FileName &directoryName);

protected:
  class Token;
  class TableLine;

  static std::vector<Token> tokenize(const std::string &text);
  static void generateTableLine(XmlNodePtr xmlNode, UInt depth, std::vector<TableLine> &tableLines);
  static std::vector<TableLine> generateTable(Config &config);
  static Bool isClassType(const std::string &type);
};

/***********************************************/

#endif /* __GROOPS__ */
