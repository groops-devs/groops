/***********************************************/
/**
* @file schema.h
*
* @brief XSD (XML schema) tree.
*
* @author Torsten Mayer-Guerr
* @author Sebastian Strasser
* @date 2006-06-10
*/
/***********************************************/

#ifndef __GROOPSGUI__SCHEMA__
#define __GROOPSGUI__SCHEMA__

#include <QStringList>
#include "base/importGroops.h"
#include "base/xml.h"

/***** TYPES ***********************************/

class   XsdElement;
class   XsdComplex;
typedef std::shared_ptr<XsdElement> XsdElementPtr;
typedef std::shared_ptr<XsdComplex> XsdComplexPtr;

/***** CLASS ***********************************/

class XsdElement
{
public:
  XsdElement() : optional(true), unbounded(false), complex(nullptr) {}

  QStringList   names; // with possible renames
  QString       type;
  QString       annotation;
  QString       defaultValue;
  QStringList   tags;
  bool          optional;
  bool          unbounded;
  XsdComplexPtr complex;
};

/***** CLASS ***********************************/

class XsdComplex
{
public:
  XsdComplex() : isReady(false) {}

  QString                    type;    // "choice" or "sequence"
  std::vector<XsdElementPtr> elements;
  bool                       isReady; // internal: complex pointer is set for all children?

  /** @brief Return type for @p name from schema considering renames, or nullptr if not found. */
  XsdElementPtr getXsdElement(QString name) const;
};

/***** CLASS ***********************************/

/** @brief XSD (XML schema) tree. */
class Schema
{
  XsdElementPtr readElement(XmlNodePtr xmlNode, const std::map<QString, QString> &renames, bool isComplexType);
  XsdComplexPtr readComplex(XmlNodePtr xmlNode, const std::map<QString, QString> &renames);
  void          setComplexPointer(XsdElementPtr element);

public:
  std::vector<XsdElementPtr> complexType;
  XsdElementPtr rootElement;

  /** @brief reads a GROOPS XML schema.
  * @return is a valid GROOPS XML schema. */
  bool readFile(QString fileName);
};

/***********************************************/

#endif
