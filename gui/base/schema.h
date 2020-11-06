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

  QString         name;
  QString         type;
  QString         annotation;
  QStringList     tags;
  QString         defaultValue;
  Bool            optional;
  Bool            unbounded;
  XsdComplexPtr   complex;
};

/***** CLASS ***********************************/

class XsdComplex
{
public:
  XsdComplex() : isReady(false) {}

  QString                    type;
  std::vector<XsdElementPtr> element;
  Bool                       isReady; // internal: complex pointer is set for all children?
  std::map<QString, QString> renames; // oldName --> newName

  /** @brief Return type for @p name from schema considering renames, or nullptr if not found. */
  XsdElementPtr getXsdElement(QString name) const;
};

/***** CLASS ***********************************/

/**
* @brief XSD (XML schema) tree.
*/
class Schema
{
  std::vector<XsdComplexPtr> complexList;
  std::vector<XsdElementPtr> complexType;

  XsdElementPtr readElement(XmlNodePtr xmlNode);
  XsdComplexPtr readComplex(XmlNodePtr xmlNode);
  void          setComplexPointer(XsdElementPtr element);

public:
  XsdElementPtr rootElement;

  Schema() {}
  Schema(const QString &fileName);

  std::vector<XsdElementPtr> programList() const;
  QString programType() const;

  /** Returns true if @p fileName is a valid GROOPS XML schema, false otherweise. */
  static bool validateSchema(QString fileName);
};

/***********************************************/

#endif
