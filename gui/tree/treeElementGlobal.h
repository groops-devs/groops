/***********************************************/
/**
* @file treeElementGlobal.h
*
* @brief The global element.
*
* @author Torsten Mayer-Guerr
* @author Sebastian Strasser
* @date 2006-07-10
*/
/***********************************************/

#ifndef __GROOPSGUI__TREEELEMENTGLOBAL__
#define __GROOPSGUI__TREEELEMENTGLOBAL__

#include "base/importGroops.h"
#include "tree/treeElement.h"
#include "tree/treeElementComplex.h"

/***** CLASS ***********************************/

class TreeElementGlobal : public TreeElementComplex
{
  Q_OBJECT

public:
  TreeElementGlobal(Tree *tree, TreeElementComplex *parentElement, XsdElementPtr xsdElement,
                    XmlNodePtr xmlNode);
 ~TreeElementGlobal();

  /** @brief Generate XML-tree.
  * recursively called for all children. */
  XmlNodePtr createXmlTree(bool /*createRootEvenIfEmpty*/) const override;

  // Inform about changes in variables
  // ---------------------------------
  /** @brief must be called if a variable is changed within this scope. */
  void updateParserResultsInScope() override;

  /** @brief must be called if a variable is added/removed/renamed within this scope.
   * recursively called for all children. */
  void updateLinksInScope() override;

  /** @brief inform this element about changed variables.
   * recursively called for all children. */
  void updateParserResults(VariableList &varList) override;

  /** @brief inform this element about changed variables.
   * recursively called for all children. */
  void updateLinks(QMap<QString, QString> &labelTypes) override;

  // Add/Remove elements
  // -------------------
  /** @brief Is it possible to insert an element with @a type before @a targetElement?
  * New elements can be added directly to the global section. */
  bool canAddChild(TreeElement *targetElement, const QString &type, const QString &label) const override;

  /** @brief add child before @a targetElement.
  * Is undoable. the new element is selected.
  * @return success */
  TreeElement *addChild(TreeElement *targetElement, const QString &type, const QString &label, XmlNodePtr xmlNode) override;
};

/***********************************************/

#endif /* __GROOPSGUI__TREEELEMENTGLOBAL__ */
