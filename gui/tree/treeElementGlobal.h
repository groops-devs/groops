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

  VariableList varList;
  XsdElementPtr findXsdElement(const QString &type) const;
  bool checkLabel(const QString &oldLabel, const QString &defaultLabel, QString &label) const;

public:
  TreeElementGlobal(Tree *tree, TreeElementComplex *parentElement, XsdElementPtr xsdElement,
                    XmlNodePtr xmlNode);
 ~TreeElementGlobal();

  /** @brief Generate XML-tree.
  * recursively called for all children. */
  XmlNodePtr createXmlTree(bool /*createRootEvenIfEmpty*/) const override;

  // special functions
  // -----------------
  /** @brief calls @a element->informAboutLink(elementInGlobal) for all elements in global. */
  void informAboutGlobalElements(TreeElement *element, bool recursively) const;

  /** @brief recomputes the variable list.
  * Must be called whenever a child is changed. */
  void updateVariableList();

  /** @brief GROOPS VariableList for parsing expressions. */
  const VariableList &variableList() const {return varList;}

   /** @brief Returns list of element names for all child elements in global. */
   QStringList names() const;

  /** @brief Returns list of global types.
  * Used by addGlobalDialog*/
  XsdComplexPtr xsdComplex() const {return xsdElement->complex;}

  // ========================================================

  /** @brief can @a element be set in the global section? */
  bool canSetGlobal(TreeElement *element) const;

  /** @brief Creates a new @a element in the global section.
  * Is undoable.
  * @return success? */
  bool setGlobal(TreeElement *element);

  // Called by the add button
  TreeElement *addNewChild(QString type=QString(), QString label=QString());

  // ========================================================

  // Add/Remove elements
  // -------------------
  /** @brief Is it possible to insert an element with @a type before @a targetElement?
  * New elements can be added directly to the global section. */
  bool canAddChild(TreeElement *targetElement, const QString &type) const override;

  /** @brief Is it possible to remove this element from tree?
  * All elements except the add element in the global section can be deleted. */
  bool canRemoveChild(TreeElement *element) const override;

  /** @brief add child before @a targetElement.
  * Is undoable. the new element is selected.
  * @return success */
  TreeElement *addChild(TreeElement *targetElement, const QString &type, XmlNodePtr xmlNode) override;

  /** @brief remove child @a element.
  * Is undoable. if successful the element is deleted.
  * @return success */
  bool removeChild(TreeElement *element) override;

  // Rename (label) global elements
  // ------------------------------
private:
  class UndoCommandRename;
  friend class UndoCommandRename;

public:
  /** @brief Is it possible to rename the element? */
  bool canRenameChild(TreeElement *element) const;

  /** @brief Rename the element.
  * Is undoable.
  * All elements informed by @a renamedLink() and @a updateParserResults(). */
  void renameChild(TreeElement *element, const QString &label);
};

/***********************************************/

#endif /* __GROOPSGUI__TREEELEMENTGLOBAL__ */
