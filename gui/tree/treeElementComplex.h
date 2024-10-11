/***********************************************/
/**
* @file treeElementComplex.h
*
* @brief Abstract element with children (sequence or choice).
*
* @author Torsten Mayer-Guerr
* @author Sebastian Strasser
* @date 2006-10-07
*/
/***********************************************/

#ifndef __GROOPSGUI__TREEELEMENTCOMPLEX__
#define __GROOPSGUI__TREEELEMENTCOMPLEX__

#include <QJsonObject>
#include "base/importGroops.h"
#include "tree/treeElement.h"

/***** TYPES ***********************************/

class TreeElementAdd;

/***** CLASS ***********************************/

class TreeElementComplex : public TreeElement
{
  friend class TreeElementGlobal;

  QVector<QVector<TreeElement*>> children_;              // list of children for each choice
  QVector<XsdElementPtr>         xsdElementList;         // xsd for each choice
  QVector<QJsonObject>           overrideDefaultObjects; // defaultValue for each choice
  bool                           recieveAutoComments;

protected:
  QJsonObject defaultObject;

public:
  TreeElementComplex(Tree *tree, TreeElementComplex *parentElement, XsdElementPtr xsdElement, const QString &defaultOverride,
                     XmlNodePtr xmlNode, bool recieveAutoComments=false);
 ~TreeElementComplex();

protected:
  /** @brief add XML-nodes of the children to @a xmlNode. */
  void createXmlChildren(XmlNodePtr xmlNode, bool createRootEvenIfEmpty=false) const;

  // init the children lists
  // -----------------------
  /** @brief append a additional choice/sequence. */
  int addChoice(const QString &value, XsdElementPtr xsdComplex, const QJsonObject &overrideDefaultObject);

  /** @brief create elements for a choice with @a index. */
  void createChildrenElements(int index, XmlNodePtr xmlNode);

  /** @brief changes the current index.
  * calls TreeElement::selectIndex
  * creates the children elemens and items.
  * if an @a index is choosen without children all items are removed.
  * Set items expanded. */
  void setSelectedIndex(int index) override;

  // ========================================================

public:
  const QVector<TreeElement*> &children() const {static QVector<TreeElement*> empty; return isLinked() ? empty : children_[selectedIndex()];}

  // ========================================================

  // Inform about changes in variables
  // ---------------------------------
private:
  bool                   initializedVariables;
  VariableListPtr        varList;
  QMap<QString, QString> labelTypes;

public:
  /** @brief must be called if a variable is changed within this scope. */
  virtual void updateParserResultsInScope();

  /** @brief must be called if a variable is added/removed/renamed within this scope.
   * recursively called for all children. */
  virtual void updateLinksInScope();

  /** @brief inform this element about changed variables.
  * recursively called for all children. Returns the input @a varList. */
  VariableListPtr updateParserResults(VariableListPtr varList, Bool /*addVariableInReturn*/) override;

  /** @brief inform this element about changed variables.
  * recursively called for all children. */
  void updateLinks(QMap<QString, QString> &labelTypes) override;

  /** @brief must be called if a variable is renamed.
   * recursively called for all children.
   * @return next element can also be renamed. */
  bool renamedLink(const QString &oldLabel, const QString &newLabel) override;

  // ========================================================

  // Overwrite content
  // -----------------
private:
  class UndoCommandOverwriteChildren;
  friend class UndoCommandOverwriteChildren;

public:
  /** @brief copy the children of @a xmlNode into selected choice.
  * Is undoable. */
  void overwriteChildren(XmlNodePtr xmlNode);

  // ========================================================

private:
  TreeElement *skipCommentElements(TreeElement *targetElement, int index=-1) const;
  void childSetPushAutoComments(bool on, int index=-1);

  // Add/Remove elements
  // -------------------
protected:
  class UndoCommandRemoveAddChild : public TreeElement::UndoCommand
  {
    bool                toAdd;
    TreeElementComplex *parent;
    TreeElement        *targetElement;
    TreeElementAdd     *addElement;

  public:
    UndoCommandRemoveAddChild(TreeElement *newElement, TreeElement *targetElement, TreeElementComplex *parent, bool isAdd)
      : UndoCommand(newElement, (isAdd ? "add" : "remove")), toAdd(isAdd), parent(parent), targetElement(targetElement), addElement(newElement->elementAdd()) {}
   ~UndoCommandRemoveAddChild() {if(toAdd) delete treeElement;}

    void redo();
    void undo() {redo();}
  };
  friend class UndoCommandRemoveAddChild;

public:
  /** @brief Is it possible to insert an element with @a type before @a targetElement? */
  virtual bool canAddChild(TreeElement *targetElement, const QString &type, const QString &label) const;

  /** @brief Is it possible to remove this element from tree? */
  virtual bool canRemoveChild(TreeElement *element) const;

  /** @brief add child before/after @a targetElement.
  * Is undoable.
  * the new element is selected.
  * @return successfully created element */
  virtual TreeElement *addChild(TreeElement *targetElement, const QString &type, const QString &label, XmlNodePtr xmlNode);

  /** @brief remove child @a element.
  * Is undoable.
  * if successful the element is deleted.
  * @return success */
  virtual bool removeChild(TreeElement *element);

  // ========================================================

  // move elements
  // --------------
private:
  class UndoCommandMoveChild;
  friend class UndoCommandMoveChild;

public:
  /** @brief Is it possible to move an element before @a targetElement?
  * Only possible within the same unbounded list. */
  virtual bool canMoveChild(TreeElement *targetElement, TreeElement *element) const;

  /** @brief move child before @a targetElement.
  * Is undoable. */
  virtual void moveChild(TreeElement *targetElement, TreeElement *element);

  // ========================================================

  // TreeItem: visual element (one line in the tree widget)
  // ------------------------------------------------------
private:
  void createChildrenItems();
  void removeChildrenItems();

public:
  /** @brief create visual element.
  * Recursively for all children. */
  TreeItem *createItem(TreeItem *parent, TreeItem *after) override;

  /** @brief remove visual element.
  * Recursively for all children. */
  void removeItem() override;
};

/***********************************************/

#endif /* __GROOPSGUI__TREEELEMENTCOMPLEX__ */
