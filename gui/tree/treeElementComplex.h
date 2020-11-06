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

/***** CLASS ***********************************/

class TreeElementComplex : public TreeElement
{
  friend class TreeElementGlobal;

  QVector<XsdElementPtr>           xsdElementList;         // xsd for each choice
  QVector<QJsonObject>             overrideDefaultObjects; // defaultValue for each choice
  QVector<QVector<TreeElement*>>   children;               // list of children for each choice
  Bool                             recieveAutoComments;
  Bool                            _isElementGlobal;

  Bool hasChildren(int index) const;

protected:
  QJsonObject defaultObject;

public:
  TreeElementComplex(Tree *tree, TreeElementComplex *parentElement, XsdElementPtr xsdElement, const QString &defaultOverride,
                     XmlNodePtr xmlNode, Bool recieveAutoComments=false, Bool isElementGlobal=false);
  virtual ~TreeElementComplex() override;

/** @brief is this the global element. */
Bool isElementGlobal() const {return _isElementGlobal;}

/** @brief Is it a complex element?. */
virtual Bool isComplex() const override {return true;}

/** @brief Is the selection at @p index unknown, i.e. does it not have a valid schema? */
virtual Bool isSelectionUnknown(int index) const override;

protected:
/** @brief add XML-nodes of the children to @a xmlNode. */
void getChildrenXML(XmlNodePtr xmlNode, bool withEmptyNodes = false) const;

// init the children lists
// -----------------------
/** @brief append a additional choice/sequence. */
int addChoice(const QString &value, XsdElementPtr xsdComplex, const QJsonObject &overrideDefaultObject);

/** @brief create elements for a choice with @a index. */
void createChildrenElements(int index, XmlNodePtr xmlNode=XmlNodePtr(nullptr));

protected:
/** @brief This event handler is called by setSelectedValue whenever the value is changed.
* creates the children elemens and items.
* if an @a index is choosen without children all items are removed.
* Set items expanded. */
virtual void newSelectedIndex(int index) override;

// ========================================================

// Management of the content
// -------------------------
public:
/** @brief inform this element about a new/changed link.
* recursively called for all children. */
virtual void newLink(TreeElement *elementInGlobal) override;

/** @brief inform this element about a removed link.
* recursively called for all children. */
virtual void removeLink(TreeElement *elementInGlobal) override;

/** @brief Rename a link in this element.
* recursively called for all children. */
virtual void renameLink(const QString &oldLabel, const QString &newLabel) override;

/** @brief In case this element is unknown, track or untrack it depending on @p track.
* recursively called for all children. */
virtual void trackUnknown(bool track) override;

/** @brief collects all elements linked to label.
* recursively called for all children. */
virtual void getLinkedList(const QString &label, QList<TreeElement*> &list) override;

/** @brief collects all elements using @p loop.
* recursively called for all children. */
virtual void getLoopList(const QString &label, QList<TreeElement*> &list) override;

/** @brief collects all elements using @p condition.
* recursively called for all children. */
virtual void getConditionList(const QString &label, QList<TreeElement*> &list) override;

/** @brief Update expression of associated TreeItem.
* recursively called for all children. */
virtual void updateExpression() override;

// ========================================================

// Overwrite content
// -----------------
private:
  class UndoCommandOverwrite;
  friend class UndoCommandOverwrite;

public:
/** @brief Copy the content of @a element into this.
* Is undoable.
* @return success */
Bool overwrite(TreeElement *element, const QString &type, XmlNodePtr xmlNode);

// ========================================================

// Add/Remove elements
// -------------------
private:
  class UndoCommandRemoveAddChild;
  friend class UndoCommandRemoveAddChild;

public:
/** @brief Is it possible to insert an element with @a type before @a beforeElement? */
virtual Bool canAddChild(TreeElement *beforeElement, const QString &type) const;

/** @brief Is it possible to remove this element from tree? */
virtual Bool canRemoveChild(TreeElement *element) const;

/** @brief add child before @a beforeElement.
* Is undoable.
* the new element is selected.
* @return success */
virtual Bool addChild(TreeElement *beforeElement, const QString &type, XmlNodePtr xmlNode, bool moved = false);

/** @brief remove child @a element.
* Is undoable.
* if successful the element is deleted.
* @return success */
virtual Bool removeChild(TreeElement *element);

// ========================================================

// Move elements
// -------------
private:
  class UndoCommandMoveChild;
  friend class UndoCommandMoveChild;

public:
/** @brief Is it possible to move this element within the tree? */
virtual Bool canMoveChild(TreeElement *element) const;

/** @brief Move child @a element.
* Is undoable.
* @return success */
virtual Bool moveChild(TreeElement *element, TreeElement *beforeElement);

// ========================================================

public:
/** @brief Collects all programs in a list.
* recursively called for all children. */
virtual void getProgramList(QList<TreeElement*> &list) override;

/** @brief Return the number of children elements. */
virtual int childrenCount() const;

/** @brief Return the child element at @p index. */
virtual TreeElement* childAt(int index) const;

// ========================================================

// TreeItem: visual element (one line in the tree widget)
// ------------------------------------------------------
private:
  void createChildrenItems();
  void removeChildrenItems();

public:
/** @brief create visual element.
* Recursively for all children. */
virtual TreeItem *createItem(TreeItem *parent, TreeItem *after) override;

/** @brief remove visual element.
* Recursively for all children. */
virtual void removeItem() override;

};

/***********************************************/

class TreeElementComplex::UndoCommandRemoveAddChild : public TreeElement::UndoCommand
{
  Bool                isAdd;
  Bool                isRemoved;
  TreeElementComplex *parent;
  TreeElement        *beforeElement;
  TreeElement        *newElement;
  QList<TreeElement*> elementsWithRemovedLink;

public:
  UndoCommandRemoveAddChild(TreeElement *_newElement, TreeElement *_beforeElement, TreeElementComplex *_parent, Bool _isAdd);
 ~UndoCommandRemoveAddChild();

  void addChild();
  void removeChild();

  void redo() {if(isAdd)  addChild(); else removeChild();}
  void undo() {if(!isAdd) addChild(); else removeChild();}
};

/***********************************************/

class TreeElementComplex::UndoCommandMoveChild : public TreeElement::UndoCommand
{
  TreeElement        *element;
  TreeElement        *beforeElement;
  TreeElementComplex *parent;

public:
  UndoCommandMoveChild(TreeElement *_element, TreeElement *_beforeElement, TreeElementComplex *_parent)
    : UndoCommand(_element, "move"),
      element(_element), beforeElement(_beforeElement), parent(_parent) {}

  void redo();
  void undo() {redo();}
};

/***********************************************/

#endif /* __GROOPSGUI__TREEELEMENTCOMPLEX__ */
