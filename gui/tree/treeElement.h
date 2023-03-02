/***********************************************/
/**
* @file treeElement.h
*
* @brief Node of the tree.
*
* @author Torsten Mayer-Guerr
* @author Sebastian Strasser
* @date 2006-10-07
*/
/***********************************************/

#ifndef __GROOPSGUI__TREEELEMENT__
#define __GROOPSGUI__TREEELEMENT__

#include <QtDebug>
#include <QPointer>
#include <QList>
#include <QStringList>
#include <QComboBox>
#include <QUndoCommand>
#include "base/importGroops.h"
#include "base/xml.h"
#include "base/schema.h"

/***** TYPES ***********************************/

class  Tree;
class  TreeElement;
class  TreeElementComplex;
class  TreeElementProgram;
class  TreeElementAdd;
class  TreeItem;
class  QComboBox;

/***** CLASS ***********************************/

class TreeElement : public QObject
{
  Q_OBJECT

protected:
  class UndoCommand : public QUndoCommand
  {
  protected:
    TreeElement *treeElement;
    Tree        *tree;

  public:
    UndoCommand(TreeElement *treeElement, const QString &name);
    virtual ~UndoCommand() {}
  };

private:
  QString _name, _label, _defaultOverride, _originalName;

public:
TreeElement(Tree *tree, TreeElementComplex *parentElement, XsdElementPtr xsdElement, const QString &defaultOverride, XmlNodePtr xmlNode);
virtual ~TreeElement();

/** @brief Tree of TreeElements from XSD-Schema. *
* recursively called for all children. */
static TreeElement *newTreeElement(Tree *tree, TreeElementComplex *parentElement, XsdElementPtr xsdElement, const QString &defaultOverride, XmlNodePtr xmlNode, Bool fromFile);

  Tree               *tree;
  TreeElementComplex *parentElement;
  XsdElementPtr       xsdElement;

  QString   name()         const {return (_label.isEmpty()) ? _name : _label;}
  QString   originalName() const {return _originalName;} // original name in case element has been renamed
  QString   label()        const {return _label;}
  QString   loop()         const {return _loop;}
  QString   condition()    const {return _condition;}
  QString   type()         const {return (xsdElement) ? QString (xsdElement->type)       : QString();}
  QString   annotation()   const {return (xsdElement) ? QString (xsdElement->annotation) : QString();}
  QString   defaultValue() const {return !_defaultOverride.isEmpty() ? _defaultOverride : ((xsdElement) ? QString (xsdElement->defaultValue) : QString());}
  Bool      optional()     const {return (xsdElement) ? xsdElement->optional  : true;}
  Bool      unbounded()    const {return (xsdElement) ? xsdElement->unbounded : false;}
  Bool      disabled()     const {return _disabled;}
  Bool      isUnknown()    const {return type().isEmpty();}
  Bool      isProgram()    const {return type() == "programType" || type() == "programmeType";}

  void setLabel(const QString &label) {_name = type(); _label = label;}
  void setOriginalName(const QString &name) {_originalName = name;}

// ========================================================

public:
/** @brief Generate XML-tree.
* recursively called for all children. */
virtual XmlNodePtr getXML(Bool withEmptyNodes=false) const = 0;

protected:
/** @brief basis of a XML node (starting tag).
* sets the attributes ("label", "comment", "link"). */
XmlNodePtr getBaseXML() const;

// ========================================================

// Management of the content
// -------------------------
private:
  int         _selectedIndex;
  int         _valueCount;   // _selectedIndex >= _valueCount is a link
  QStringList _valueList;    // History or Choices

  class UndoCommandChangeSelectedIndex;
  friend class UndoCommandChangeSelectedIndex;

protected:
  virtual QString parseExpression(const QString &value) const;

public:
/** @brief Can values be edited?. */
virtual Bool isEditable() const {return false;}

/** @brief Is it a complex element?. */
virtual Bool isComplex() const {return false;}

/** @brief Has the element been renamed in the schema? */
virtual Bool isRenamed() const {return !originalName().isEmpty();}

/** @brief Is the selection at @p index unknown, i.e. does it not have a valid schema? */
virtual Bool isSelectionUnknown(int /*index*/) const {return false;}

/** @brief Is the selection at @p index renamed, i.e. has the selection been renamed in the schema? */
virtual Bool isSelectionRenamed(int /*index*/) const {return false;}

/** @brief the selected index. */
int selectedIndex() const {return _selectedIndex;}

/** @brief the selected value. */
const QString selectedValue() const {return ((_selectedIndex >= 0 && _selectedIndex < _valueList.size()) ? _valueList[_selectedIndex] : "");}

/** @brief list of values. */
const QStringList &valueList() const {return _valueList;}

/** @brief Number of values (excluding links). */
int valueCount() const {return _valueCount;}

/** @brief the selected result.
* The value of a link or the result of an expression. */
const QString selectedResult() const {return parseExpression(selectedValue());}

/** @brief is the selected value a link? */
Bool isLinked() const {return _selectedIndex >= _valueCount;}

/** @brief find index of a value.
* returns -1 if not found. */
int findValueIndex(const QString &value) const;

/** @brief find index of a link.
* returns -1 if not found. */
int findLinkIndex(const QString &value) const;

/** @brief changes the current index.
* - calls setSelectedIndex
* Is undoable. */
void changeSelectedIndex(int index);

/** @brief changes the current index.
* - calls changeSelectedIndex
* Is undoable. */
void changeSelectedValue(const QString &value);

/** @brief Rename a link in this element.
* recursively called for all children. */
virtual void renameLink(const QString &oldLabel, const QString &newLabel);

/** @brief In case this element is unknown, track or untrack it depending on @p track. */
virtual void trackUnknown(Bool track);

/** @brief In case this element is renamed, track or untrack it depending on @p track. */
virtual void trackRenamed(Bool track);

protected:
/** @brief changes the current index.
* if selected index is changed:
* - create auto comment
* - update item and comboBox
* - call tree->setChanged()
* - call newSelectedIndex() */
virtual void setSelectedIndex(int index);

/** @brief event handler of current index change.
* This event handler is called by setSelectedValue whenever the value is changed. */
virtual void newSelectedIndex(int /*index*/) {}

/** @brief insert a new value into the valueList.
* the selectedIndex/selectedValue is not changed.
* @return the index of the new value */
int insertNewValue(const QString &value, Bool prepend=false);

/** @brief collects all elements linked to label.
* recursively called for all children. */
virtual void getLinkedList(const QString &label, QList<TreeElement*> &list);

private:
/** @brief inform this element about a new link.
* recursively called for all children. */
virtual void newLink(TreeElement *elementInGlobal);

/** @brief inform this element about a removed link.
* recursively called for all children. */
virtual void removeLink(TreeElement *elementInGlobal);

/** @brief Update expression of associated TreeItem.
* recursively called for all children. */
virtual void updateExpression();

  friend class TreeElementComplex;
  friend class TreeElementGlobal;

// ========================================================

// comments
// --------
private:
  QString _comment;
  QString _autoComment;
  Bool    _pushComment;

  class UndoCommandSetComment;
  friend class UndoCommandSetComment;

  // @brief Receive a new auto-comment.
  // auto-comment is send by setSelectedValue() to the parent.
  void setAutoComment(const QString &text);

public:
/** @brief comment of this element. */
const QString &comment() const {return (_comment.isEmpty()) ? _autoComment : _comment;}

/** @brief sets a new comment.
* - item->setText() is called
* Is undoable. */
void setComment(const QString &text);

/** @brief Should this element send its value as auto-comment to the parent element?
* Auto-comment is send by setSelectedValue() to the parent.
* An Auto-comment never overwrites a real comment. */
void setPushAutoComments(Bool on=true);

// ========================================================

// Enable/Diable elements
// ----------------------
private:
  Bool _disabled;

  class UndoCommandSetEnabled;
  friend class UndoCommandSetEnabled;

public:
/** @brief Is it possible to disable the element?
* optional or unbounded elements can be disabled. */
Bool canDisabled() const;

/** @brief disable/enable the element.
* Is undoable. */
void setDisabled(Bool state=true);

// ========================================================

// Set/remove loop for elements
// ----------------------------
private:
  QString _loop;

  class UndoCommandSetLoop;
  friend class UndoCommandSetLoop;

public:
/** @brief Is it possible to set a loop for the element?
* Loops can only be set for unbounded elements, but not for the add element or global elements. */
Bool canSetLoop() const;

/** @brief Set/remove loop for element.
* Is undoable. */
void setLoop(const QString &loop);

/** @brief collects all elements using @p loop.
* recursively called for all children. */
virtual void getLoopList(const QString &loop, QList<TreeElement*> &list);

// ========================================================

// Set/remove conditions for elements
// ----------------------------
private:
  QString _condition;

  class UndoCommandSetCondition;
  friend class UndoCommandSetCondition;

public:
/** @brief Is it possible to set a condition for the element?
* Conditions can only be set for unbounded elements, but not for the add element or global elements. */
Bool canSetCondition() const;

/** @brief Set/remove condition for element.
* Is undoable. */
void setCondition(const QString &condition);

/** @brief collects all elements using @p condition.
* recursively called for all children. */
virtual void getConditionList(const QString &condition, QList<TreeElement*> &list);

// ========================================================

// Rename global elements
// ----------------------
private:
  class UndoCommandRename;
  friend class UndoCommandRename;

public:
/** @brief Is it possible to rename the element?
* Only global elements can be renamed. */
Bool canRename() const;

/** @brief Rename the element.
* Is undoable. */
void rename(const QString label);

// ========================================================

// Update name of elements in case of schema changes
// -------------------------------------------------
private:
  class UndoCommandUpdateName;
  friend class UndoCommandUpdateName;

public:
/** @brief Is it possible to update name of the element?
* Only elements with name changes in the schema can be updated. */
Bool canUpdateName() const;

/** @brief Update name of the element.
* Is undoable. */
void updateName();

// ========================================================

// Unbounded lists (managed in the add-element)
// --------------------------------------------
private:
  TreeElementAdd *_elementAdd;

public:
// Darf nur von TreeElementComplex verwendet werden
// Setzt entsprechend den Zaehler im addElement hoch oder runter
void setElementAdd(TreeElementAdd *elementAdd);

/** @brief If this is unbounded() returns the add element.
* the add element contains the counter of elements in the unbounded list.
* Overwritten by ElementAdd. */
virtual TreeElementAdd *elementAdd() const {return _elementAdd;}

/** @brief Is this the add-element?.
* Overwritten by ElementAdd. */
virtual Bool isElementAdd() const {return false;}

// ========================================================

public:
/** @brief collects all program in a list.
* recursively called for all children. */
virtual void getProgramList(QList<TreeElement*> &) {}

// ========================================================

// TreeItem: visual element (one line in the tree widget)
// ------------------------------------------------------
private:
  TreeItem *_item;

public:
/** @brief pointer to visual element. */
TreeItem *item() const {return _item;}

/** @brief create visual element. */
virtual TreeItem *createItem(TreeItem *parent, TreeItem *after);

/** @brief remove visual element. */
virtual void removeItem();

// ========================================================

private:
  QPointer<QComboBox> comboBox;

  class UndoCommandEdit;
  friend class UndoCommandEdit;

protected:
/** @brief create a combo box.
* The combo box is filled with values and links.
* Connections are set. */
QComboBox *createComboBox(Bool isEditable);

public:
/** @brief Set tooltip for combo box if text does not fit completely. */
void comboBoxSetToolTip();

/** @brief Create the editor element.
* Is called from item when the it gets the focus.
* The editor is deleted by the caller. */
virtual QWidget *createEditor() {return nullptr;}

/** @brief Interacts with the element (open dialog, ...).
 * Used by TreeElementFileName, TreeElementProgram, and TreeElementTime. */
virtual void interact() {}

/** @brief event handler called by item when it gets selected.
* Used by TreeElementChoice. */
virtual void startSelected()  {}

/** @brief event handler called by item when it losts the selection.
* Used by TreeElementChoice. */
virtual void stopSelected() {}

public slots:
  void comboBoxActivated(int index);
  void comboBoxTextEdited(const QString &text);
  void columnResized(int column, int oldSize, int newSize);
};

/***********************************************/

#endif
