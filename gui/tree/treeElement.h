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
class  TreeElementLoopCondition;
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
  QString _name, _schemaName, _label;

public:
  TreeElement(Tree *tree, TreeElementComplex *parentElement, XsdElementPtr xsdElement, const QString &defaultOverride, XmlNodePtr xmlNode);
  virtual ~TreeElement();

  /** @brief Tree of TreeElements from XSD-Schema. *
  * recursively called for all children. */
  static TreeElement *newTreeElement(Tree *tree, TreeElementComplex *parentElement, XsdElementPtr xsdElement,
                                     const QString &defaultOverride, XmlNodePtr xmlNode, bool fillWithDefaults);

  Tree               *tree;
  TreeElementComplex *parentElement;
  XsdElementPtr       xsdElement;
  QString             defaultValue;

  virtual QString name()              const {return (_label.isEmpty()) ? _name : _label;}
  virtual QString xmlName()           const {return _name;}
  virtual QString label()             const {return _label;}
  virtual QString annotation()        const {return (xsdElement) ? QString (xsdElement->annotation) : QString();}
  virtual QString type()              const {return (xsdElement) ? QString (xsdElement->type)       : QString();}
  virtual bool    optional()          const {return (xsdElement) ? xsdElement->optional  : true;}
  virtual bool    unbounded()         const {return (xsdElement) ? xsdElement->unbounded : false;}
  virtual bool    isRenamedInSchema() const {return (xmlName() != _schemaName);}

  /** @brief Is the selection without coresponding schema? (choice only) */
  virtual bool isSelectionUnknown(int /*index*/) const  {return false;}

  /** @brief Is the selection renamed (new name in the schema)? (choice only)  */
  virtual bool isSelectionRenamedInSchema(int /*index*/) const {return false;}

  // ========================================================

public:
  /** @brief Generate XML-tree.
  * recursively called for all children. */
  virtual XmlNodePtr createXmlTree(bool createRootEvenIfEmpty=false) const = 0;

protected:
  /** @brief basis of a XML node (starting tag).
  * sets the attributes ("label", "comment", "link", ...). */
  XmlNodePtr createXmlBaseNode() const;

  // ========================================================

  // Management of the content
  // -------------------------
protected:
  int         _selectedIndex;
  int         _valueCount;        // _selectedIndex >= _valueCount is a link
  int         _brokenLinkIndex;
  QStringList _valueList;         // History or Choices

  class UndoCommandChangeSelectedIndex;
  friend class UndoCommandChangeSelectedIndex;

public:
  /** @brief Can values be edited?. */
  virtual bool isEditable() const {return false;}

  /** @brief the selected index. */
  int selectedIndex() const {return _selectedIndex;}

  /** @brief the selected value. */
  const QString &selectedValue() const {return _valueList[_selectedIndex];}

  /** @brief the selected value as result of parser. */
  virtual QString selectedResult() const {return QString();}

  /** @brief is the selected value a link? */
  bool isLinked() const {return _selectedIndex >= _valueCount;}

  /** @brief is the selected value a broken link? */
  bool isBrokenLinked() const {return _selectedIndex == _brokenLinkIndex;}

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

  /** @brief changes the current index.
  * - calls changeSelectedIndex
  * Is undoable. */
  void changeSelectedLink(const QString &link);

protected:
  /** @brief changes the current index.
  * if selected index is changed:
  * - create auto comment
  * - update item and comboBox
  * - call tree->setChanged() */
  virtual void setSelectedIndex(int index);

  /** @brief insert a new value into the valueList.
  * the selectedIndex/selectedValue is not changed.
  * @return the index of the new value */
  int insertNewValue(const QString &value, bool prepend=false);

  // Inform about changes in variables
  // ---------------------------------
public:
  /** @brief inform this element about a renamed link.
   * recursively called for all children.
   * @return next element can also be renamed. */
  virtual bool renamedLink(const QString &oldLabel, const QString &newLabel);

protected:
  /** @brief inform this element about changed variables.
  * recursively called for all children. */
  virtual void updateParserResults(VariableList &varList);

  /** @brief inform this element about changed variables.
  * recursively called for all children. */
  virtual void updateLinks(QMap<QString, QString> &labelTypes);

  friend class TreeElementComplex;
  friend class TreeElementGlobal;
  friend class TreeElementLoopCondition;

  // ========================================================

  // Overwrite
  // ---------
protected:
  /** @brief Copy loop, conditions, comments and links from @a xmlNode.
  * Is undoable.
  * if @a contentOnly loop,conditions and comments are not copied.
  * @return is not link -> do more */
  bool baseOverwrite(XmlNodePtr xmlNode, bool contentOnly=false);

public:
  /** @brief Is it possible to overweite the element? */
  virtual bool canOverwrite(const QString &/*type*/) {return false;}

  /** @brief Copy the content of @a xmlNode into this.
  * Is undoable.
  * @return success */
  virtual bool overwrite(const QString &/*type*/, XmlNodePtr /*xmlNode*/, bool /*contentOnly*/=false) {return false;}

  // ========================================================

  // Set/remove loop for elements
  // ----------------------------
private:
  class UndoCommandSetLoop;
  friend class UndoCommandSetLoop;

public:
  TreeElementLoopCondition *loop;

  /** @brief Is it possible to set a loop for the element?
  * Loops can only be set for unbounded elements, but not for the add element or global elements. */
  virtual bool canSetLoop() const;

  /** @brief Set loop for element.
  * Is undoable. */
  void setLoop();

  // ========================================================

  // Set/remove conditions for elements
  // ----------------------------------
private:
  class UndoCommandSetCondition;
  friend class UndoCommandSetCondition;

public:
  TreeElementLoopCondition *condition;

  /** @brief Is it possible to set a condition for the element?
  * Conditions can only be set for unbounded elements, but not for the add element or global elements. */
  virtual bool canSetCondition() const;

  /** @brief Set condition for element.
  * Is undoable. */
  void setCondition();

  /** @brief Remove loop/condition for element.
  * Is undoable. */
  void removeLoopOrCondition(TreeElement *element);

  // ========================================================

  // Enable/Diable elements
  // ----------------------
private:
  bool _disabled;

  class UndoCommandSetEnabled;
  friend class UndoCommandSetEnabled;

public:
  /** @brief Is the element disabled? */
  bool disabled() const {return _disabled;}

  /** @brief Is it possible to disable the element?
   * optional or unbounded elements can be disabled. */
  virtual bool canDisabled() const;

  /** @brief disable/enable the element.
   * Is undoable. */
  void setDisabled(bool disabled=true);

  // ========================================================

  // Rename elements
  // ---------------
private:
  class UndoCommandRename;
  friend class UndoCommandRename;

public:
  /** @brief Is it possible to rename the element? */
  virtual bool canRename() const;

  /** @brief Rename the element.
  * Is undoable.
  * All elements informed by @a renamedLink() and @a updateParserResults(). */
  virtual bool rename(const QString &label);

  // ========================================================

  // Update name of elements in case of schema changes
  // -------------------------------------------------
protected:
  class UndoCommandUpdateName : public UndoCommand
  {
    QString name, value;

  public:
    UndoCommandUpdateName(TreeElement *treeElement, const QString &value)
    : UndoCommand(treeElement, "update name"), name(treeElement->_schemaName), value(value) {}

    void redo();
    void undo() {redo();}
  };
  friend class UndoCommandUpdateName;

public:
  /** @brief Is it possible to update name of the element?
  * Only elements with name changes in the schema can be updated. */
  virtual bool canUpdateName() const;

  /** @brief Update name of the element.
  * Is undoable. */
  virtual void updateName();

  // ========================================================

  // comments
  // --------
private:
  QString _comment;
  QString _autoComment;
  bool    _pushComment;

  class UndoCommandSetComment;
  friend class UndoCommandSetComment;

  // @brief Receive a new auto-comment.
  // auto-comment is send by setSelectedValue() to the parent.
  void setAutoComment(const QString &text);

public:
  /** @brief comment of this element. */
  const QString &comment() const {return (_comment.isEmpty()) ? _autoComment : _comment;}

  /** @brief Can element be commented. */
  virtual bool canComment() const {return true;}

  /** @brief sets a new comment.
  * - item->setText() is called
  * Is undoable. */
  void setComment(const QString &text);

  /** @brief Should this element send its value as auto-comment to the parent element?
  * Auto-comment is send by setSelectedValue() to the parent.
  * An Auto-comment never overwrites a real comment. */
  void setPushAutoComments(bool on=true);

  // ========================================================

  // Unbounded lists (managed in the add-element)
  // --------------------------------------------
private:
  TreeElementAdd *_elementAdd;

public:
  /** @brief Add this element to an unbounded list.
  * Must only be used in TreeElementComplex.
  * Adjust the element counter in addElement. */
  void setElementAdd(TreeElementAdd *elementAdd);

  /** @brief If this is unbounded() returns the add element.
  * the add element contains the counter of elements in the unbounded list.
  * Overwritten by ElementAdd. */
  virtual TreeElementAdd *elementAdd() const {return _elementAdd;}

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

private:
  QPointer<QComboBox> comboBox;
  void comboBoxSetToolTip();

protected:
  class UndoCommandEdit : public UndoCommand
  {
    bool    isCreated, isLink;
    QString newValue, oldValue;

  public:
    UndoCommandEdit(TreeElement *treeElement, const QString &text)
    : UndoCommand(treeElement, "edit"), isCreated(treeElement->findValueIndex(text)<0), isLink(treeElement->isLinked()),
      newValue(text), oldValue(treeElement->selectedValue())
      {setText("edit "+treeElement->name()+": "+text);}

    int  id() const {return 999;}
    void redo();
    void undo();
    bool mergeWith(const QUndoCommand *command);
  };
  friend class UndoCommandEdit;

  /** @brief create a combo box.
  * The combo box is filled with values and links.
  * Connections are set. */
  QComboBox *createComboBox(bool isEditable);

public:
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
