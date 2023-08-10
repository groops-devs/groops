/***********************************************/
/**
* @file treeElement.cpp
*
* @brief Node of the tree.
*
* @author Torsten Mayer-Guerr
* @author Sebastian Strasser
* @date 2006-10-07
*/
/***********************************************/

#include <QtDebug>
#include <QComboBox>
#include <QLineEdit>
#include <QHeaderView>
#include "base/importGroops.h"
#include "base/xml.h"
#include "base/schema.h"
#include "tree/tree.h"
#include "tree/treeItem.h"
#include "tree/treeElementComplex.h"
#include "tree/treeElementSequence.h"
#include "tree/treeElementChoice.h"
#include "tree/treeElementSimple.h"
#include "tree/treeElementBool.h"
#include "tree/treeElementFileName.h"
#include "tree/treeElementTime.h"
#include "tree/treeElementGlobal.h"
#include "tree/treeElementProgram.h"
#include "tree/treeElementAdd.h"
#include "tree/treeElement.h"

/***********************************************/

TreeElement::UndoCommand::UndoCommand(TreeElement *_treeElement, const QString &name)
  : QUndoCommand(name+" "+_treeElement->name()+" "+_treeElement->selectedValue()),
    treeElement(_treeElement),
    tree(_treeElement->tree) {}

/***********************************************/
/***********************************************/

TreeElement *TreeElement::newTreeElement(Tree *tree, TreeElementComplex *parentElement, XsdElementPtr xsdElement, const QString &defaultOverride,
                                         XmlNodePtr xmlNode, bool fillWithDefaults)
{
  try
  {
    if(xsdElement->complex)
    {
      if((xsdElement->complex->type=="sequence") && xsdElement->names.contains("global"))
        return new TreeElementGlobal(tree, parentElement, xsdElement, xmlNode);
      if((xsdElement->type=="programType") || (xsdElement->type== "programmeType"))
        return new TreeElementProgram(tree, parentElement, xsdElement, defaultOverride, xmlNode, fillWithDefaults);
      if(xsdElement->complex->type=="sequence")
        return new TreeElementSequence(tree, parentElement, xsdElement, defaultOverride, xmlNode, fillWithDefaults);
      if(xsdElement->complex->type=="choice")
        return new TreeElementChoice(tree, parentElement, xsdElement, defaultOverride, xmlNode, fillWithDefaults);
      throw(Exception("unknown complex element"));
    }

    if(xsdElement->type=="boolean")
      return new TreeElementBool(tree, parentElement, xsdElement, defaultOverride, xmlNode, fillWithDefaults);
    if(xsdElement->type=="filename")
      return new TreeElementFileName(tree, parentElement, xsdElement, defaultOverride, xmlNode, fillWithDefaults);
    if(xsdElement->type=="time")
      return new TreeElementTime(tree, parentElement, xsdElement, defaultOverride, xmlNode, fillWithDefaults);
    return new TreeElementSimple(tree, parentElement, xsdElement, defaultOverride, xmlNode, fillWithDefaults);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/
/***********************************************/

TreeElement::TreeElement(Tree *tree, TreeElementComplex *parentElement, XsdElementPtr xsdElement, const QString &defaultOverride, XmlNodePtr xmlNode)
{
  try
  {
    this->tree             = tree;
    this->parentElement    = parentElement;
    this->xsdElement       = xsdElement;
    if(xsdElement)
      this->defaultValue  = QString (xsdElement->defaultValue);
    if(!defaultOverride.isEmpty())
      this->defaultValue  = defaultOverride;
    this->_selectedIndex   = -1;
    this->_valueCount      = 0;
    this->_pushComment     = false;
    this->_disabled        = false;
    this->_elementAdd      = nullptr;
    this->_item            = nullptr;
    this->comboBox         = nullptr;

    if(xsdElement)
      _name = _schemaName = xsdElement->names.front();

    // evaluation of the XML nodes
    // ---------------------------
    if(xmlNode)
    {
      _name = xmlNode->getName();

      readAttribute(xmlNode, "label",     _label);
      readAttribute(xmlNode, "loop",      _loop);
      readAttribute(xmlNode, "condition", _condition);
      readAttribute(xmlNode, "comment",   _comment);

      XmlAttrPtr attr = xmlNode->getAttribute("disabled");
      if(attr && (attr->text == "1"))
        _disabled = true;

      // is this a link?
      attr = xmlNode->getAttribute("link");
      if(attr)
      {
        _valueList.insert(_valueCount, attr->text);
        _selectedIndex = _valueCount;
      }
    }

    if(_name.isEmpty())
      _name = "comment";
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

TreeElement::~TreeElement()
{
  removeItem();
  setElementAdd(nullptr);
}

/***********************************************/

XmlNodePtr TreeElement::createXmlBaseNode() const
{
  try
  {
    XmlNodePtr xmlNode(new XmlNode(_name));
    if(!label().isEmpty())     writeAttribute(xmlNode, "label",     label());
    if(!loop().isEmpty())      writeAttribute(xmlNode, "loop",      loop());
    if(!condition().isEmpty()) writeAttribute(xmlNode, "condition", condition());
    if(!comment().isEmpty())   writeAttribute(xmlNode, "comment",   comment());
    if(isLinked())             writeAttribute(xmlNode, "link",      selectedValue());
    if(disabled())             writeAttribute(xmlNode, "disabled",  "1");
    return xmlNode;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/
/***********************************************/

int TreeElement::findValueIndex(const QString &value) const
{
  try
  {
    for(int i=0; i<_valueCount; i++)
      if(value == _valueList[i])
        return i;
    return -1;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

int TreeElement::findLinkIndex(const QString &value) const
{
  try
  {
    for(int i=_valueCount; i<_valueList.size(); i++)
      if(value == _valueList[i])
        return i;
    return -1;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

class TreeElement::UndoCommandChangeSelectedIndex : public TreeElement::UndoCommand
{
  int index;
public:
  UndoCommandChangeSelectedIndex(TreeElement *treeElement, int index) : UndoCommand(treeElement, "change"), index(index)
      {setText("change "+treeElement->name()+" from "+treeElement->selectedValue()+" to "+treeElement->_valueList[index]);}
  void redo();
  void undo() {redo();}
};

/***********************************************/

void TreeElement::UndoCommandChangeSelectedIndex::redo()
{
  try
  {
    if(treeElement->item())
      tree->setSelectedItem(treeElement->item());

    int indexOld = treeElement->selectedIndex();
    treeElement->setSelectedIndex(index);
    index = indexOld;
    tree->treeChanged();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

void TreeElement::changeSelectedIndex(int index)
{
  try
  {
    if((index<0) || (index>=_valueList.size()))
      throw(Exception("index out of range"));
    if(index == _selectedIndex)
      return;
    tree->undoStack->push(new UndoCommandChangeSelectedIndex(this, index));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

void TreeElement::changeSelectedValue(const QString &value)
{
  changeSelectedIndex(insertNewValue(value, true/*prepend*/));
}

/***********************************************/
/***********************************************/

void TreeElement::setSelectedIndex(int index)
{
  try
  {
    if((index < 0) || (index >= _valueList.size()))
      throw(Exception("index out of range"));

    if(index == _selectedIndex)
      return;
    _selectedIndex = index;

    // change of global variable?
    // not called during first initialization (rootElement is not set yet)
    if(dynamic_cast<TreeElementGlobal*>(parentElement) && tree->rootElement)
      dynamic_cast<TreeElementGlobal*>(parentElement)->updateVariableList();

    // send auto comment
    if(_pushComment && parentElement)
      parentElement->setAutoComment(selectedValue());

    if(item())
      item()->updateValue();

    if(comboBox)
    {
      comboBox->setCurrentIndex(_selectedIndex);
      comboBoxSetToolTip();
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

int TreeElement::insertNewValue(const QString &value, bool prepend)
{
  try
  {
    // is the value already in the list?
    int index = findValueIndex(value);
    if(index >= 0)
      return index;

    if(value.isEmpty())
      prepend = true;

    // insert at index;
    index = (prepend) ? 0 : _valueCount;
    if((index == 0) && (_valueCount > 0) && _valueList[0].isEmpty())
      index++;
    _valueList.insert(index, value);
    _valueCount++;
    if(_selectedIndex >= index)
      _selectedIndex++;

    // update comboBox
    if(comboBox)
      comboBox->insertItem(index, value);

    // empty values only allowed at the beginning
    for(int i=1; i<_valueList.size(); i++)
      if(_valueList[i].isEmpty())
      {
        if(_selectedIndex >= i)
          _selectedIndex--;
        _valueList.removeAt(i);
        _valueCount--;
        if(comboBox)
          comboBox->removeItem(i);
        i--;
      }

    if(comboBox && (comboBox->count() != _valueList.size()))
      throw(Exception("comboBox not consistent"));

    return index;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/
/***********************************************/

class TreeElement::UndoCommandAddRemoveLink : public TreeElement::UndoCommand
{
  bool    isAdd;
  QString label;

public:
  UndoCommandAddRemoveLink(TreeElement *element, const QString &label, bool isAdd)
  : UndoCommand(element, (isAdd ? "add" : "remove")), isAdd(isAdd), label(label)  {}
  ~UndoCommandAddRemoveLink() {}

  void redo();
  void undo() {redo();}
};

/***********************************************/

void TreeElement::UndoCommandAddRemoveLink::redo()
{
  if(isAdd)
  {
    isAdd = false;
    treeElement->_valueList.push_back(label);
    if(treeElement->comboBox)
      treeElement->comboBox->insertItem(treeElement->comboBox->count(), QIcon(":/icons/scalable/link.svg"), label);
  }
  else
  {
    isAdd = true;
    int index = treeElement->findLinkIndex(label);
    if(index < treeElement->_selectedIndex)
      treeElement->_selectedIndex--;
    treeElement->_valueList.removeAt(index);
    if(treeElement->comboBox)
      treeElement->comboBox->removeItem(index);
  }
}

/***********************************************/
/***********************************************/

void TreeElement::informAboutLink(TreeElement *element, bool /*recursively*/)
{
  if((element != this) && (element->type() == type()) && (findLinkIndex(element->label()) < 0))
    _valueList.push_back(element->label());
}

/***********************************************/

void TreeElement::addedLink(TreeElement *element)
{
  try
  {
    // is this link already in the list?
    if((element != this) && (element->type() == type()) && (findLinkIndex(element->label()) < 0))
      tree->undoStack->push(new UndoCommandAddRemoveLink(this, element->label(), true/*add*/));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

void TreeElement::removedLink(TreeElement *element)
{
  try
  {
    if(element==this)
      return;

    // remove loop
    if(element->label() == _loop)
      setLoop("");

    // remove condition
    if(element->label() == _condition)
      setCondition("");

    // is this link in the list?
    int index = findLinkIndex(element->label());
    if(index < 0)
      return;

    if(index == _selectedIndex)
      overwrite(element->type(), element->createXmlTree(), true/*contentOnly*/);

    tree->undoStack->push(new UndoCommandAddRemoveLink(this, element->label(), false/*add*/));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

void TreeElement::renamedLink(const QString &oldLabel, const QString &newLabel)
{
  try
  {
    // is this link set as loop?
    if(loop() == oldLabel)
    {
      _loop = newLabel;
      if(item())
        item()->updateName();
    }

    // is this link set as condition?
    if(condition() == oldLabel)
    {
      _condition = newLabel;
      if(item())
        item()->updateName();
    }

    // is this link in the list?
    int index = findLinkIndex(oldLabel);
    if(index<0)
      return;

    // replace label
    _valueList.replace(index, newLabel);
    if(item())
      item()->updateValue();
    if(comboBox && index < comboBox->count())
      comboBox->setItemText(index, newLabel);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/
/***********************************************/

class TreeElement::UndoCommandSetComment : public TreeElement::UndoCommand
{
  QString text;

public:
  UndoCommandSetComment(TreeElement *treeElement, const QString &_text)
    : UndoCommand(treeElement, "comment"),
      text(_text) {}

  void redo();
  void undo() {redo();}
};

/***********************************************/

void TreeElement::UndoCommandSetComment::redo()
{
  try
  {
    if(treeElement->item())
      tree->setSelectedItem(treeElement->item());
    std::swap(treeElement->_comment, text);
    if(treeElement->item())
      treeElement->item()->updateComment();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

void TreeElement::setComment(const QString &text)
{
  try
  {
    if(_comment != text) // is something changed?
      tree->undoStack->push(new UndoCommandSetComment(this, text));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/
/***********************************************/

void TreeElement::setAutoComment(const QString &text)
{
  try
  {
    if((_comment == _autoComment) || (_comment == text))
      _comment = QString();
    _autoComment = text;
    if(item())
      item()->updateComment();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

void TreeElement::setPushAutoComments(bool on)
{
  try
  {
    _pushComment = on;
    if(_pushComment && parentElement)
      parentElement->setAutoComment(selectedValue());
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/
/***********************************************/

bool TreeElement::baseOverwrite(XmlNodePtr xmlNode, bool contentOnly)
{
  if(!contentOnly)
  {
    QString loop, condition, comment;
    readAttribute(xmlNode, "loop",      loop);
    readAttribute(xmlNode, "condition", condition);
    readAttribute(xmlNode, "comment",   comment);
    setLoop(loop);
    setCondition(condition);
    setComment(comment);
  }

  QString link;
  readAttribute(xmlNode, "link", link);
  if(!link.isEmpty())
    changeSelectedIndex(findLinkIndex(link));
  return link.isEmpty();
}

/***********************************************/
/***********************************************/

class TreeElement::UndoCommandSetLoop : public TreeElement::UndoCommand
{
  QString loop;

public:
  UndoCommandSetLoop(TreeElement *treeElement, QString loop)
    : UndoCommand(treeElement, (loop.isEmpty() ? "remove loop "+treeElement->loop()+" from" : "set loop "+loop+" for")),
      loop(loop) {}

  void redo();
  void undo() {redo();}
};

/***********************************************/

void TreeElement::UndoCommandSetLoop::redo()
{
  try
  {
    if(treeElement->item())
      tree->setSelectedItem(treeElement->item());

    QString oldLoop = treeElement->loop();
    treeElement->_loop = loop;
    loop = oldLoop;
    if(treeElement->item())
      treeElement->item()->updateName();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

bool TreeElement::canSetLoop() const
{
  if(dynamic_cast<TreeElementGlobal*>(parentElement))
    return false;
  return unbounded();
}

/***********************************************/

void TreeElement::setLoop(const QString &loop)
{
  try
  {
    if((loop != _loop) && (canSetLoop() || loop.isEmpty()))
      tree->undoStack->push(new UndoCommandSetLoop(this, loop));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/
/***********************************************/

class TreeElement::UndoCommandSetCondition : public TreeElement::UndoCommand
{
  QString condition;

public:
  UndoCommandSetCondition(TreeElement *treeElement, QString condition)
    : UndoCommand(treeElement, (condition.isEmpty() ? "remove condition "+treeElement->condition()+" from" : "set condition "+condition+" for")),
      condition(condition) {}

  void redo();
  void undo() {redo();}
};

/***********************************************/

void TreeElement::UndoCommandSetCondition::redo()
{
  try
  {
    if(treeElement->item())
      tree->setSelectedItem(treeElement->item());

    QString oldCondition = treeElement->condition();
    treeElement->_condition = condition;
    condition = oldCondition;
    if(treeElement->item())
      treeElement->item()->updateName();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

bool TreeElement::canSetCondition() const
{
  if(dynamic_cast<TreeElementGlobal*>(parentElement))
    return false;
  return optional() || unbounded();
}

/***********************************************/

void TreeElement::setCondition(const QString &condition)
{
  try
  {
    if((condition != _condition) && (canSetCondition() || condition.isEmpty()))
      tree->undoStack->push(new UndoCommandSetCondition(this, condition));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/
/***********************************************/

class TreeElement::UndoCommandSetEnabled : public TreeElement::UndoCommand
{
  bool disabled;

public:
  UndoCommandSetEnabled(TreeElement *treeElement, bool disabled)
  : UndoCommand(treeElement, (disabled ? "disable" : "enable")),
  disabled(disabled) {}

  void redo();
  void undo() {redo();}
};

/***********************************************/

void TreeElement::UndoCommandSetEnabled::redo()
{
  try
  {
    if(treeElement->item())
      tree->setSelectedItem(treeElement->item());

    std::swap(disabled, treeElement->_disabled);
    if(treeElement->item())
      treeElement->item()->updateIcon();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

bool TreeElement::canDisabled() const
{
  if(disabled())
    return true;  // a disabled element can always be enabled
  return unbounded() || optional();
}

/***********************************************/

void TreeElement::setDisabled(bool disabled)
{
  try
  {
    if((disabled != _disabled) && canDisabled())
      tree->undoStack->push(new UndoCommandSetEnabled(this, disabled));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/
/***********************************************/

void TreeElement::UndoCommandUpdateName::redo()
{
  try
  {
    if(!name.isEmpty())
      std::swap(name,  treeElement->_name);
    if(!value.isEmpty())
      std::swap(value, treeElement->_valueList[treeElement->selectedIndex()]);

    if(treeElement->item())
    {
      tree->setSelectedItem(treeElement->item());
      treeElement->item()->lostCurrent();
      treeElement->item()->updateName();
      treeElement->item()->updateIcon();
      treeElement->item()->updateValue();
      treeElement->item()->becomeCurrent(); // rebuild comboBox
    }

    tree->treeChanged();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

bool TreeElement::canUpdateName() const
{
  return isRenamedInSchema();
}

/***********************************************/

void TreeElement::updateName()
{
  try
  {
    if(TreeElement::canUpdateName())
      tree->undoStack->push(new UndoCommandUpdateName(this, QString()));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/
/***********************************************/

void TreeElement::setElementAdd(TreeElementAdd *elementAdd)
{
  try
  {
    if(_elementAdd)
      _elementAdd->unboundedCount--;
    _elementAdd = elementAdd;
    if(_elementAdd)
      _elementAdd->unboundedCount++;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/
/***********************************************/

TreeItem *TreeElement::createItem(TreeItem *parent, TreeItem *after)
{
  try
  {
    if(!_item)
      _item = TreeItem::newTreeItem(this, parent, after);
    else
    {
      _item->updateIcon();
      _item->updateValue();
      _item->updateComment();
    }
    return _item;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

void TreeElement::removeItem()
{
  try
  {
    delete _item;
    _item = nullptr;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/
/***********************************************/

QComboBox *TreeElement::createComboBox(bool isEditable)
{
  try
  {
    // is there a choice?
    if(!isEditable && (_valueList.size() <= 1))
      return nullptr;

    comboBox = new QComboBox(tree);
    comboBox->setSizePolicy(QSizePolicy(QSizePolicy::Expanding, QSizePolicy::Fixed));
    comboBox->setEditable(isEditable);
    comboBox->setCompleter(nullptr);
    comboBox->setInsertPolicy(QComboBox::NoInsert);
    comboBox->installEventFilter(tree);
    comboBox->setFocusPolicy(Qt::StrongFocus);

    // fill comboBox
    for(int i=0; i<_valueCount; i++)
      comboBox->insertItem(comboBox->count(), _valueList[i]);
    for(int i=_valueCount; i<_valueList.size(); i++)
      comboBox->insertItem(comboBox->count(), QIcon(":/icons/scalable/link.svg"), _valueList[i]);
    comboBox->setCurrentIndex(_selectedIndex);

    connect(comboBox, SIGNAL(activated(int)), this, SLOT(comboBoxActivated(int)));
    if(comboBox->lineEdit())
    {
      connect(comboBox->lineEdit(), SIGNAL(textEdited(const QString &)), this, SLOT(comboBoxTextEdited(const QString &)));
      connect(tree, SIGNAL(sectionResized(int, int, int)), this, SLOT(columnResized(int, int, int)));
    }

    return comboBox;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

void TreeElement::comboBoxSetToolTip()
{
  try
  {
    if(!comboBox || !comboBox->lineEdit())
      return;

    QLineEdit *lineEdit = comboBox->lineEdit();
    if(lineEdit->width() < QFontMetrics(lineEdit->font()).boundingRect(lineEdit->text()).width()+7)
      lineEdit->setToolTip(lineEdit->text());
    else
      lineEdit->setToolTip("");
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

void TreeElement::comboBoxActivated(int index)
{
  try
  {
    changeSelectedIndex(index);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/
/***********************************************/

void TreeElement::UndoCommandEdit::redo()
{
  try
  {
    if(treeElement->item())
      tree->setSelectedItem(treeElement->item());
    treeElement->setSelectedIndex(treeElement->insertNewValue(newValue, true/*prepend*/));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

void TreeElement::UndoCommandEdit::undo()
{
  try
  {
    if(treeElement->item())
      tree->setSelectedItem(treeElement->item());

    if(isLink)
      treeElement->setSelectedIndex(treeElement->findLinkIndex(oldValue));
    else
      treeElement->setSelectedIndex(treeElement->insertNewValue(oldValue, true/*prepend*/));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

bool TreeElement::UndoCommandEdit::mergeWith(const QUndoCommand *command)
{
  try
  {
    const UndoCommandEdit *newCommand = dynamic_cast<const UndoCommandEdit*>(command);
    if(treeElement != newCommand->treeElement)
      return false;

    // remove old text
    if(isCreated && (newCommand->newValue != newValue))
    {
      int index = treeElement->findValueIndex(newValue);
      if(index>=0)
      {
        treeElement->_valueList.removeAt(index);
        treeElement->_valueCount--;
        if(treeElement->_selectedIndex >= index)
          treeElement->_selectedIndex--;

        // update comboBox
        auto comboBox = treeElement->comboBox;
        if(comboBox)
        {
          int cursorPosition = comboBox->lineEdit()->cursorPosition();
          comboBox->removeItem(index);
          comboBox->lineEdit()->setCursorPosition(cursorPosition);

          if(comboBox->count() != treeElement->_valueList.size())
            throw(Exception("comboBox not consistent"));
          if(comboBox->currentIndex() != treeElement->_selectedIndex)
            comboBox->setCurrentIndex(treeElement->_selectedIndex);
        }
      }
    }

    newValue  = newCommand->newValue;
    isCreated = newCommand->isCreated;
    setText("edit "+treeElement->name()+": "+newValue);
    return true;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

void TreeElement::comboBoxTextEdited(const QString &text)
{
  try
  {
    if(isLinked() || (text != selectedValue()))
      tree->undoStack->push(new UndoCommandEdit(this, text));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

void TreeElement::columnResized(int column, int /*oldSize*/, int /*newSize*/)
{
  try
  {
    if((column == 1) && comboBox && comboBox->lineEdit())
      comboBoxSetToolTip();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/
/***********************************************/
