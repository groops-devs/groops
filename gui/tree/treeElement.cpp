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
                                         XmlNodePtr xmlNode, Bool fromFile)
{
  try
  {
    if((tree==nullptr)||(xsdElement==nullptr))
      throw(Exception("Null Pointer"));

    if(xsdElement->complex)
    {
      if(xsdElement->name=="global")
        return new TreeElementGlobal(tree, parentElement, xsdElement, xmlNode);
      else if(xsdElement->name=="programme" || xsdElement->name=="program")
        return new TreeElementProgram(tree, parentElement, xsdElement, defaultOverride, xmlNode, fromFile);
      else if(xsdElement->complex->type=="sequence")
        return new TreeElementSequence(tree, parentElement, xsdElement, defaultOverride, xmlNode, fromFile);
      else if(xsdElement->complex->type=="choice")
        return new TreeElementChoice(tree, parentElement, xsdElement, defaultOverride, xmlNode, fromFile);
      else
        throw(Exception("unknown complex element"));
    }

    if(xsdElement->type=="boolean")
      return new TreeElementBool(tree, parentElement, xsdElement, defaultOverride, xmlNode);
    else if(xsdElement->type=="filename")
      return new TreeElementFileName(tree, parentElement, xsdElement, defaultOverride, xmlNode, fromFile);
    else if(xsdElement->type=="time")
      return new TreeElementTime(tree, parentElement, xsdElement, defaultOverride, xmlNode, fromFile);

    return new TreeElementSimple(tree, parentElement, xsdElement, defaultOverride, xmlNode, fromFile);
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
    this->_defaultOverride = defaultOverride;
    this->_selectedIndex   = -1;
    this->_valueCount      = 0;
    this->_pushComment     = false;
    this->_disabled        = false;
    this->_elementAdd      = nullptr;
    this->_item            = nullptr;
    this->comboBox         = nullptr;

    if(xsdElement != nullptr)
      _name = xsdElement->name;

    // evaluation of the XML nodes
    // ---------------------------
    if(xmlNode != nullptr)
    {
      if(_name.isEmpty())
        _name = xmlNode->getName();

      readAttribute(xmlNode, "label",     _label);
      readAttribute(xmlNode, "loop",      _loop);
      readAttribute(xmlNode, "condition", _condition);
      readAttribute(xmlNode, "comment",   _comment);

      XmlAttrPtr attr = xmlNode->getAttribute("disabled");
      if((attr!=nullptr) && (attr->getText() == "1"))
        _disabled = true;

      // is this a link?
      attr = xmlNode->getAttribute("link");
      if(attr!=nullptr)
      {
        _valueList.insert(_valueCount, attr->getText());
        _selectedIndex = _valueCount;
      }
    }

    if(_name.isEmpty())
      _name = "unknown";
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

TreeElement::~TreeElement()
{
  try
  {
    removeItem();
    setElementAdd(nullptr);
  }
  catch(std::exception &e)
  {
    qDebug() << QString::fromStdString("Exception in destructor at "+_GROOPS_ERRORLINE+"\n"+e.what());
  }
}

/***********************************************/

XmlNodePtr TreeElement::getBaseXML() const
{
  try
  {
    if(_name.isEmpty())
      return XmlNodePtr(nullptr);
    XmlNodePtr xmlNode(new XmlNode(isRenamed() ? originalName() : _name));
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

QString TreeElement::parseExpression(const QString &value) const
{
  if(xsdElement==nullptr || isElementAdd() || value.isEmpty())
    return QString();

  QString result;
  try
  {
    if(parentElement && parentElement == tree->elementGlobal())
      for(auto &&var : tree->setVarList())
        var.second->setValue(var.second->getText()); // reset as TEXT in case variables are still set CIRCULAR

    Bool resolved;
    QString text = (findLinkIndex(value) >= 0 /*&& !isComplex()*/) ? "{"+value+"}" : value;
    text = QString::fromStdString(StringParser::parse(name().toStdString(), text.toStdString(), tree->varList(), resolved));
    result = text;

    if(QStringList({"bool", "int", "uint", "double", "angle", "time", "expression"}).contains(type())) // only numerical values
    {
      Double d = Expression::parse(text.toStdString())->evaluate(tree->varList());
      result.setNum(d, 'f', 7).remove(QRegExp("0+$")).remove(QRegExp("\\.$")); // %.7f with trailing zeros removed
    }

    return result;
  }
  catch(...)
  {
    return result;
  }
}

/***********************************************/

int TreeElement::findValueIndex(const QString &value) const
{
  try
  {
    for(int i=_valueCount; i-->0; ) // start from behind to find potential renames first
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
  Bool    isLink;
  QString text;

public:
  UndoCommandChangeSelectedIndex(TreeElement *treeElement, int index)
    : UndoCommand(treeElement, "change"),
      isLink(index >= treeElement->_valueCount),
      text(treeElement->_valueList[index]) {setText("change "+treeElement->name()+" from "+treeElement->selectedValue()+" to "+treeElement->_valueList[index]);}

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

    // swap node
    Bool    tmpIsLink = treeElement->isLinked();
    QString tmpText   = treeElement->selectedValue();

    if(isLink)
      treeElement->setSelectedIndex(treeElement->findLinkIndex(text));
    else
      treeElement->setSelectedIndex(treeElement->insertNewValue(text));

    isLink = tmpIsLink;
    text   = tmpText;

    treeElement->trackRenamed(treeElement->isSelectionRenamed(treeElement->selectedIndex()));
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

    tree->pushUndoCommand(new UndoCommandChangeSelectedIndex(this, index));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

void TreeElement::changeSelectedValue(const QString &value)
{
  try
  {
    changeSelectedIndex(insertNewValue(value, true/*prepend*/));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/
/***********************************************/

void TreeElement::setSelectedIndex(int index)
{
  try
  {
    if((index<0) || (index>=_valueList.size()))
      throw(Exception("index out of range"));

    if(comboBox)
      comboBox->setCurrentIndex(index);

    if(index == _selectedIndex)
      return;
    _selectedIndex = index;

    // send auto comment
    if(_pushComment && parentElement)
      parentElement->setAutoComment(selectedValue());

    if(item())
      item()->updateValue();

    // call the handler
    newSelectedIndex(index);

    // update variable list
    if(!label().isEmpty())
    {
      if(isLinked())
        tree->setVarList().addVariable(ExpressionVariablePtr(new ExpressionVariable(label().toStdString(), "{"+selectedValue().toStdString()+"}")));
      else
        tree->setVarList().addVariable(ExpressionVariablePtr(new ExpressionVariable(label().toStdString(), selectedValue().toStdString())));
      if(tree->rootElement())
      {
        tree->rootElement()->newLink(this);
        tree->updateExpressions(this);
      }
    }

    if(comboBox && (comboBox->count() != _valueList.size()))
      throw(Exception("comboBox not consistent"));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

int TreeElement::insertNewValue(const QString &value, Bool prepend)
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
    if(_valueCount > 0 && _valueList[0].isEmpty())
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

void TreeElement::newLink(TreeElement *element)
{
  try
  {
    if((!element) || (element==this) || (element->type()!=type() && !selectedValue().contains(element->name())))
      return;

    // is this link already in the list?
    if(findLinkIndex(element->label()) >= 0)
      return;

    // new link
    if(element->type()==type())
    {
      _valueList.push_back(element->label());
      if(comboBox)
        comboBox->insertItem(comboBox->count(), QIcon(":/icons/scalable/link.svg"), element->label());
    }

    if(comboBox && (comboBox->count() != _valueList.size()))
      throw(Exception("comboBox not consistent"));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

void TreeElement::removeLink(TreeElement *element)
{
  try
  {
    if((!element) || (element==this))
      return;

    // is this link in the list?
    int index = findLinkIndex(element->label());
    if(index<0)
      return;

    if(_selectedIndex >= index)
      _selectedIndex--;

    _valueList.removeAt(index);
    if(comboBox)
      comboBox->removeItem(index);

    if(comboBox && (comboBox->count() != _valueList.size()))
      throw(Exception("comboBox not consistent"));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

void TreeElement::renameLink(const QString &oldLabel, const QString &newLabel)
{
  try
  {
    if(oldLabel == newLabel)
      return;

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

void TreeElement::trackUnknown(Bool track)
{
  try
  {
    if(!isUnknown() || !tree)
      return;

    if(track)
      tree->trackUnknownElement(this);
    else
      tree->untrackUnknownElement(this);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

void TreeElement::trackRenamed(Bool track)
{
  try
  {
    if(!tree)
      return;

    if(track)
      tree->trackRenamedElement(this);
    else
      tree->untrackRenamedElement(this);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

void TreeElement::getLinkedList(const QString &label, QList<TreeElement*> &list)
{
  try
  {
    if(isLinked() && (selectedValue() == label))
      list.push_back(this);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

void TreeElement::getLoopList(const QString &loop, QList<TreeElement*> &list)
{
  try
  {
    if(_loop == loop)
      list.push_back(this);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

void TreeElement::getConditionList(const QString &condition, QList<TreeElement*> &list)
{
  try
  {
    if(_condition == condition)
      list.push_back(this);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

void TreeElement::updateExpression()
{
  try
  {
    if(item())
      item()->updateValue();
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

    // swap text
    QString tmp = treeElement->_comment;
    treeElement->_comment = text;
    text = tmp;

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
    if(_comment == text) // is something changed?
      return;

    tree->pushUndoCommand(new UndoCommandSetComment(this, text));
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
    if((_comment == _autoComment)||(_comment == text))
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

void TreeElement::setPushAutoComments(Bool on)
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

class TreeElement::UndoCommandSetEnabled : public TreeElement::UndoCommand
{
  bool state;

public:
  UndoCommandSetEnabled(TreeElement *treeElement, bool _state)
    : UndoCommand(treeElement, (_state ? "disable" : "enable")),
      state(_state) {}

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

    treeElement->_disabled = state;
    state = !state;
    if(treeElement->item())
      treeElement->item()->updateDisabled();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

Bool TreeElement::canDisabled() const
{
  if(_disabled)
    return true;  // a disabled element can always be enabled
  if(isElementAdd())
    return false; // the add button cannot be disabled
  if(optional())
    return true;
  return unbounded();
}

/***********************************************/

void TreeElement::setDisabled(Bool state)
{
  try
  {
    if(state == _disabled)
      return;
    if(!canDisabled())
      return;
    tree->pushUndoCommand(new UndoCommandSetEnabled(this, state));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/
/***********************************************/

class TreeElement::UndoCommandSetLoop : public TreeElement::UndoCommand
{
  QString loop;

public:
  UndoCommandSetLoop(TreeElement *treeElement, QString _loop)
    : UndoCommand(treeElement, (_loop.isEmpty() ? "remove loop " + treeElement->loop() + " from" : "set loop " + _loop + " for")),
      loop(_loop) {}

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

Bool TreeElement::canSetLoop() const
{
  if(isElementAdd())
    return false; // loops cannot be set for the add button
  if(parentElement == tree->elementGlobal())
    return false;
  return unbounded();
}

/***********************************************/

void TreeElement::setLoop(const QString &loop)
{
  try
  {
    if(loop == _loop)
      return;
    if(!canSetLoop() && !loop.isEmpty())
      return;
    tree->pushUndoCommand(new UndoCommandSetLoop(this, loop));
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
  UndoCommandSetCondition(TreeElement *treeElement, QString _condition)
    : UndoCommand(treeElement, (_condition.isEmpty() ? "remove condition " + treeElement->condition() + " from" : "set condition " + _condition + " for")),
      condition(_condition) {}

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

Bool TreeElement::canSetCondition() const
{
  if(isElementAdd())
    return false; // conditions cannot be set for the add button
  if(parentElement == tree->elementGlobal())
    return false;
  return optional() || unbounded();
}

/***********************************************/

void TreeElement::setCondition(const QString &condition)
{
  try
  {
    if(condition == _condition)
      return;
    if(!canSetCondition() && !condition.isEmpty())
      return;
    tree->pushUndoCommand(new UndoCommandSetCondition(this, condition));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/
/***********************************************/

class TreeElement::UndoCommandRename : public TreeElement::UndoCommand
{
  QString label;

public:
  UndoCommandRename(TreeElement *treeElement, QString _label)
    : UndoCommand(treeElement, "rename")
  {
    label = _label;
    setText("rename "+treeElement->name()+" to "+_label);
  }

  void redo();
  void undo() {redo();}
};

/***********************************************/

void TreeElement::UndoCommandRename::redo()
{
  try
  {
    if(treeElement->item())
      tree->setSelectedItem(treeElement->item());

    if(treeElement->parentElement == tree->elementGlobal())
    {
      QString oldLabel = treeElement->label();
      tree->rootElement()->renameLink(oldLabel, label);
      treeElement->setLabel(label);
      label = oldLabel;

      // replace varList entry
      tree->setVarList().erase(oldLabel.toStdString());
      if(treeElement->isLinked())
        tree->setVarList().addVariable(ExpressionVariablePtr(new ExpressionVariable(treeElement->label().toStdString(), "{"+treeElement->selectedValue().toStdString()+"}")));
      else
        tree->setVarList().addVariable(ExpressionVariablePtr(new ExpressionVariable(treeElement->label().toStdString(), treeElement->selectedValue().toStdString())));

      // replace link entries
      for(auto &&entry : tree->setVarList())
        if(entry.second->getText() == "{"+oldLabel.toStdString()+"}")
          entry.second->setValue("{"+treeElement->label().toStdString()+"}");

      if(tree->rootElement())
        tree->rootElement()->updateExpression();
    }

    if(treeElement->item())
    {
      treeElement->item()->updateName();
      treeElement->item()->updateDisabled();
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

Bool TreeElement::canRename() const
{
  if(isElementAdd())
    return false; // the add button cannot be renamed
  if(parentElement == tree->elementGlobal())
    return true; // global elements can always be renamed
  return false;
}

/***********************************************/

void TreeElement::rename(QString label)
{
  try
  {
    if(!canRename())
      return;
    if(label != this->label() && parentElement == tree->elementGlobal())
      tree->pushUndoCommand(new UndoCommandRename(this, label));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/
/***********************************************/

class TreeElement::UndoCommandUpdateName : public TreeElement::UndoCommand
{
  QString label;
  Bool renameSelection;

public:
  UndoCommandUpdateName(TreeElement *treeElement, Bool renameSelection)
    : UndoCommand(treeElement, "update name")
  {
    this->renameSelection = renameSelection;

    if(renameSelection && treeElement->isSelectionRenamed(treeElement->selectedIndex()))
    {
      TreeElementChoice *element = dynamic_cast<TreeElementChoice*>(treeElement);
      label = element->renamedSelection();
      setText("update "+treeElement->name()+" choice name from "+treeElement->selectedValue()+" to "+label);
    }
    else if(treeElement->originalName() != treeElement->name())
    {
      label = treeElement->originalName();
      if(treeElement->parentElement->isElementGlobal())
        setText("update "+treeElement->name()+" type from "+treeElement->originalName()+" to "+treeElement->type());
      else
        setText("update name of "+treeElement->originalName()+" to "+treeElement->name());
    }
  }

  void redo();
  void undo() {redo();}
};

/***********************************************/

void TreeElement::UndoCommandUpdateName::redo()
{
  try
  {
    if(treeElement->item())
      tree->setSelectedItem(treeElement->item());

    if(renameSelection)
    {
      QString oldSelection = treeElement->selectedValue();
      treeElement->_valueList[treeElement->selectedIndex()] = label;
      label = oldSelection;
      if(treeElement->item())
      {
        treeElement->item()->lostCurrent();
        treeElement->item()->updateValue();
        treeElement->item()->becomeCurrent();
      }
    }
    else // rename element name
    {
      if(label != treeElement->name()) // renamed element
        treeElement->setOriginalName(!treeElement->isRenamed() ? label : "");

      if(treeElement->item())
      {
        treeElement->item()->updateName();
        treeElement->item()->updateDisabled();
      }
    }

    treeElement->trackRenamed(treeElement->isRenamed() || treeElement->isSelectionRenamed(treeElement->selectedIndex()));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

Bool TreeElement::canUpdateName() const
{
  if(isElementAdd())
    return false; // the add button cannot be renamed
  if(isRenamed())
    return true; // renamed elements can be renamed
  if(isSelectionRenamed(selectedIndex()))
    return true; // renamed selections can be renamed

  return false;
}

/***********************************************/

void TreeElement::updateName()
{
  try
  {
    if(!canUpdateName())
      return;
    if(isRenamed())
      tree->pushUndoCommand(new UndoCommandUpdateName(this, false));
    if(isSelectionRenamed(selectedIndex()))
      tree->pushUndoCommand(new UndoCommandUpdateName(this, true));
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
      _item->updateDisabled();
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

QComboBox *TreeElement::createComboBox(Bool isEditable)
{
  try
  {
    // is there a choice?
    if((_valueList.size()<=1) && (!isEditable))
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
      connect(tree->header(), SIGNAL(sectionResized(int, int, int)), this, SLOT(columnResized(int, int, int)));
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
    if(!comboBox || (comboBox && !comboBox->lineEdit()))
      return;

    QLineEdit* lineEdit = comboBox->lineEdit();
    QFontMetrics fm(lineEdit->font());

    if(lineEdit->width() < fm.boundingRect(lineEdit->text()).width()+7)
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
    comboBoxSetToolTip();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/
/***********************************************/

class TreeElement::UndoCommandEdit : public TreeElement::UndoCommand
{
  Bool    isCreated;
  QString newValue;
  Bool    isLink;
  QString oldValue;

public:
  UndoCommandEdit(TreeElement *treeElement, const QString &text)
    : UndoCommand(treeElement, "edit"),
      isCreated(treeElement->findValueIndex(text)<0),
      newValue(text),
      isLink(treeElement->isLinked()),
      oldValue(treeElement->selectedValue())
    {setText("edit "+treeElement->name()+": "+text);}

  int  id() const {return 999;}
  void redo();
  void undo();
  bool mergeWith(const QUndoCommand *command);
};

/***********************************************/

void TreeElement::UndoCommandEdit::redo()
{
  try
  {
    if(treeElement->item())
      tree->setSelectedItem(treeElement->item());

    treeElement->setSelectedIndex(treeElement->insertNewValue(newValue, true/*prepend*/));

    tree->updateExpressions(treeElement);
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
      treeElement->setSelectedIndex(treeElement->insertNewValue(oldValue));

    tree->updateExpressions(treeElement);
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
        if(treeElement->comboBox)
        {
          int cursorPosition = treeElement->comboBox->lineEdit()->cursorPosition();
          treeElement->comboBox->removeItem(index);
          treeElement->comboBox->lineEdit()->setCursorPosition(cursorPosition);

          if(treeElement->comboBox->count() != treeElement->_valueList.size())
            throw(Exception("comboBox not consistent"));
          if(treeElement->comboBox->currentIndex() != treeElement->_selectedIndex)
            treeElement->comboBox->setCurrentIndex(treeElement->_selectedIndex);
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
    comboBoxSetToolTip();

    if(isLinked() || (text != selectedValue()))
      tree->pushUndoCommand(new UndoCommandEdit(this, text));
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
    if(column == 1 && comboBox && comboBox->lineEdit())
      comboBoxSetToolTip();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/
/***********************************************/
