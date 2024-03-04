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
#include "tree/treeElementLoopCondition.h"
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
      this->defaultValue = QString (xsdElement->defaultValue);
    if(!defaultOverride.isEmpty())
      this->defaultValue   = defaultOverride;
    this->_selectedIndex   = -1;
    this->_brokenLinkIndex = -1;
    this->_valueCount      = 0;
    this->_pushComment     = false;
    this->loop             = nullptr;
    this->condition        = nullptr;
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

      QString loopLabel;
      readAttribute(xmlNode, "loop", loopLabel);
      if(!loopLabel.isEmpty())
      {
        XmlNodePtr xmlLoop;
        if(loopLabel == "_localLoop_")
          xmlLoop = xmlNode->getChild("loopType");
        else
        {
          xmlLoop = XmlNode::create("loopType");
          writeAttribute(xmlLoop, "link", loopLabel);
        }
        loop = new TreeElementLoopCondition(tree, this, tree->xsdElementLoop, xmlLoop);
      }

      QString conditionLabel;
      readAttribute(xmlNode, "condition", conditionLabel);
      if(!conditionLabel.isEmpty())
      {
        XmlNodePtr xmlCondition;
        if(conditionLabel == "_localCondition_")
          xmlCondition = xmlNode->getChild("conditionType");
        else
        {
          xmlCondition = XmlNode::create("conditionType");
          writeAttribute(xmlCondition, "link", conditionLabel);
        }
        condition = new TreeElementLoopCondition(tree, this, tree->xsdElementCondition, xmlCondition);
      }

      readAttribute(xmlNode, "label",   _label);
      readAttribute(xmlNode, "comment", _comment);
      XmlAttrPtr attr = xmlNode->getAttribute("disabled");
      if(attr && (attr->text == "1"))
        _disabled = true;

      // is this a link?
      QString link;
      readAttribute(xmlNode, "link", link);
      if(!link.isEmpty())
      {
        _valueList.insert(_valueCount, link);
        _selectedIndex = _valueCount;
      }
    }
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
  delete loop;
  delete condition;
  setElementAdd(nullptr);
}

/***********************************************/

XmlNodePtr TreeElement::createXmlBaseNode() const
{
  try
  {
    XmlNodePtr xmlNode = XmlNode::create(xmlName());
    if(!label().isEmpty())     writeAttribute(xmlNode, "label",     label());
    if(disabled())             writeAttribute(xmlNode, "disabled",  "1");
    if(isLinked())             writeAttribute(xmlNode, "link",      selectedValue());

    if(loop)
    {
      if(loop->isLinked())
        writeAttribute(xmlNode, "loop", loop->selectedValue());
      else
      {
        writeAttribute(xmlNode, "loop", "_localLoop_");
        XmlNodePtr xmlLoop = loop->createXmlTree(true);
        xmlNode->addChild(xmlLoop);
      }
    }

    if(condition)
    {
      if(condition->isLinked())
        writeAttribute(xmlNode, "condition", condition->selectedValue());
      else
      {
        writeAttribute(xmlNode, "condition", "_localCondition_");
        XmlNodePtr xmlCondition = condition->createXmlTree(true);
        xmlNode->addChild(xmlCondition);
      }
    }

    if(!comment().isEmpty())
      writeAttribute(xmlNode, "comment",   comment());
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

void TreeElement::changeSelectedLink(const QString &link)
{
  try
  {
    int selectedIndex = findLinkIndex(link);
    if(selectedIndex < 0)
    {
      selectedIndex = _valueList.size();
      _valueList.push_back(link);
      if(comboBox)
        comboBox->insertItem(comboBox->count(), QIcon(":/icons/scalable/link.svg"), link);
    }
    changeSelectedIndex(selectedIndex);
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
    if((index < 0) || (index >= _valueList.size()))
      throw(Exception("index out of range"));

    if(index == _selectedIndex)
      return;
    _selectedIndex = index;

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
    if(_brokenLinkIndex >= index)
      _brokenLinkIndex++;

    // update comboBox
    if(comboBox)
      comboBox->insertItem(index, value);

    // empty values only allowed at the beginning
    for(int i=1; i<_valueList.size(); i++)
      if(_valueList[i].isEmpty())
      {
        if(_selectedIndex >= i)
          _selectedIndex--;
        if(_brokenLinkIndex >= i)
          _brokenLinkIndex--;
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

void TreeElement::updateParserResults(VariableList &varList)
{
  if(loop)      {VariableList varListLocal = varList; loop->updateParserResults(varListLocal);}
  if(condition) {VariableList varListLocal = varList; condition->updateParserResults(varListLocal);}
}

/***********************************************/

void TreeElement::updateLinks(QMap<QString, QString> &labelTypes)
{
  try
  {
    if(loop)      {auto labelTypesLocal = labelTypes; loop->updateLinks(labelTypesLocal);}
    if(condition) {auto labelTypesLocal = labelTypes; condition->updateLinks(labelTypesLocal);}

    if(!label().isEmpty() && !disabled())
      labelTypes[label()] = type();

    QString link;
    if(isLinked())
      link = selectedValue();

    // remove all links
    while(_valueList.size() > _valueCount)
      _valueList.removeAt(_valueCount);
    while(comboBox && (comboBox->count() > _valueCount))
      comboBox->removeItem(_valueCount);

    // refill links
    for(auto iter = labelTypes.keyValueBegin(); iter != labelTypes.keyValueEnd(); iter++)
      if((iter->first != label()) && (iter->second == type())) // should not linked to self
      {
        _valueList.push_back(iter->first);
        if(comboBox)
          comboBox->insertItem(comboBox->count(), QIcon(":/icons/scalable/link.svg"), iter->first);
      }

    // restore selected link
    _brokenLinkIndex = -1;
    if(!link.isEmpty())
    {
      int selectedIndex = findLinkIndex(link);
      if(selectedIndex < 0)
      {
        selectedIndex = _valueList.size();
        _brokenLinkIndex = selectedIndex;
        _valueList.push_back(link);
        if(comboBox)
          comboBox->insertItem(comboBox->count(), QIcon(":/icons/scalable/link-broken.svg"), link);
      }
      setSelectedIndex(selectedIndex);
      if(comboBox && (comboBox->currentIndex() != selectedIndex))
      {
        comboBox->setCurrentIndex(selectedIndex);
        comboBoxSetToolTip();
      }

      if(item())
        item()->updateValue();
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}


/***********************************************/
/***********************************************/

bool TreeElement::renamedLink(const QString &oldLabel, const QString &newLabel)
{
  try
  {
    // Is this a new variable with same name (hide scope of renamed)
    if(label() == oldLabel)
      return false;

    if(loop)      loop->renamedLink(oldLabel, newLabel);
    if(condition) condition->renamedLink(oldLabel, newLabel);

    // is this link in the list?
    int index = findLinkIndex(oldLabel);
    if(index < 0)
      return true;

    // replace label
    _valueList.replace(index, newLabel);
    if(item())
      item()->updateValue();
    if(comboBox && (index < comboBox->count()))
      comboBox->setItemText(index, newLabel);

    return true;
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
    QString loopLabel, conditionLabel, comment;
    readAttribute(xmlNode, "loop",      loopLabel);
    readAttribute(xmlNode, "condition", conditionLabel);
    readAttribute(xmlNode, "comment",   comment);
    setComment(comment);

    if(!loopLabel.isEmpty())
    {
      XmlNodePtr xmlLoop;
      if(loopLabel == "_localLoop_")
        xmlLoop = xmlNode->getChild("loopType");
      else
      {
        xmlLoop = XmlNode::create("loopType");
        writeAttribute(xmlLoop, "link", loopLabel);
      }
      if(canSetLoop())
        loop = new TreeElementLoopCondition(tree, this, tree->xsdElementLoop, xmlLoop);
    }

    if(!conditionLabel.isEmpty())
    {
      XmlNodePtr xmlCondition;
      if(conditionLabel == "_localCondition_")
        xmlCondition = xmlNode->getChild("conditionType");
      else
      {
        xmlCondition = XmlNode::create("conditionType");
        writeAttribute(xmlCondition, "link", conditionLabel);
      }
      if(canSetCondition())
        condition = new TreeElementLoopCondition(tree, this, tree->xsdElementCondition, xmlCondition);
    }
  }

  QString link;
  readAttribute(xmlNode, "link", link);
  if(!link.isEmpty())
    changeSelectedLink(link);
  return link.isEmpty();
}

/***********************************************/
/***********************************************/

class TreeElement::UndoCommandSetLoop : public TreeElement::UndoCommand
{
  TreeElementLoopCondition *loop;

public:
  UndoCommandSetLoop(TreeElement *treeElement, TreeElementLoopCondition *loop)
    : UndoCommand(treeElement, (loop ? "set loop for" : "remove loop from")),
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

    std::swap(loop, treeElement->loop);

    if(treeElement->loop && treeElement->parentElement)
    {
      treeElement->parentElement->updateParserResultsInScope();
      treeElement->parentElement->updateLinksInScope();
      tree->treeChanged();
    }

    if(treeElement->item())
    {
      if(treeElement->loop)
      {
        treeElement->loop->createItem(treeElement->item(), nullptr/*after*/);
        tree->setSelectedItem(treeElement->loop->item());
      }
      else
        loop->removeItem();

      treeElement->item()->updateIcon();
      treeElement->item()->updateName();
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

bool TreeElement::canSetLoop() const
{
  if(loop || !label().isEmpty())
    return false;
  return unbounded();
}

/***********************************************/

void TreeElement::setLoop()
{
  try
  {
    if(canSetLoop())
      tree->undoStack->push(new UndoCommandSetLoop(this, new TreeElementLoopCondition(tree, this, tree->xsdElementLoop, nullptr)));
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
  TreeElementLoopCondition *condition;

public:
  UndoCommandSetCondition(TreeElement *treeElement, TreeElementLoopCondition *condition)
    : UndoCommand(treeElement, (condition ?  "set condition for" : "remove condition from")),
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

    std::swap(condition, treeElement->condition);

    if(treeElement->condition && treeElement->parentElement)
    {
      treeElement->parentElement->updateParserResultsInScope();
      treeElement->parentElement->updateLinksInScope();
      tree->treeChanged();
    }

    if(treeElement->item())
    {
      if(treeElement->condition)
      {
        treeElement->condition->createItem(treeElement->item(), (treeElement->loop ? treeElement->loop->item() : nullptr)/*after*/);
        tree->setSelectedItem(treeElement->condition->item());
      }
      else
        condition->removeItem();

      treeElement->item()->updateIcon();
      treeElement->item()->updateName();
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

bool TreeElement::canSetCondition() const
{
  if(condition || !label().isEmpty())
    return false;
  return optional() || unbounded();
}

/***********************************************/

void TreeElement::setCondition()
{
  try
  {
    if(canSetCondition())
      tree->undoStack->push(new UndoCommandSetCondition(this, new TreeElementLoopCondition(tree, this, tree->xsdElementCondition, nullptr)));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

void TreeElement::removeLoopOrCondition(TreeElement *element)
{
  try
  {
    if(element == loop)
      tree->undoStack->push(new UndoCommandSetLoop(this, nullptr));
    if(element == condition)
      tree->undoStack->push(new UndoCommandSetCondition(this, nullptr));
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

    if(!treeElement->label().isEmpty() && treeElement->parentElement)
    {
      treeElement->parentElement->updateParserResultsInScope();
      treeElement->parentElement->updateLinksInScope();
    }

    if(treeElement->item())
      treeElement->item()->updateIcon();

    tree->treeChanged();
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

class TreeElement::UndoCommandRename : public TreeElement::UndoCommand
{
  QString label;

public:
  UndoCommandRename(TreeElement *treeElement, QString label)
  : UndoCommand(treeElement, "rename"), label(label)
  {
    setText("rename "+treeElement->name()+" to "+label);
  }

  void redo();
  void undo() {redo();}
};

/***********************************************/

void TreeElement::UndoCommandRename::redo()
{
  try
  {
    std::swap(label, treeElement->_label);

    // renames links in other elements
    auto parent = treeElement->parentElement;
    if(parent)
    {
      const QVector<TreeElement*> &children = parent->children();
      bool stopped = false;
      for(int i=children.indexOf(treeElement)+1; i<children.size(); i++)
        if(!children.at(i)->renamedLink(label, treeElement->_label))
          {stopped = true; break;}

      // variables in global section are globally valid -> repeat with parent of global
      if(!stopped && dynamic_cast<TreeElementGlobal*>(parent) && parent->parentElement)
      {
        parent = parent->parentElement; // parent of global -> groops
        const QVector<TreeElement*> &children = parent->children();
        for(int i=children.indexOf(treeElement->parentElement)+1; i<children.size(); i++) // start after global
          if(!children.at(i)->renamedLink(label, treeElement->_label))
            break;
      }

      parent->updateParserResultsInScope();
      parent->updateLinksInScope();
      tree->treeChanged();
    }

    if(treeElement->item())
    {
      tree->setSelectedItem(treeElement->item());
      treeElement->item()->updateName();
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

bool TreeElement::canRename() const
{
  return !label().isEmpty();
}

/***********************************************/

bool TreeElement::rename(const QString &label)
{
  try
  {
    if(!canRename() || label.isEmpty() || (label == this->label()))
      return false;
    tree->undoStack->push(new UndoCommandRename(this, label));
    return true;
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

    if(loop)      loop->createItem(_item, nullptr/*after*/);
    if(condition) condition->createItem(_item, (loop ? loop->item() : nullptr)/*after*/);

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
    if(!_item)
      return;
    if(loop)      loop->removeItem();
    if(condition) condition->removeItem();
    auto parentItem = dynamic_cast<QTreeWidgetItem*>(_item)->parent();
    if(parentItem)
      parentItem->removeChild(_item);
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
      comboBox->insertItem(comboBox->count(), ((i == _brokenLinkIndex) ? QIcon(":/icons/scalable/link-broken.svg") : QIcon(":/icons/scalable/link.svg")), _valueList[i]);
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
    tree->treeChanged();
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
    tree->treeChanged();
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
        if(treeElement->_brokenLinkIndex >= index)
          treeElement->_brokenLinkIndex--;

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
