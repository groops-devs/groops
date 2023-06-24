/***********************************************/
/**
* @file treeElementComplex.cpp
*
* @brief Abstract element with children (sequence or choice).
*
* @author Torsten Mayer-Guerr
* @author Sebastian Strasser
* @date 2006-10-07
*/
/***********************************************/

#include <QtDebug>
#include <QJsonDocument>
#include <QJsonObject>
#include <QJsonArray>
#include "base/importGroops.h"
#include "base/xml.h"
#include "base/schema.h"
#include "tree/tree.h"
#include "tree/treeItem.h"
#include "tree/treeElement.h"
#include "tree/treeElementGlobal.h"
#include "tree/treeElementAdd.h"
#include "tree/treeElementUnknown.h"
#include "tree/treeElementComplex.h"

/***********************************************/

TreeElementComplex::TreeElementComplex(Tree *tree, TreeElementComplex *parentElement,
                                       XsdElementPtr xsdElement, const QString &defaultOverride,
                                       XmlNodePtr xmlNode, Bool recieveAutoComments, Bool isElementGlobal)
                   : TreeElement(tree, parentElement, xsdElement, defaultOverride, xmlNode)
{
  this->recieveAutoComments = recieveAutoComments;
  this->_isElementGlobal    = isElementGlobal;

  // convert default string to QJsonObject
  if(!defaultValue().isEmpty())
  {
    QJsonDocument doc = QJsonDocument::fromJson(defaultValue().toUtf8());
    if(!doc.isNull() && doc.isObject())
      defaultObject = doc.object();
    else
      defaultObject = QJsonDocument::fromJson(QString("{\""+defaultValue()+"\" : {}}").toUtf8()).object();
  }
}

/***********************************************/

TreeElementComplex::~TreeElementComplex()
{
  try
  {
    for(int index=0; index<children.size(); index++)
      for(int i=0; i<children[index].size(); i++)
        delete children[index][i];
  }
  catch(std::exception &e)
  {
    qDebug() << QString::fromStdString("Exception in destructor at "+_GROOPS_ERRORLINE+"\n"+e.what());
  }
}

/***********************************************/

Bool TreeElementComplex::isSelectionUnknown(int index) const
{
  try
  {
    if(index >= 0 && !isLinked() && index < xsdElementList.size() && (!optional() || index > 0) && !xsdElementList.at(index))
      return true;
    return false;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

void TreeElementComplex::getChildrenXML(XmlNodePtr xmlNode, bool withEmptyNodes) const
{
  try
  {
    XmlNodePtr xmlChild;
    if((xmlNode) && hasChildren(selectedIndex()))
      for(int i=0; i<children[selectedIndex()].size(); i++)
      {
        bool getEmptyNodes = (withEmptyNodes && dynamic_cast<TreeElementComplex*>(children[selectedIndex()][i]) == nullptr) ? true : false; // only for noncomplex elements
        if((xmlChild = children[selectedIndex()][i]->getXML(getEmptyNodes)))
          xmlNode->addChild(xmlChild);
      }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/
/***********************************************/

Bool TreeElementComplex::hasChildren(int index) const
{
  try
  {
    if((index<0) || (index>=children.size()))
      return false;
    return (!children[index].empty());
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

int TreeElementComplex::addChoice(const QString &value, XsdElementPtr xsdComplex, const QJsonObject &defaultObject)
{
  try
  {
    int index = insertNewValue(value, false);
    if(index!=xsdElementList.size())
      throw(Exception("lists not consistent"));
    xsdElementList.push_back(xsdComplex);
    overrideDefaultObjects.push_back(defaultObject);
    children.push_back(QVector<TreeElement*>());
    return index;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

void TreeElementComplex::createChildrenElements(int index, XmlNodePtr xmlNode)
{
  try
  {
    if((index<0) || (index>=xsdElementList.size()))
      return;

    if(xmlNode && hasChildren(index))
      throw(Exception("cannot create"));

    // children not created yet?
    if(!hasChildren(index))
    {
      XsdComplexPtr xsdComplex = nullptr;
      if(xsdElementList[index])
        xsdComplex = xsdElementList[index]->complex;
      if(xsdComplex)
      {
        if(xsdComplex->type != "sequence")
          throw(Exception("complex must be sequence"));

        for(UInt k=0; k<xsdComplex->element.size(); k++)
        {
          XsdElementPtr xsdElement = xsdComplex->element.at(k);
          QJsonValue overrideDefaultValue = overrideDefaultObjects[index].value(xsdElement->name);
          if((overrideDefaultValue.isUndefined() || overrideDefaultValue.isNull()) && xsdElement->unbounded && xsdElement->defaultValue.startsWith("["))
            overrideDefaultValue = QJsonDocument::fromJson(xsdElement->defaultValue.toUtf8()).array(); // catch arrays in default

          QStringList names({xsdElement->name});
          auto findRename = [&](const QString name){ return std::find_if(xsdComplex->renames.begin(), xsdComplex->renames.end(), [&](auto pair){ return pair.second == name; }); };
          auto iterRenames = findRename(names.back());
          while(iterRenames != xsdComplex->renames.end())
          {
            names.push_back(iterRenames->first);
            iterRenames = findRename(names.back());
          }

          // lambda function to add new tree element as child
          // -------------------------------------------------
          auto addNewChildElement = [&] (TreeElementAdd *elementAdd, const QJsonValue &defaultValue, Bool fromFile)
          {
            XmlNodePtr xmlChild;
            if(xmlNode!=nullptr)
              xmlChild = xmlNode->getChild(names);

            QString defaultStr;
            switch(defaultValue.type())
            {
              case QJsonValue::String: defaultStr = defaultValue.toString();                  break;
              case QJsonValue::Bool:   defaultStr = (defaultValue.toBool() ? "1" : "0");      break;
              case QJsonValue::Double: defaultStr = QString::number(defaultValue.toDouble()); break;
              case QJsonValue::Object: if(!defaultValue.toObject().isEmpty()) defaultStr = QJsonDocument(defaultValue.toObject()).toJson(QJsonDocument::Compact); break;
              case QJsonValue::Array:  if(!defaultValue.toArray().isEmpty())  defaultStr = QJsonDocument(defaultValue.toArray()).toJson(QJsonDocument::Compact);  break;
              default: break;
            }

            TreeElement *treeElement;
            try
            {
              treeElement = TreeElement::newTreeElement(tree, this, xsdElement, defaultStr, xmlChild, fromFile);
            }
            catch(...)
            {
              // add empty element from schema and add xmlNode as an unknown element afterwards
              treeElement = TreeElement::newTreeElement(tree, this, xsdElement, "", XmlNodePtr(nullptr), fromFile);
              xmlNode->addChild(xmlChild);
            }
            if(treeElement && xmlChild && xsdElement->name != xmlChild->getName())
              treeElement->setOriginalName(xmlChild->getName());
            if(treeElement && elementAdd)
              treeElement->setElementAdd(elementAdd);
            if(treeElement)
              treeElement->trackRenamed(true);
            if(treeElement)
              children[index].push_back(treeElement);
            else
              throw(Exception("cannot add child element"));
          };
          // -------------------------------------------------

          if(!xsdElement->unbounded)
            addNewChildElement(nullptr, overrideDefaultValue, xmlNode!=nullptr);
          else
          {
            // the add element
            TreeElementAdd *elementAdd = new TreeElementAdd(tree, this, xsdElement, !isElementGlobal()/*==visible*/);
            if(xmlNode)
            {
              UInt childCount = xmlNode->getChildCount(names);
              for(UInt i=0; i<childCount; i++)
                addNewChildElement(elementAdd, overrideDefaultValue, true);
            }
            else if(overrideDefaultValue.isArray()) // create multiple elements from default array
            {
              for(auto arrayElement : overrideDefaultValue.toArray())
                addNewChildElement(elementAdd, arrayElement, false);
            }
            else if(!xsdElement->optional) // not optional -> create at least one element
              addNewChildElement(elementAdd, overrideDefaultValue, false);

            children[index].push_back(elementAdd);
          } // if(unbounded)
        } // for(xsdElements)
      } // if(xsdComplex)

      // are there unknown XML nodes?
      // ----------------------------
      if(xmlNode && xmlNode->hasChildren())
      {
        // invisible add element
        TreeElementAdd *elementAdd = new TreeElementAdd(tree, this, XsdElementPtr(nullptr), false/*visible*/);
        while(xmlNode->hasChildren())
        {
          TreeElementUnknown *treeElement = new TreeElementUnknown(tree, this, xmlNode->getNextChild());
          treeElement->setElementAdd(elementAdd);
          children[index].push_back(treeElement);
          treeElement->trackUnknown(true);
        }
        children[index].push_back(elementAdd);
      }

      // inform the new elements about all links
      // ---------------------------------------
      if(tree->elementGlobal())
        for(int i=0; i<children[index].size(); i++)
          tree->elementGlobal()->informAboutGlobalElements(children[index][i]);
    } // if(!hasChildren)

    // set autocomment for the first element
    // -------------------------------------
    if(recieveAutoComments && hasChildren(index))
      children[index][0]->setPushAutoComments(true);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

void TreeElementComplex::newSelectedIndex(int index)
{
  try
  {
    removeChildrenItems();
    if((index<0) || (index>=xsdElementList.size()))
      return;

    if(!hasChildren(index))
      createChildrenElements(index);

    // create items
    if(item())
    {
      createChildrenItems();
      item()->setExpanded(true);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/
/***********************************************/

void TreeElementComplex::newLink(TreeElement *elementInGlobal)
{
  try
  {
    TreeElement::newLink(elementInGlobal);
    for(int index=0; index<children.size(); index++)
      for(int i=0; i<children[index].size(); i++)
        if(children[index][i])
          children[index][i]->newLink(elementInGlobal);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

void TreeElementComplex::removeLink(TreeElement *elementInGlobal)
{
  try
  {
    TreeElement::removeLink(elementInGlobal);
    for(int index=0; index<children.size(); index++)
      for(int i=0; i<children[index].size(); i++)
        if(children[index][i])
          children[index][i]->removeLink(elementInGlobal);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

void TreeElementComplex::renameLink(const QString &oldLabel, const QString &newLabel)
{
  try
  {
    if(oldLabel == newLabel)
      return;

    TreeElement::renameLink(oldLabel, newLabel);
    for(int index=0; index<children.size(); index++)
      for(int i=0; i<children[index].size(); i++)
        if(children[index][i])
          children[index][i]->renameLink(oldLabel, newLabel);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

void TreeElementComplex::trackUnknown(Bool track)
{
  try
  {
    TreeElement::trackUnknown(track);
    for(int index=0; index<children.size(); index++)
      for(int i=0; i<children[index].size(); i++)
        if(children[index][i])
          children[index][i]->trackUnknown(track);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

void TreeElementComplex::getLinkedList(const QString &label, QList<TreeElement*> &list)
{
  try
  {
    TreeElement::getLinkedList(label, list);
    for(int index=0; index<children.size(); index++)
      for(int i=0; i<children[index].size(); i++)
        if(children[index][i])
          children[index][i]->getLinkedList(label, list);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

void TreeElementComplex::getLoopList(const QString &loop, QList<TreeElement*> &list)
{
  try
  {
    TreeElement::getLoopList(loop, list);
    for(int index=0; index<children.size(); index++)
      for(int i=0; i<children[index].size(); i++)
        if(children[index][i])
          children[index][i]->getLoopList(loop, list);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

void TreeElementComplex::getConditionList(const QString &condition, QList<TreeElement*> &list)
{
  try
  {
    TreeElement::getConditionList(condition, list);
    for(int index=0; index<children.size(); index++)
      for(int i=0; i<children[index].size(); i++)
        if(children[index][i])
          children[index][i]->getConditionList(condition, list);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

void TreeElementComplex::updateExpression()
{
  try
  {
    TreeElement::updateExpression();
    for(int index=0; index<children.size(); index++)
      for(int i=0; i<children[index].size(); i++)
        if(children[index][i])
          children[index][i]->updateExpression();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/
/***********************************************/

class TreeElementComplex::UndoCommandOverwrite : public TreeElement::UndoCommand
{
  TreeElementComplex *parent;
  int                 index;

public:
  UndoCommandOverwrite(TreeElement *treeElement, TreeElementComplex *_parent, int _index)
    : UndoCommand(treeElement, "overwrite"), parent(_parent), index(_index) {}
 ~UndoCommandOverwrite();

  void redo();
  void undo() {redo();}
};

/***********************************************/

TreeElementComplex::UndoCommandOverwrite::~UndoCommandOverwrite()
{
  delete treeElement;
}

/***********************************************/

void TreeElementComplex::UndoCommandOverwrite::redo()
{
  try
  {
    int indexChoice = parent->selectedIndex();

    TreeElement *oldElement = parent->children[indexChoice][index];
    parent->children[indexChoice][index] = treeElement;

    // copy the pointer to the add element
    treeElement->setElementAdd(oldElement->elementAdd());
    oldElement->setElementAdd(nullptr);

    // set autocomment for the first element
    if(parent->recieveAutoComments && (index==0))
      treeElement->setPushAutoComments(true);

    // new link content
    if((!treeElement->label().isEmpty()) && tree->rootElement())
      tree->rootElement()->newLink(treeElement);

    // create item
    if(oldElement->item())
    {
      treeElement->createItem(parent->item(), oldElement->item());
      oldElement->removeItem();
      tree->setSelectedItem(treeElement->item());
    }

    // update varList
    if(!treeElement->label().isEmpty() && treeElement->isLinked())
      tree->varList[treeElement->label()] = "{"+treeElement->selectedValue()+"}";

    treeElement = oldElement;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

Bool TreeElementComplex::overwrite(TreeElement *element, const QString &type, XmlNodePtr xmlNode)
{
  try
  {
    if((element->type() != type) || (!xmlNode))
      return false;

    // find element in the list
    int indexChoice = selectedIndex();
    int index = 0;
    if(indexChoice >= children.size())
      return false;
    for(; index<children[indexChoice].size(); index++)
      if(children[indexChoice][index] == element)
        break;
    if(index >= children[indexChoice].size())
      throw(Exception("element not in the list"));

    // Create new element
    xmlNode->setName(element->name());
    xmlNode->getAttribute("label"); // remove label
    if(!element->label().isEmpty())
      writeAttribute(xmlNode, "label", element->label());
    xmlNode->setName(element->name());
    TreeElement *newElement = TreeElement::newTreeElement(tree, this, element->xsdElement, "", xmlNode, true);
    if(tree->elementGlobal())
      tree->elementGlobal()->informAboutGlobalElements(newElement); // inform the new elements about all links

    tree->pushUndoCommand(new UndoCommandOverwrite(newElement, this, index));
    return true;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

TreeElementComplex::UndoCommandRemoveAddChild::UndoCommandRemoveAddChild(TreeElement *_newElement, TreeElement *_beforeElement, TreeElementComplex *_parent, Bool _isAdd)
     : UndoCommand(_newElement, (_isAdd ? "add" : "remove")),
       isAdd(_isAdd),
       isRemoved(_isAdd ? true : false),
       parent(_parent),
       beforeElement(_beforeElement),
       newElement(_newElement)
{
}

/***********************************************/

TreeElementComplex::UndoCommandRemoveAddChild::~UndoCommandRemoveAddChild()
{
  if(isRemoved)
    delete newElement;
}

/***********************************************/

void TreeElementComplex::UndoCommandRemoveAddChild::addChild()
{
  try
  {
    int indexChoice = parent->selectedIndex();

    TreeElementAdd *elementAdd = beforeElement ? beforeElement->elementAdd() : new TreeElementAdd(tree, tree->elementGlobal(), newElement->xsdElement, false/*visible*/);
    newElement->setElementAdd(elementAdd); // insert into the unbounded list
    if(tree->elementGlobal())
      tree->elementGlobal()->informAboutGlobalElements(newElement); // inform the new element about all links

    // disable auto comment of the first element (temporary)
    // -----------------------------------------------------
    if(parent->hasChildren(indexChoice))
      parent->children[indexChoice][0]->setPushAutoComments(false);

    // find element in the list and insert
    // -----------------------------------
    int idx=0;
    for(; idx<parent->children[indexChoice].size(); idx++)
      if(parent->children[indexChoice][idx] == beforeElement)
        break;
    if(beforeElement && idx >= parent->children[indexChoice].size())
      throw(Exception("element not in the list"));
    if(parent && parent == tree->elementGlobal() && idx == parent->children[indexChoice].size())
      idx--;
    if(!newElement->label().isEmpty() && parent != tree->elementGlobal()) // global element is added to non-global element
    {
      newElement->setLabel("");
      newElement->_name = newElement->xsdElement->name;
    }
    parent->children[indexChoice].insert(idx, newElement);
    isRemoved = false;

    // enable auto comment of the first element
    // ----------------------------------------
    if(parent->recieveAutoComments && parent->hasChildren(indexChoice))
      parent->children[indexChoice][0]->setPushAutoComments(true);

    // update varList
    // -------------
    if(tree->rootElement() && !newElement->label().isEmpty())
    {
      if(newElement->isLinked())
        tree->varList[newElement->label()] = "{"+newElement->selectedValue()+"}";
      else
        tree->varList[newElement->label()] = newElement->selectedValue();
      tree->rootElement()->newLink(treeElement);
    }

    // set visible and get the focus
    // -----------------------------
    if(parent->item())
    {
      // find precursor item
      TreeItem *after = nullptr;
      for(int i=0; i<parent->children[indexChoice].size(); i++)
      {
        if(parent->children[indexChoice][i] == newElement)
          break;
        if(parent->children[indexChoice][i]->item())
          after = parent->children[indexChoice][i]->item();
      }
      TreeItem *newItem = newElement->createItem(parent->item(), after);
      tree->setSelectedItem(newItem);
    }

    // if global element is added, update expressions of all elements in the tree
    if(!newElement->label().isEmpty() && tree->rootElement())
      tree->rootElement()->updateExpression();

    newElement->trackUnknown(true);
    newElement->trackRenamed(true);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

void TreeElementComplex::UndoCommandRemoveAddChild::removeChild()
{
  try
  {
    int indexChoice = parent->selectedIndex();

    // find element in the list & remove
    int index = 0;
    for(; index<parent->children[indexChoice].size(); index++)
      if(parent->children[indexChoice][index] == newElement)
        break;
    if(index >= parent->children[indexChoice].size())
      throw(Exception("element not in the list"));

    parent->children[indexChoice].remove(index);
    isRemoved = true;
    newElement->removeItem();
    newElement->setElementAdd(nullptr);

    // update varList
    // --------------
    if(!newElement->label().isEmpty())
    {
      // only delete varList entry if element was actually removed and not e.g. moved via drag&drop (= add copy, then remove)
      auto iter = std::find_if(parent->children[indexChoice].begin(), parent->children[indexChoice].end(), [&](TreeElement* el) { return el->label() == newElement->label(); });
      if(iter == parent->children[indexChoice].end())
        tree->varList.remove(newElement->label());
    }

    if(parent->recieveAutoComments && parent->hasChildren(indexChoice))
      parent->children[indexChoice][0]->setPushAutoComments(true);

    if(beforeElement && beforeElement->item())
      tree->setSelectedItem(beforeElement->item());

    // if global element is removed, update expressions of all elements in the tree
    if(!newElement->label().isEmpty() && tree->rootElement())
      tree->rootElement()->updateExpression();

    newElement->trackUnknown(false);
    newElement->trackRenamed(false);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

void TreeElementComplex::UndoCommandMoveChild::redo()
{
  try
  {
    const int indexChoice   = parent->selectedIndex();
    const int childrenCount = parent->children[indexChoice].size();

    auto findChildIndex = [&] (const TreeElement* el)
    {
      for(int i = 0; i < childrenCount; i++)
        if(parent->children[indexChoice][i] == el)
          return i;
      return -1;
    };

    int idxSource = findChildIndex(element);
    int idxTarget = findChildIndex(beforeElement);

    if(idxSource == idxTarget)
      throw(Exception("cannot move element in place"));

    if(idxSource < 0 || idxTarget < 0 || idxSource >= childrenCount || idxTarget >= childrenCount)
      throw(Exception("element not in the list"));

    // move treeElement
    parent->children[indexChoice].remove(idxSource);
    parent->children[indexChoice].insert(idxTarget, element);

    // move treeItem
    if(element->item() && element->parentElement && element->parentElement->item())
    {
      bool isExpanded = element->item()->isExpanded();
      QTreeWidgetItem* child = element->parentElement->item()->takeChild(idxSource);
      element->parentElement->item()->insertChild(idxTarget, child);
      tree->setSelectedItem(element->item());
      element->item()->setExpanded(isExpanded);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

Bool TreeElementComplex::canAddChild(TreeElement *beforeElement, const QString &type) const
{
  try
  {
    if(!beforeElement)
      throw(Exception("try to add before null element"));
    if(beforeElement->parentElement != this)
      throw(Exception("try to add by an external child"));
    if(!hasChildren(selectedIndex()))
      return false;
    if(type.isEmpty() || (type != beforeElement->type()))
      return false;
    // is element in an unbounded list?
    if(!beforeElement->elementAdd())
      return false;
    return true;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

Bool TreeElementComplex::canRemoveChild(TreeElement *element) const
{
  try
  {
    if(!element)
      return false;
    if(element->parentElement != this)
      throw(Exception("try to remove an external child"));
    if(element->parentElement == tree->elementGlobal())
      return element->parentElement->canRemoveChild(element);
    if(!hasChildren(selectedIndex()))
      return false;
    if(element->isUnknown()) // unknown elements can always be removed
      return true;
    if(element->isElementAdd()) // the add element cannot be removed
      return false;
    if(element->elementAdd()) // unbounded list?
      return (element->optional() || (element->elementAdd()->unboundedCount>1));
    return false;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

Bool TreeElementComplex::canMoveChild(TreeElement *element) const
{
  try
  {
    if(!element)
      return false;
    if(element->parentElement != this)
      throw(Exception("try to remove an external child"));
    if(!hasChildren(selectedIndex()))
      return false;
    if(element->isElementAdd()) // the add element cannot be removed
      return false;
    if(element->parentElement == tree->elementGlobal())
      return true;
    if(element->elementAdd()) // unbounded list?
      return (element->optional() || (element->elementAdd()->unboundedCount>1));
    return false;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}


/***********************************************/

Bool TreeElementComplex::addChild(TreeElement *beforeElement, const QString &type, XmlNodePtr xmlNode, bool /*moved*/)
{
  try
  {
    if(!TreeElementComplex::canAddChild(beforeElement, type))
      return false;

    // create & init new element
    TreeElement *newElement = TreeElement::newTreeElement(tree, this, beforeElement->xsdElement, beforeElement->defaultValue(), xmlNode, xmlNode!=nullptr);

    tree->pushUndoCommand(new UndoCommandRemoveAddChild(newElement, beforeElement, this, true/*isAdd*/));
    return true;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

Bool TreeElementComplex::removeChild(TreeElement *element)
{
  try
  {
    if(!TreeElementComplex::canRemoveChild(element))
      return false;

    // find element in the list
    int index = 0;
    for(; index<children[selectedIndex()].size(); index++)
      if(children[selectedIndex()][index] == element)
        break;

    if(index >= children[selectedIndex()].size())
      throw(Exception("element not in the list"));
    TreeElement *beforeElement = index+1 < children[selectedIndex()].size() ? children[selectedIndex()][index+1] : nullptr;

    tree->pushUndoCommand(new UndoCommandRemoveAddChild(element, beforeElement, this, false/*isAdd*/));
    return true;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

Bool TreeElementComplex::moveChild(TreeElement *element, TreeElement *beforeElement)
{
  try
  {
    if(!TreeElementComplex::canMoveChild(element))
      return false;

    if(!isElementGlobal() && element->elementAdd() != beforeElement->elementAdd())
      return false;

    // find element in the list
    int indexSource = -1;
    int indexTarget = -1;
    for(int i = 0; i<children[selectedIndex()].size(); i++)
    {
      if(children[selectedIndex()][i] == element)       indexSource = i;
      if(children[selectedIndex()][i] == beforeElement) indexTarget = i;
    }

    if(indexSource == indexTarget || indexSource < 0 || indexTarget < 0)
      return false;

    if(indexSource >= children[selectedIndex()].size() || indexTarget >= children[selectedIndex()].size())
      throw(Exception("element not in the list"));

    tree->pushUndoCommand(new UndoCommandMoveChild(element, beforeElement, this));
    return true;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/
/***********************************************/

void TreeElementComplex::getProgramList(QList<TreeElement*> &list)
{
  try
  {
    if(hasChildren(selectedIndex()))
      for(int i=0; i<children[selectedIndex()].size(); i++)
        children[selectedIndex()][i]->getProgramList(list);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

int TreeElementComplex::childrenCount() const
{
  try
  {
    if(hasChildren(selectedIndex()))
      return children[selectedIndex()].size();

    return 0;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

TreeElement* TreeElementComplex::childAt(int index) const
{
  try
  {
    if(hasChildren(selectedIndex()))
      if(index >= 0 && index < children[selectedIndex()].size())
        return children[selectedIndex()].at(index);

    return nullptr;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/
/***********************************************/

void TreeElementComplex::createChildrenItems()
{
  try
  {
    TreeItem *after = nullptr;
    if(item() && hasChildren(selectedIndex()))
      for(int i=0; i<children[selectedIndex()].size(); i++)
      {
        TreeItem *child = children[selectedIndex()][i]->createItem(item(), after);
        if(child)
          after = child;
      }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

void TreeElementComplex::removeChildrenItems()
{
  try
  {
    if(item())
      item()->setExpanded(false); // hide all children

    for(int index=0; index<children.size(); index++)
      for(int i=0; i<children[index].size(); i++)
        children[index][i]->removeItem();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

TreeItem *TreeElementComplex::createItem(TreeItem *parent, TreeItem *after)
{
  try
  {
    TreeItem *item = TreeElement::createItem(parent, after);
    createChildrenItems();
    return item;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

void TreeElementComplex::removeItem()
{
  try
  {
    removeChildrenItems();
    TreeElement::removeItem();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/
/***********************************************/
