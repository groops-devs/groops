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
#include "tree/treeElementComment.h"
#include "tree/treeElementUnknown.h"
#include "tree/treeElementComplex.h"

/***********************************************/

TreeElementComplex::TreeElementComplex(Tree *tree, TreeElementComplex *parentElement,
                                       XsdElementPtr xsdElement, const QString &defaultOverride,
                                       XmlNodePtr xmlNode, bool recieveAutoComments)
                   : TreeElement(tree, parentElement, xsdElement, defaultOverride, xmlNode),
                     recieveAutoComments(recieveAutoComments)
{
  // convert default string to QJsonObject
  if(!defaultValue.isEmpty())
  {
    QJsonDocument doc = QJsonDocument::fromJson(defaultValue.toUtf8());
    if(!doc.isNull() && doc.isObject())
      defaultObject = doc.object();
    else
      defaultObject = QJsonDocument::fromJson(QString("{\""+defaultValue+"\" : {}}").toUtf8()).object();
  }
}

/***********************************************/

TreeElementComplex::~TreeElementComplex()
{
  for(auto &childdrenAtIndex : children_)
    for(auto &child : childdrenAtIndex)
      delete child;
}

/***********************************************/

void TreeElementComplex::createXmlChildren(XmlNodePtr xmlNode, bool createRootEvenIfEmpty) const
{
  try
  {
    XmlNodePtr xmlChild;
    if(xmlNode && !isLinked())
      for(auto &child : children_[selectedIndex()])
        if((xmlChild = child->createXmlTree(createRootEvenIfEmpty)))
          xmlNode->addChild(xmlChild);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/
/***********************************************/

int TreeElementComplex::addChoice(const QString &value, XsdElementPtr xsdComplex, const QJsonObject &defaultObject)
{
  try
  {
    int index = insertNewValue(value, false);
    if(index != xsdElementList.size())
      throw(Exception("lists not consistent"));
    xsdElementList.push_back(xsdComplex);
    overrideDefaultObjects.push_back(defaultObject);
    children_.push_back(QVector<TreeElement*>());
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
    if((index < 0) || (index >= xsdElementList.size()))
      return;

    // children not created yet?
    if(xmlNode && children_[index].size())
      throw(Exception("cannot create"));
    if(!children_[index].size())
    {
      if(xsdElementList[index] && xsdElementList[index]->complex)
      {
        auto xsdComplex = xsdElementList[index]->complex;
        if(xsdComplex->type != "sequence")
          throw(Exception("complex must be sequence"));

        for(auto &xsdElement : xsdComplex->elements)
        {
          QJsonValue overrideDefaultValue = overrideDefaultObjects[index].value(xsdElement->names.front());
          if((overrideDefaultValue.isUndefined() || overrideDefaultValue.isNull()) && xsdElement->unbounded && xsdElement->defaultValue.startsWith("["))
            overrideDefaultValue = QJsonDocument::fromJson(xsdElement->defaultValue.toUtf8()).array(); // catch arrays in default

          // lambda function to add new tree element as child
          // -------------------------------------------------
          auto addNewChildElement = [&](TreeElementAdd *elementAdd, const QJsonValue &defaultValue, bool fillWithDefaults)
          {
            XmlNodePtr xmlChild;
            if(xmlNode)
            {
              // insert possible comment elements
              while(xmlNode->hasChildren() && (xmlNode->peekNextChild()->getName() == "COMMENT"))
                children_[index].push_back(new TreeElementComment(tree, this, xmlNode->getNextChild()->getText()));
              xmlChild = xmlNode->getChild(xsdElement->names);
            }

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
              treeElement = TreeElement::newTreeElement(tree, this, xsdElement, defaultStr, xmlChild, fillWithDefaults);
            }
            catch(std::exception &/*e*/)
            {
              // add empty element from schema and add xmlNode as an unknown element afterwards
              treeElement = TreeElement::newTreeElement(tree, this, xsdElement, "", XmlNodePtr(nullptr), fillWithDefaults);
              xmlNode->addChild(xmlChild);
            }
            if(!treeElement)
              throw(Exception("cannot add child element"));
            treeElement->setElementAdd(elementAdd);
            children_[index].push_back(treeElement);
          };
          // -------------------------------------------------

          if(!xsdElement->unbounded)
            addNewChildElement(nullptr, overrideDefaultValue, !xmlNode);
          else
          {
            // the add element
            TreeElementAdd *elementAdd = new TreeElementAdd(tree, this, xsdElement);
            if(xmlNode)
            {
              UInt childCount = xmlNode->getChildCount(xsdElement->names);
              for(UInt i=0; i<childCount; i++)
                addNewChildElement(elementAdd, overrideDefaultValue, false);
            }
            else if(overrideDefaultValue.isArray()) // create multiple elements from default array
            {
              for(const auto &arrayElement : overrideDefaultValue.toArray())
                addNewChildElement(elementAdd, arrayElement, true);
            }
            else if(!xsdElement->optional) // not optional -> create at least one element
              addNewChildElement(elementAdd, overrideDefaultValue, true);

            children_[index].push_back(elementAdd);
          } // if(unbounded)
        } // for(xsdElements)
      } // if(xsdComplex)

      // insert possible comment elements
      // --------------------------------
      while(xmlNode && xmlNode->hasChildren() && (xmlNode->peekNextChild()->getName() == "COMMENT"))
        children_[index].push_back(new TreeElementComment(tree, this, xmlNode->getNextChild()->getText()));

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
          children_[index].push_back(treeElement);
        }
        children_[index].push_back(elementAdd);
      }

      // inform the new elements about all links
      // ---------------------------------------
      if(tree->elementGlobal && tree->rootElement)
      {
        for(auto &child : children_[index])
          tree->elementGlobal->informAboutGlobalElements(child, false/*recursively*/);
        for(auto &child : children_[index])
          child->updateParserResults(tree->elementGlobal->variableList(), false/*recursively*/);
      }
    } // if(!hasChildren)

    childSetPushAutoComments(true, index); // set autocomment for the first element
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

void TreeElementComplex::setSelectedIndex(int index)
{
  try
  {
    TreeElement::setSelectedIndex(index);

    removeChildrenItems();
    if(isLinked())
      return;

    createChildrenElements(index, XmlNodePtr(nullptr));

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

void TreeElementComplex::informAboutLink(TreeElement *elementInGlobal, bool recursively)
{
  try
  {
    TreeElement::informAboutLink(elementInGlobal, recursively);
    if(recursively)
      for(auto &childdrenAtIndex : children_)
        for(auto &child : childdrenAtIndex)
          child->informAboutLink(elementInGlobal, recursively);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

void TreeElementComplex::updateParserResults(const VariableList &varList, bool recursively)
{
  try
  {
    if(recursively)
      for(auto &childdrenAtIndex : children_)
        for(auto &child : childdrenAtIndex)
          child->updateParserResults(varList, recursively);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

void TreeElementComplex::addedLink(TreeElement *elementInGlobal)
{
  try
  {
    TreeElement::addedLink(elementInGlobal);
    for(auto &childdrenAtIndex : children_)
      for(auto &child : childdrenAtIndex)
        child->addedLink(elementInGlobal);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

void TreeElementComplex::removedLink(TreeElement *elementInGlobal)
{
  try
  {
    TreeElement::removedLink(elementInGlobal);
    for(auto &childdrenAtIndex : children_)
      for(auto &child : childdrenAtIndex)
        child->removedLink(elementInGlobal);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

void TreeElementComplex::renamedLink(const QString &oldLabel, const QString &newLabel)
{
  try
  {
    TreeElement::renamedLink(oldLabel, newLabel);
    for(auto &childdrenAtIndex : children_)
      for(auto &child : childdrenAtIndex)
        child->renamedLink(oldLabel, newLabel);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/
/***********************************************/

class TreeElementComplex::UndoCommandOverwriteChildren : public TreeElement::UndoCommand
{
  TreeElementComplex *treeElement;
  XmlNodePtr          xmlNode;

public:
  UndoCommandOverwriteChildren(TreeElementComplex *treeElement, XmlNodePtr xmlNode)
    : UndoCommand(treeElement, "overwrite"), treeElement(treeElement), xmlNode(xmlNode) {}
 ~UndoCommandOverwriteChildren() {}

  void redo();
  void undo() {redo();}
};

/***********************************************/

void TreeElementComplex::UndoCommandOverwriteChildren::redo()
{
  try
  {
    XmlNodePtr xmlNodeOld = XmlNode::create("parent");
    treeElement->createXmlChildren(xmlNodeOld);

    treeElement->removeChildrenItems();
    for(auto &child : treeElement->children_[treeElement->selectedIndex()])
      delete child;
    treeElement->children_[treeElement->selectedIndex()].clear();

    treeElement->createChildrenElements(treeElement->selectedIndex(), xmlNode);
    treeElement->createChildrenItems();
    tree->treeChanged();

    xmlNode = xmlNodeOld;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

void TreeElementComplex::overwriteChildren(XmlNodePtr xmlNode)
{
  tree->undoStack->push(new UndoCommandOverwriteChildren(this, xmlNode));
}

/***********************************************/
/***********************************************/

TreeElement *TreeElementComplex::skipCommentElements(TreeElement *targetElement, int index) const
{
  if(index < 0)
    index = selectedIndex();
  auto end  = children_[index].end();
  auto iter = std::find(children_[index].begin(), end, targetElement);
  while((iter != end) && ((*iter)->type() == "COMMENT")) // skip comment elements inbetween
    iter++;
  return (iter != end) ? *iter : nullptr;
}

/***********************************************/

void TreeElementComplex::childSetPushAutoComments(bool on, int index)
{
  if(index < 0)
    index = selectedIndex();
  if(!children_[index].size())
    return;
  TreeElement *child = skipCommentElements(children_[index][0], index);
  if(child)
    child->setPushAutoComments(on && recieveAutoComments);
}

/***********************************************/

void TreeElementComplex::UndoCommandRemoveAddChild::redo()
{
  try
  {
    // disable auto comment of the first element (temporary)
    parent->childSetPushAutoComments(false);

    auto &children = parent->children_[parent->selectedIndex()];
    if(toAdd) // add child
    {
      // find element in the list and insert
      int index = std::distance(children.begin(), std::find(children.begin(), children.end(), targetElement));
      children.insert(index, treeElement);
      treeElement->setElementAdd(addElement); // insert into the unbounded list

      // set visible and get the focus
      if(parent->item())
      {
        // find precursor item
        TreeItem *afterItem = nullptr;
        for(int i=0; (i<children.size()) && (children[i] != treeElement); i++)
          if(children[i]->item())
            afterItem = children[i]->item();
        tree->setSelectedItem(treeElement->createItem(parent->item(), afterItem));
      }
    }
    else // remove child
    {
      // find element in the list & remove
      int index = std::distance(children.begin(), std::find(children.begin(), children.end(), treeElement));
      children.remove(index);
      addElement = treeElement->elementAdd();
      treeElement->setElementAdd(nullptr);
      treeElement->removeItem();

      if(targetElement->item())
        tree->setSelectedItem(targetElement->item());
    }

    if(dynamic_cast<TreeElementGlobal*>(parent))
      dynamic_cast<TreeElementGlobal*>(parent)->updateVariableList();

    // enable auto comment of the first element
    parent->childSetPushAutoComments(true);

    tree->treeChanged();

    toAdd = !toAdd;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

bool TreeElementComplex::canAddChild(TreeElement *targetElement, const QString &type) const
{
  try
  {
    if(!targetElement || (targetElement->parentElement != this))
      throw(Exception("invalid call"));
    if(type.isEmpty())
      return false;
    if(type == "COMMENT")
      return true;
    targetElement = skipCommentElements(targetElement);
    return (targetElement->type() == type) && targetElement->unbounded();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

bool TreeElementComplex::canRemoveChild(TreeElement *element) const
{
  try
  {
    if(!element)
      return false;
    if(element->parentElement != this)
      throw(Exception("try to remove an external child"));
    return element->unbounded() && (element->optional() || (element->elementAdd() && element->elementAdd()->unboundedCount > 1));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

TreeElement *TreeElementComplex::addChild(TreeElement *targetElement, const QString &type, XmlNodePtr xmlNode)
{
  try
  {
    if(!canAddChild(targetElement, type))
      return nullptr;

    if(type == "COMMENT")
    {
      TreeElement *newElement = new TreeElementComment(tree, this, xmlNode ? xmlNode->getText() : QString());
      tree->undoStack->push(new UndoCommandRemoveAddChild(newElement, targetElement, this, true/*isAdd*/));
      return newElement;
    }

    // create & init new element
    TreeElement *targetElement2 = skipCommentElements(targetElement);
    TreeElement *newElement = TreeElement::newTreeElement(tree, this, targetElement2->xsdElement, targetElement2->defaultValue, xmlNode, false/*fillWithDefaults*/);
    newElement->_name  = newElement->_schemaName;
    newElement->_label = "";
    newElement->setElementAdd(targetElement2->elementAdd());
    tree->elementGlobal->informAboutGlobalElements(newElement, false/*recursively*/); // inform the new element about all links
    newElement->updateParserResults(tree->elementGlobal->variableList(), false/*recursively*/);

    tree->undoStack->push(new UndoCommandRemoveAddChild(newElement, targetElement, this, true/*isAdd*/));
    return newElement;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

bool TreeElementComplex::removeChild(TreeElement *element)
{
  try
  {
    if(!canRemoveChild(element))
      return false;

    // find element in the list
    auto &children = children_[selectedIndex()];
    int index = std::distance(children.begin(), std::find(children.begin(), children.end(), element));
    if(index+1 >= children_[selectedIndex()].size())
      throw(Exception("element not in the list or add element is missing"));

    tree->undoStack->push(new UndoCommandRemoveAddChild(element, children[index+1], this, false/*isAdd*/));
    return true;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/
/***********************************************/

class TreeElementComplex::UndoCommandMoveChild : public TreeElement::UndoCommand
{
  TreeElement *targetElement;

public:
  UndoCommandMoveChild(TreeElement *treeElement, TreeElement *targetElement)
  : UndoCommand(treeElement, "move"), targetElement(targetElement) {}
  ~UndoCommandMoveChild() {}

  void redo();
  void undo() {redo();}
};

/***********************************************/

void TreeElementComplex::UndoCommandMoveChild::redo()
{
  try
  {
    TreeElementComplex *parent = treeElement->parentElement;
    parent->childSetPushAutoComments(false); // disable auto comment of the first element (temporary)
    treeElement->removeItem();

    // find element in the list & remove
    auto &children = parent->children_[parent->selectedIndex()];
    int index = std::distance(children.begin(), std::find(children.begin(), children.end(), treeElement));
    TreeElement *targetElementOld = children.at(index+1);
    children.remove(index);
    index = std::distance(children.begin(), std::find(children.begin(), children.end(), targetElement));
    children.insert(index, treeElement);
    targetElement = targetElementOld;

    parent->childSetPushAutoComments(true); // enable auto comment of the first element

    // set visible and get the focus
    // -----------------------------
    if(parent->item())
    {
      // find precursor item
      TreeItem *after = nullptr;
      for(int i=0; i<children.size(); i++)
      {
        if(children[i] == treeElement)
          break;
        if(children[i]->item())
          after = children[i]->item();
      }
      tree->setSelectedItem(treeElement->createItem(parent->item(), after));
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

bool TreeElementComplex::canMoveChild(TreeElement *targetElement, TreeElement *element) const
{
  if((targetElement == element) || (element->parentElement != this) || (targetElement->parentElement != this))
    return false;
  if(dynamic_cast<TreeElementAdd*>(element))
    return false;
  if(dynamic_cast<TreeElementComment*>(element))
    return true;
  return element->elementAdd() && (element->elementAdd() == skipCommentElements(targetElement)->elementAdd());
}

/***********************************************/

void TreeElementComplex::moveChild(TreeElement *targetElement, TreeElement *element)
{
  try
  {
    if(canMoveChild(targetElement, element))
      tree->undoStack->push(new UndoCommandMoveChild(element, targetElement));
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
    if(item() && !isLinked())
      for(auto &childElement : children_[selectedIndex()])
      {
        TreeItem *child = childElement->createItem(item(), after);
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

    for(auto &childdrenAtIndex : children_)
      for(auto &child : childdrenAtIndex)
        child->removeItem();
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
