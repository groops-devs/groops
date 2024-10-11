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
#include "tree/treeElementLoopCondition.h"
#include "tree/treeElementComment.h"
#include "tree/treeElementUnknown.h"
#include "tree/treeElementComplex.h"

/***********************************************/

TreeElementComplex::TreeElementComplex(Tree *tree, TreeElementComplex *parentElement,
                                       XsdElementPtr xsdElement, const QString &defaultOverride,
                                       XmlNodePtr xmlNode, bool recieveAutoComments)
                   : TreeElement(tree, parentElement, xsdElement, defaultOverride, xmlNode),
                     recieveAutoComments(recieveAutoComments),
                     initializedVariables(false)
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
        {
          // write _local loops and conditions before element
          XmlAttrPtr attr = xmlChild->findAttribute("loop");
          if(attr && (attr->text == "_localLoop_"))
            xmlNode->addChild(xmlChild->getChild("loopType"));

          attr = xmlChild->findAttribute("condition");
          if(attr && (attr->text == "_localCondition_"))
            xmlNode->addChild(xmlChild->getChild("conditionType"));

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

    // children already created?
    if(children_[index].size())
    {
      if(xmlNode)
        throw(Exception("cannot create"));
      childSetPushAutoComments(true, index); // set autocomment for the first element
      return;
    }

    auto defaultStr = [](const QJsonValue &defaultValue) -> QString
    {
      switch(defaultValue.type())
      {
        case QJsonValue::String: return defaultValue.toString();
        case QJsonValue::Bool:   return (defaultValue.toBool() ? "1" : "0");
        case QJsonValue::Double: return QString::number(defaultValue.toDouble());
        case QJsonValue::Object: if(!defaultValue.toObject().isEmpty()) return QJsonDocument(defaultValue.toObject()).toJson(QJsonDocument::Compact); break;
        case QJsonValue::Array:  if(!defaultValue.toArray().isEmpty())  return QJsonDocument(defaultValue.toArray()).toJson(QJsonDocument::Compact);  break;
        default: break;
      }
      return QString();
    };

    // lambda function to add comments & variables inbetween
    // -----------------------------------------------------
    auto insertCommentsAndVariables = [&]()
    {
      if(!xmlNode)
        return;
      XmlNodePtr xmlChild;
      XmlAttrPtr attrLabel;
      std::list<XmlNodePtr> xmlLocals;
      for(;;)
      {
        if(!(xmlChild = xmlNode->peekNextChild()))
          break;
        if(xmlChild->getName() == "COMMENT")
          children_[index].push_back(new TreeElementComment(tree, this, xmlNode->getNextChild()->getText()));
        else if((attrLabel = xmlChild->findAttribute("label")))
        {
          xmlChild = xmlNode->getNextChild(); // remove node
          if(attrLabel->text.startsWith("_local"))
            xmlLocals.push_front(xmlChild); // _local variables (loops/conditions) are added to next element
          else
          {
            if(xmlLocals.size())
              throw(Exception("variables cannot have loops or conditions"));
            children_[index].push_back(TreeElement::newTreeElement(tree, this, tree->xsdElement(xmlChild->getName()), "", xmlChild, false/*fillWithDefaults*/));
          }
        }
        else
          break;
      }
      // add _local variables to next xml element
      if(xmlChild && xmlLocals.size())
        for(auto &xmlLocal : xmlLocals)
          xmlChild->addChild(xmlLocal, true);
    };
    // -------------------------------------------------

    if(xsdElementList[index] && xsdElementList[index]->complex)
    {
      auto xsdComplex = xsdElementList[index]->complex;
      if(xsdComplex->type != "sequence")
        throw(Exception("complex must be sequence"));

      for(auto &xsdElement : xsdComplex->elements)
      {
        QJsonValue overrideDefaultValue = overrideDefaultObjects[index].value(xsdElement->names.front());

        // lambda function to add new tree element as child
        // -------------------------------------------------
        auto addNewChildElement = [&](TreeElementAdd *elementAdd, const QJsonValue &defaultValue)
        {
          insertCommentsAndVariables();
          XmlNodePtr xmlChild;
          if(xmlNode)
            xmlChild = xmlNode->getChild(xsdElement->names);

          TreeElement *treeElement;
          try
          {
            treeElement = TreeElement::newTreeElement(tree, this, xsdElement, defaultStr(defaultValue), xmlChild, !xmlNode/*fillWithDefaults*/);
          }
          catch(std::exception &/*e*/)
          {
            // add empty element from schema and add xmlNode as an unknown element afterwards
            treeElement = TreeElement::newTreeElement(tree, this, xsdElement, "", nullptr, !xmlNode/*fillWithDefaults*/);
            xmlNode->addChild(xmlChild);
          }
          treeElement->setElementAdd(elementAdd);
          children_[index].push_back(treeElement);
        };
        // -------------------------------------------------

        if(!xsdElement->unbounded)
          addNewChildElement(nullptr, overrideDefaultValue);
        else
        {
          if((overrideDefaultValue.isUndefined() || overrideDefaultValue.isNull()) && xsdElement->defaultValue.startsWith("["))
            overrideDefaultValue = QJsonDocument::fromJson(xsdElement->defaultValue.toUtf8()).array(); // catch arrays in default

          // the add element
          TreeElementAdd *elementAdd = new TreeElementAdd(tree, this, xsdElement, defaultStr(defaultValue));
          if(xmlNode)
          {
            UInt childCount = xmlNode->getChildCount(xsdElement->names);
            for(UInt i=0; i<childCount; i++)
              addNewChildElement(elementAdd, overrideDefaultValue);
          }
          else if(overrideDefaultValue.isArray()) // create multiple elements from default array
          {
            for(const auto &arrayElement : overrideDefaultValue.toArray())
              addNewChildElement(elementAdd, arrayElement);
          }
          else if(!xsdElement->optional) // not optional -> create at least one element
            addNewChildElement(elementAdd, overrideDefaultValue);

          children_[index].push_back(elementAdd);
        } // if(unbounded)
      } // for(xsdElements)
    } // if(xsdComplex)

    insertCommentsAndVariables();

    // are there unknown XML nodes?
    if(xmlNode && xmlNode->hasChildren())
      while(xmlNode->hasChildren())
        children_[index].push_back(new TreeElementUnknown(tree, this, xmlNode->getNextChild()));

    // invisible add element to be able to remove links and comments at end of list
    children_[index].push_back(new TreeElementAdd(tree, this, nullptr, "", false/*visible*/));

    updateParserResultsInScope();
    updateLinksInScope();
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

void TreeElementComplex::updateParserResultsInScope()
{
  try
  {
    if(!initializedVariables)
      return;
    auto varList = this->varList; // without added local variables
    for(auto &child : children_[selectedIndex()])
      varList = child->updateParserResults(varList, true);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

VariableListPtr TreeElementComplex::updateParserResults(VariableListPtr varList, Bool /*addVariableInReturn*/)
{
  try
  {
    this->varList = varList;
    initializedVariables = true;
    TreeElement::updateParserResults(varList, false);
    for(auto &childdrenAtIndex : children_)
    {
      auto varListLocal = this->varList; // without added local variables
      for(auto &child : childdrenAtIndex)
        varListLocal = child->updateParserResults(varListLocal, true);
    }
    return varList;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

void TreeElementComplex::updateLinksInScope()
{
  try
  {
    if(!initializedVariables)
      return;
    auto labelTypesLocal = this->labelTypes; // without added local variables
    for(auto &child : children_[selectedIndex()])
      child->updateLinks(labelTypesLocal);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

void TreeElementComplex::updateLinks(QMap<QString, QString> &labelTypes)
{
  try
  {
    this->labelTypes = labelTypes;
    initializedVariables = true;
    TreeElement::updateLinks(labelTypes);
    for(auto &childdrenAtIndex : children_)
    {
      auto labelTypesLocal = this->labelTypes; // without added local variables
      for(auto &child : childdrenAtIndex)
        child->updateLinks(labelTypesLocal);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

bool TreeElementComplex::renamedLink(const QString &oldLabel, const QString &newLabel)
{
  try
  {
    if(!TreeElement::renamedLink(oldLabel, newLabel))
      return false;
    labelTypes[newLabel] = labelTypes.value(oldLabel);
    if(oldLabel != newLabel)
      labelTypes.remove(oldLabel);
    for(auto &childdrenAtIndex : children_)
      for(auto &child : childdrenAtIndex)
        if(!child->renamedLink(oldLabel, newLabel))
          break;
    return true;
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
  while((iter != end) && (((*iter)->type() == "COMMENT") || !(*iter)->label().isEmpty())) // skip comment elements inbetween
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
        if(parent->loop      && parent->loop->item())      afterItem = parent->loop->item();
        if(parent->condition && parent->condition->item()) afterItem = parent->condition->item();
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

    parent->updateParserResultsInScope();
    parent->updateLinksInScope();

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

bool TreeElementComplex::canAddChild(TreeElement *targetElement, const QString &type, const QString &label) const
{
  try
  {
    if(!targetElement || (targetElement->parentElement != this))
      throw(Exception("invalid call"));
    if(type.isEmpty())
      return false;
    if(type == "COMMENT")
      return true;
    if(!label.isEmpty())
      return (tree->xsdElement(type) != nullptr);
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
    if(dynamic_cast<TreeElementAdd*>(element))
      return false;
    return element->unbounded() && (element->optional() || (element->elementAdd() && element->elementAdd()->unboundedCount > 1));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

TreeElement *TreeElementComplex::addChild(TreeElement *targetElement, const QString &type, const QString &label, XmlNodePtr xmlNode)
{
  try
  {
    if(!canAddChild(targetElement, type, label))
      return nullptr;

    TreeElement *newElement;
    if(type == "COMMENT")
    {
      newElement = new TreeElementComment(tree, this, xmlNode ? xmlNode->getText() : QString());
    }
    else if(!label.isEmpty())
    {
      newElement = TreeElement::newTreeElement(tree, this, tree->xsdElement(type), "", xmlNode, !xmlNode/*fillWithDefaults*/);
      newElement->_name  = newElement->_schemaName;
      newElement->_label = label;
    }
    else
    {
      TreeElement *targetElement2 = skipCommentElements(targetElement);
      newElement = TreeElement::newTreeElement(tree, this, targetElement2->xsdElement, targetElement2->defaultValue, xmlNode, !xmlNode/*fillWithDefaults*/);
      newElement->_name  = newElement->_schemaName;
      newElement->_label = "";
      newElement->setElementAdd(targetElement2->elementAdd());
    }

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

    parent->updateParserResultsInScope();
    parent->updateLinksInScope();
    parent->childSetPushAutoComments(true); // enable auto comment of the first element
    tree->treeChanged();

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
    TreeItem *after = (item()->childCount()) ? dynamic_cast<TreeItem*>(item()->child(item()->childCount()-1)) : nullptr;
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
