/***********************************************/
/**
* @file treeElementGlobal.cpp
*
* @brief The global element.
*
* @author Torsten Mayer-Guerr
* @author Sebastian Strasser
* @date 2006-07-10
*/
/***********************************************/

#include <QtDebug>
#include <QInputDialog>
#include "base/importGroops.h"
#include "tree/tree.h"
#include "tree/treeElement.h"
#include "tree/treeElementAdd.h"
#include "tree/treeElementComplex.h"
#include "tree/treeElementUnknown.h"
#include "tree/treeElementGlobal.h"

/***********************************************/

TreeElementGlobal::TreeElementGlobal(Tree *tree, TreeElementComplex *parentElement,
                                     XsdElementPtr xsdElement, XmlNodePtr xmlNode)
  : TreeElementComplex(tree, parentElement, xsdElement, "", xmlNode, false, true)
{
  try
  {
    tree->_elementGlobal = this;
    addChoice("", xsdElement, QJsonObject());
    createChildrenElements(0, xmlNode);
    if(isLinked())
      throw(Exception("global cannot be linked"));
    setSelectedIndex(0);

    // the add element
    TreeElementAdd *elementAdd = new TreeElementAdd(tree, this, xsdElement, true/*visible*/);
    for(auto &&child : children[0])
      child->setElementAdd(elementAdd);
    children[0].push_back(elementAdd);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

TreeElementGlobal::~TreeElementGlobal()
{
  try
  {
    tree->_elementGlobal = nullptr;
  }
  catch(std::exception &e)
  {
    qDebug() << QString::fromStdString("Exception in destructor at "+_GROOPS_ERRORLINE+"\n"+e.what());
  }
}

/***********************************************/

XmlNodePtr TreeElementGlobal::getXML(Bool /*withEmptyNodes*/) const
{
  try
  {
    XmlNodePtr xmlNode = getBaseXML();
    if((xmlNode==nullptr) || isLinked())
      return xmlNode;
    getChildrenXML(xmlNode, true/*withEmptyNodes*/);
    return xmlNode;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

void TreeElementGlobal::informAboutGlobalElements(TreeElement *element) const
{
  try
  {
    for(int i=0; i<children[0].size(); i++)
      if(!children[0][i]->isElementAdd()) // Do not send the Add element
        element->newLink(children[0][i]);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

XsdComplexPtr TreeElementGlobal::xsdComplexGlobalTypes() const
{
  try
  {
    return xsdElement->complex;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

QStringList TreeElementGlobal::getChildrenNames() const
{
  try
  {
    QStringList childrenNames;
    for(int i = 0; i < children.size(); i++)
      for(int j = 0; j < children.at(i).size(); j++)
        if(children.at(i).at(j)->canDisabled())
          childrenNames.push_back(children.at(i).at(j)->name());

    return childrenNames;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

void TreeElementGlobal::createChildrenElements(int index, XmlNodePtr xmlNode)
{
  try
  {
    if((index<0) || (index>=xsdElementList.size()))
      return;

    if(xmlNode && hasChildren(index))
      throw(Exception("cannot create"));

    // children not created yet?
    if(!hasChildren(index) && xmlNode)
    {
      XsdComplexPtr xsdComplex = xsdComplexGlobalTypes();
      XmlNodePtr xmlChild = xmlNode->getNextChild();
      while(xmlChild)
      {
        XsdElementPtr xsdElement = xsdComplex->getXsdElement(xmlChild->getName());

        if(!xsdElement)
        {
          TreeElementUnknown *treeElement = new TreeElementUnknown(tree, this, xmlChild);
          treeElement->setElementAdd(elementAdd());
          children[index].push_back(treeElement);
          treeElement->trackUnknown(true);
          xmlChild = xmlNode->getNextChild();
          continue;
        }

        TreeElement* treeElement;
        try
        {
          treeElement = TreeElement::newTreeElement(tree, this, xsdElement, "", xmlChild, true);
        }
        catch(...)
        {
          xmlChild = xmlNode->getNextChild();
          continue;
        }

        if(treeElement && xmlChild && treeElement->type() != xmlChild->getName())
          treeElement->setOriginalName(xmlChild->getName());
        if(treeElement)
          treeElement->trackRenamed(true);

        // update varList
        if(treeElement->isLinked())
          tree->setVarList().addVariable(ExpressionVariablePtr(new ExpressionVariable(treeElement->label().toStdString(), "{"+treeElement->selectedValue().toStdString()+"}")));
        else
          tree->setVarList().addVariable(ExpressionVariablePtr(new ExpressionVariable(treeElement->label().toStdString(), treeElement->selectedValue().toStdString())));

        children[index].push_back(treeElement);
        xmlChild = xmlNode->getNextChild();
      }
    } // if(!hasChildren)

    std::stable_sort(children[index].begin(), children[index].end(), [](TreeElement* e1, TreeElement* e2){ return (!e1->isUnknown() && e2->isUnknown()); });

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

Bool TreeElementGlobal::canAddChild(TreeElement *beforeElement, const QString &type) const
{
  if(isProgram())
    return false;

  if(beforeElement && beforeElement->isElementAdd())
    return true;

  std::vector<XsdElementPtr> schemaTypes = xsdComplexGlobalTypes()->element;
  for(UInt i = 0; i < schemaTypes.size(); i++)
    if(schemaTypes.at(i)->type == type)
      return true;

  return false;
}

/***********************************************/

Bool TreeElementGlobal::canRemoveChild(TreeElement *element) const
{
  if(element && element->isElementAdd())
    return false;

  return true;
}


/***********************************************/

Bool TreeElementGlobal::canSetGlobal(TreeElement *element) const
{
  try
  {
    if((!element) || element->isElementAdd())
      return false;

    std::vector<XsdElementPtr> xsdElementList = xsdComplexGlobalTypes()->element;
    for(UInt i=0; i<xsdElementList.size(); i++)
      if(xsdElementList.at(i)->type == element->type())
        return true;
    return false;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

Bool TreeElementGlobal::setGlobal(TreeElement *element)
{
  try
  {
    if(!canSetGlobal(element))
      return false;

    // new global elements are added to the bottom of the list by default
    TreeElement *beforeElement = children[0].back();

    // if element is already in global, new global element is added directly in front of it
    if(element->parentElement == this)
      beforeElement = element;

    // add element to global section
    tree->undoStack()->beginMacro("set global "+element->name());
    QString label = element->name();
    if(!addChild(beforeElement, element->type(), element->getXML(), label))
    {
      // abort/end macro and overwrite it with empty command
      tree->undoStack()->endMacro();
      tree->undoStack()->undo();
      tree->undoStack()->push(new QUndoCommand);
      tree->undoStack()->undo();
      return false;
    }

    // element becomes linked
    int index = element->findLinkIndex(label);
    if(index < 0)
      throw(Exception("cannot find link"));
    element->changeSelectedIndex(index);
    element->updateExpression();

    tree->undoStack()->endMacro();

    return true;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

Bool TreeElementGlobal::addChild(TreeElement *beforeElement, const QString &type, XmlNodePtr xmlNode, bool moved)
{
  QString label;
  return addChild(beforeElement, type, xmlNode, label, moved);
}

/***********************************************/

Bool TreeElementGlobal::addChild(TreeElement *beforeElement, const QString &type, XmlNodePtr xmlNode, QString &label, bool moved)
{
  try
  {
    if(!TreeElementGlobal::canAddChild(beforeElement, type))
      throw(Exception("Cannot add child to global section"));

    QStringList existingNames = getChildrenNames();
    if(!moved || existingNames.contains(label))
    {
      if(xmlNode && label.isEmpty())
      {
        XmlAttrPtr attributeLabel = xmlNode->getAttribute("label");
        label = attributeLabel ? attributeLabel->getText() : xmlNode->getName();
      }

      bool ok;
      label = QInputDialog::getText(tree, tr("Add global element - GROOPS"), tr("Name of global element:"), QLineEdit::Normal, label, &ok);
      QRegExp regex("[a-zA-Z]([a-zA-Z0-9])*");
      while(ok && (label.isEmpty() || existingNames.contains(label) || !regex.exactMatch(label)))
        label = QInputDialog::getText(tree, tr("Add global element - GROOPS"), tr("Name already exists or is invalid (only letters and digits allowed)!\nChoose another name:"), QLineEdit::Normal, label, &ok);

      if(!ok)
        return false;
    }

    // set schema type of new element
    XsdElementPtr xsdElement;
    std::vector<XsdElementPtr> schemaTypes = xsdComplexGlobalTypes()->element;
    for(UInt i = 0; i < schemaTypes.size(); i++)
      if(schemaTypes.at(i)->type == type)
        xsdElement = schemaTypes.at(i);

    // remove loop and condition attributes
    if(xmlNode)
    {
      xmlNode->getAttribute("loop");
      xmlNode->getAttribute("condition");
    }

    // create & init new element
    TreeElement *newElement = TreeElement::newTreeElement(tree, this, xsdElement, "", xmlNode, xmlNode!=nullptr);

    if(!label.isEmpty())
      newElement->setLabel(label);

    tree->pushUndoCommand(new UndoCommandRemoveAddChild(newElement, beforeElement, this, true/*isAdd*/));
    if(newElement->isLinked())
      tree->setVarList().addVariable(ExpressionVariablePtr(new ExpressionVariable(newElement->label().toStdString(), "{"+newElement->selectedValue().toStdString()+"}")));
    else
      tree->setVarList().addVariable(ExpressionVariablePtr(new ExpressionVariable(newElement->label().toStdString(), newElement->selectedValue().toStdString())));

    return true;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

Bool TreeElementGlobal::removeChild(TreeElement *element)
{
  try
  {
    tree->undoStack()->beginMacro("remove link "+element->name());

    // count how often this element exists in global section (==> at this point there are two copies of the element if it was moved)
    QStringList existingNames = getChildrenNames();
    int count = 0;
    for(QString name : existingNames)
      if(name == element->name())
        count++;

    // only remove links if element was not moved
    if(count <= 1)
    {
      // replace selected links to global element with copy of global element
      QList<TreeElement*> list;
      tree->rootElement()->getLinkedList(element->label(), list);
      for(QList<TreeElement*>::iterator iter=list.begin(); iter!=list.end(); iter++)
      {
        XmlNodePtr xmlNode = element->getXML();
        if(xmlNode != nullptr)
        {
          xmlNode->getAttribute("label"); // remove label
          XmlNodePtr xmlNodeOld = (*iter)->getXML();
          xmlNode->addAttribute(xmlNodeOld->getAttribute("loop"));
          xmlNode->addAttribute(xmlNodeOld->getAttribute("condition"));
          (*iter)->parentElement->overwrite(*iter, element->type(), xmlNode);
        }
      }

      // remove unselected links to global element from all elements
      for(int i = children[0].size()-1; i >= 0; i--)
        if(children[0][i] == element)
          tree->rootElement()->removeLink(element);

      // remove loop from all elements in tree
      if(element->type() == "loopType")
      {
        QList<TreeElement*> list;
        tree->rootElement()->getLoopList(element->name(), list);
        for(QList<TreeElement*>::iterator iter=list.begin(); iter!=list.end(); iter++)
          if((*iter)->loop() == element->name())
            (*iter)->setLoop("");
      }

      // remove condition from all elements in tree
      if(element->type() == "conditionType")
      {
        QList<TreeElement*> list;
        tree->rootElement()->getConditionList(element->name(), list);
        for(QList<TreeElement*>::iterator iter=list.begin(); iter!=list.end(); iter++)
          if((*iter)->condition() == element->name())
            (*iter)->setCondition("");
      }
    }

    // remove global element
    if(!TreeElementComplex::removeChild(element))
      throw(Exception("Cannot remove child from global section"));

    tree->undoStack()->endMacro();

    return true;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/
/***********************************************/
