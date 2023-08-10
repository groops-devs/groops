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
#include "addGlobalDialog/addGlobalDialog.h"
#include "base/importGroops.h"
#include "tree/tree.h"
#include "tree/treeElement.h"
#include "tree/treeElementSimple.h"
#include "tree/treeElementAdd.h"
#include "tree/treeElementComplex.h"
#include "tree/treeElementComment.h"
#include "tree/treeElementUnknown.h"
#include "tree/treeItem.h"
#include "tree/treeElementGlobal.h"

/***********************************************/

TreeElementGlobal::TreeElementGlobal(Tree *tree, TreeElementComplex *parentElement,
                                     XsdElementPtr xsdElement, XmlNodePtr xmlNode)
  : TreeElementComplex(tree, parentElement, xsdElement, "", xmlNode, false/*recieveAutoComments*/)
{
  try
  {
    if(isLinked())
      throw(Exception("global cannot be linked"));

    addChoice("", xsdElement, QJsonObject());

    // add children
    if(xmlNode)
      for(XmlNodePtr xmlChild = xmlNode->getNextChild(); xmlChild; xmlChild = xmlNode->getNextChild())
      {
        XsdElementPtr xsdChild = xsdElement->complex->getXsdElement(xmlChild->getName());
        if(xsdChild)
          children_[0].push_back(TreeElement::newTreeElement(tree, this, xsdChild, "", xmlChild, false/*fillWithDefaults*/));
        else
          children_[0].push_back(new TreeElementUnknown(tree, this, xmlChild));
      }

    // sort unknown at end
    std::stable_sort(children_[0].begin(), children_[0].end(), [](auto e1, auto e2) {return (!dynamic_cast<TreeElementUnknown*>(e1) && dynamic_cast<TreeElementUnknown*>(e2));});

    // the add element
    TreeElementAdd *elementAdd = new TreeElementAdd(tree, this, xsdElement, true/*visible*/);
    for(auto &child : children_[0])
      child->setElementAdd(elementAdd);
    children_[0].push_back(elementAdd);

    setSelectedIndex(0);
    tree->elementGlobal = this;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

TreeElementGlobal::~TreeElementGlobal()
{
  tree->elementGlobal = nullptr;
}

/***********************************************/

XmlNodePtr TreeElementGlobal::createXmlTree(bool /*createRootEvenIfEmpty*/) const
{
  try
  {
    XmlNodePtr xmlNode = createXmlBaseNode();
    if(!isLinked())
      createXmlChildren(xmlNode, true/*createRootEvenIfEmpty*/);
    return xmlNode;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

XsdElementPtr TreeElementGlobal::findXsdElement(const QString &type) const
{
  auto iter = std::find_if(xsdElement->complex->elements.begin(), xsdElement->complex->elements.end(),
                           [&type](auto &x) {return x->type == type;});
  return (iter == xsdElement->complex->elements.end()) ? XsdElementPtr(nullptr) : *iter;
}

/***********************************************/
/***********************************************/

void TreeElementGlobal::informAboutGlobalElements(TreeElement *element, bool recursively) const
{
  try
  {
    for(auto elementInGlobal : children_[selectedIndex()])
      if(!dynamic_cast<TreeElementAdd*>(elementInGlobal) && !dynamic_cast<TreeElementComment*>(elementInGlobal))
        element->informAboutLink(elementInGlobal, recursively);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

QStringList TreeElementGlobal::names() const
{
  try
  {
    QStringList childrenNames;
    for(auto elementInGlobal : children_[selectedIndex()])
      if(!dynamic_cast<TreeElementAdd*>(elementInGlobal) && !dynamic_cast<TreeElementComment*>(elementInGlobal))
        childrenNames.push_back(elementInGlobal->name());
    return childrenNames;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

void TreeElementGlobal::updateVariableList()
{
  try
  {
    if(!tree->rootElement)
      return;
    varList = VariableList();
    for(auto elementInGlobal : children_[selectedIndex()])
    {
      TreeElementSimple *element = dynamic_cast<TreeElementSimple*>(elementInGlobal);
      if(element)
        varList.setVariable(element->name().toStdString(),
                            (element->isLinked() ? "{"+element->selectedValue()+"}" : element->selectedValue()).toStdString());
    }

    // inform all elements about changes
    tree->rootElement->updateParserResults(varList, true/*recursively*/);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

bool TreeElementGlobal::checkLabel(const QString &oldLabel, const QString &defaultLabel, QString &label) const
{
  bool ok = true;
  QStringList existingNames = names();
  existingNames.removeAt(existingNames.indexOf(oldLabel));
  if(label.isEmpty() || existingNames.contains(label))
  {
    if(label.isEmpty())
      label = defaultLabel; // default name
    label = QInputDialog::getText(tree, tr("Add global element - GROOPS"), tr("Name of global element:"), QLineEdit::Normal, label, &ok);
    QRegularExpression regex("^[a-zA-Z]([a-zA-Z0-9])*$");
    while(ok && (label.isEmpty() || existingNames.contains(label) || !regex.match(label).hasMatch()))
      label = QInputDialog::getText(tree, tr("Add global element - GROOPS"), tr("Name already exists or is invalid (only letters and digits allowed)!\nChoose another name:"), QLineEdit::Normal, label, &ok);
  }
  return ok;
}

/***********************************************/
/***********************************************/

bool TreeElementGlobal::canSetGlobal(TreeElement *element) const
{
  try
  {
    if(dynamic_cast<TreeElementAdd*>(element))
      return false;
    return (findXsdElement(element->type()) != nullptr);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

bool TreeElementGlobal::setGlobal(TreeElement *element)
{
  try
  {
    if(!canSetGlobal(element))
      return false;

    // new global elements are added to the bottom of the list by default
    TreeElement *targetElement = children_[selectedIndex()].back();
    // if element is already in global, new global element is added directly in front of it
    if(element->parentElement == this)
      targetElement = element;

    tree->undoStack->beginMacro("set global "+element->name());
    // add element to global section
    TreeElement *elementGlobal = addChild(targetElement, element->type(), element->createXmlTree());
    if(!elementGlobal)
    {
      // abort/end macro and overwrite it with empty command
      tree->undoStack->endMacro();
      tree->undoStack->undo();
      tree->undoStack->push(new QUndoCommand);
      tree->undoStack->undo();
      return false;
    }
    // element will be linked
    element->changeSelectedIndex(element->findLinkIndex(elementGlobal->label()));
    tree->undoStack->endMacro();
    return true;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

TreeElement *TreeElementGlobal::addNewChild(QString type, QString label)
{
  try
  {
    if(type.isEmpty() && label.isEmpty())
    {
      AddGlobalDialog dialog(this, tree);
      if(!dialog.exec())
        return nullptr;
      type  = dialog.elementType();
      label = dialog.elementName();
    }

    XsdElementPtr xsdElement = findXsdElement(type);
    if(!xsdElement)
      return nullptr;

    // create xmlNode from new element
    TreeElement *element = TreeElement::newTreeElement(tree, this, xsdElement, "", nullptr, true/*fillWithDefaults*/);
    element->_label = label;
    XmlNodePtr xmlNode = element->createXmlTree(true/*createRootEvenIfEmpty*/);
    delete element;

    return addChild(children_[selectedIndex()].back(), type, xmlNode);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

bool TreeElementGlobal::canAddChild(TreeElement */*targetElement*/, const QString &type) const
{
  return findXsdElement(type) || (type == "COMMENT");
}

/***********************************************/

bool TreeElementGlobal::canRemoveChild(TreeElement *element) const
{
  return element && !dynamic_cast<TreeElementAdd*>(element);
}

/***********************************************/

TreeElement *TreeElementGlobal::addChild(TreeElement *targetElement, const QString &type, XmlNodePtr xmlNode)
{
  try
  {
    if(!canAddChild(targetElement, type) || !xmlNode)
      return nullptr;

    if(!targetElement)
      targetElement = children_[selectedIndex()].back(); // default: add at end

    if(type == "COMMENT")
      return TreeElementComplex::addChild(targetElement, type, xmlNode);

    QString label;
    readAttribute(xmlNode, "label", label);
    if(!checkLabel(QString(), xmlNode->getName(), label))
      return nullptr;
    writeAttribute(xmlNode, "label", label);

     // create & init new element
    TreeElement *newElement = TreeElement::newTreeElement(tree, this, findXsdElement(type), "", xmlNode, false/*fillWithDefaults*/);
    newElement->_name  = newElement->_schemaName;
    newElement->_label = label;
    tree->elementGlobal->informAboutGlobalElements(newElement, false); // inform the new element about all links

    tree->undoStack->beginMacro("add link "+label);
    tree->undoStack->push(new UndoCommandRemoveAddChild(newElement, targetElement, this, true/*isAdd*/));
    tree->rootElement->addedLink(newElement);
    tree->undoStack->endMacro();

    return newElement;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

bool TreeElementGlobal::removeChild(TreeElement *element)
{
  try
  {
    tree->undoStack->beginMacro("remove link "+element->name());
    TreeElementComplex::removeChild(element);
    tree->rootElement->removedLink(element);
    tree->undoStack->endMacro();
    return true;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/
/***********************************************/

class TreeElementGlobal::UndoCommandRename : public TreeElement::UndoCommand
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

void TreeElementGlobal::UndoCommandRename::redo()
{
  try
  {
    std::swap(label, treeElement->_label);
    tree->rootElement->renamedLink(label, treeElement->_label);
    dynamic_cast<TreeElementGlobal*>(treeElement->parentElement)->updateVariableList();

    if(treeElement->item())
    {
      tree->setSelectedItem(treeElement->item());
      treeElement->item()->updateName();
      treeElement->item()->updateIcon();
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

bool TreeElementGlobal::canRenameChild(TreeElement *element) const
{
  if(dynamic_cast<const TreeElementAdd*>(element))
    return false; // the add button cannot be renamed
  return true;    // global elements can always be renamed
}

/***********************************************/

void TreeElementGlobal::renameChild(TreeElement *element, const QString &label)
{
  try
  {
    if(!element || !canRenameChild(element) || (label == element->label()))
      return;

    QString label2 = label;
    if(checkLabel(element->label(), element->label(), label2))
      tree->undoStack->push(new UndoCommandRename(element, label2));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/
/***********************************************/
