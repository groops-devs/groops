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
#include "addVariableDialog/addVariableDialog.h"
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

    tree->elementGlobal = this;

    addChoice("", nullptr/*xsdElement*/, QJsonObject());
    if(xmlNode)
      createChildrenElements(0, xmlNode);     // add children
    // the add element
    children_[0].push_back(new TreeElementAdd(tree, this, xsdElement, "", true/*visible*/));
    setSelectedIndex(0);
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
/***********************************************/

void TreeElementGlobal::updateParserResultsInScope()
{
  try
  {
    if(parentElement)
      parentElement->updateParserResultsInScope(); // global variables are valid globally
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

VariableListPtr TreeElementGlobal::updateParserResults(VariableListPtr varListOld, Bool /*addVariableInReturn*/)
{
  try
  {
    varList = std::make_shared<const VariableList>(*varListOld); // make copy
    // as the order is irrelevant -> add all variables before
    for(auto child : children())
      if(!child->label().isEmpty() && !child->disabled())
        std::const_pointer_cast<VariableList>(varList)->setVariable(child->label().toStdString(), (child->isLinked() ? "{"+child->selectedValue()+"}" : child->selectedValue()).toStdString());
    for(auto child : children_[selectedIndex()])
      child->updateParserResults(varList, false);
    return varList;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

void TreeElementGlobal::updateLinksInScope()
{
  try
  {
    if(parentElement)
      parentElement->updateLinksInScope(); // global variables are valid globally
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

void TreeElementGlobal::updateLinks(QMap<QString, QString> &labelTypes)
{
  try
  {
    // as the order is irrelevant -> add all links before
    for(auto child : children())
      if(!child->label().isEmpty() && !child->disabled())
        labelTypes[child->label()] = getLinkType(child->type());
    this->labelTypes = labelTypes;
    for(auto child : children_[selectedIndex()])
      child->updateLinks(labelTypes);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/
/***********************************************/

bool TreeElementGlobal::canAddChild(TreeElement */*targetElement*/, const QString &type, const QString &/*label*/) const
{
  return type.isEmpty() || (type == "COMMENT") || tree->xsdElement(type);
}

/***********************************************/

TreeElement *TreeElementGlobal::addChild(TreeElement *targetElement, const QString &type, const QString &label, XmlNodePtr xmlNode)
{
  try
  {
    if(!canAddChild(targetElement, type, label))
      return nullptr;

    if(type == "COMMENT")
      return TreeElementComplex::addChild(targetElement, type, "", xmlNode);

    QString newLabel = label;
    QString newType  = type;
    if(newType.isEmpty())
    {
      AddVariableDialog dialog(tree, "", "", true/*disablePlace*/, dynamic_cast<TreeElementAdd*>(targetElement)/*disableCreateLink*/, tree);
      if(!dialog.exec())
        return nullptr;
      newLabel = dialog.label();
      newType  = dialog.type();
    }

    // check label
    if(newLabel.isEmpty())
    {
      if(xmlNode)
        newLabel = xmlNode->getName(); // default name
      bool ok = true;
      newLabel = QInputDialog::getText(tree, tr("Add global element - GROOPS"), tr("Name of global element:"), QLineEdit::Normal, newLabel, &ok);
      if(!ok)
        return nullptr;
    }

    // check label
    QStringList existingNames;
    for(auto child : children())
      if(!child->label().isEmpty())
        existingNames.push_back(child->label());
    QRegularExpression regex("^[a-zA-Z]([a-zA-Z0-9])*$");
    bool ok = true;
    while(ok && (newLabel.isEmpty() || existingNames.contains(newLabel) || !regex.match(newLabel).hasMatch()))
      newLabel = QInputDialog::getText(tree, tr("Add global element - GROOPS"), tr("Name already exists or is invalid (only letters and digits allowed)!\nChoose another name:"), QLineEdit::Normal, newLabel, &ok);
    if(!ok)
      return nullptr;

    if(xmlNode)
    {
      xmlNode->setName(newType);
      // remove loop, condition, ...
      QString loopLabel, conditionLabel, tmp;
      readAttribute(xmlNode, "label",     tmp);
      readAttribute(xmlNode, "disabled",  tmp);
      readAttribute(xmlNode, "loop",      loopLabel);
      readAttribute(xmlNode, "condition", conditionLabel);
      if(loopLabel      == "_localLoop_")      xmlNode->getChild("loopType");
      if(conditionLabel == "_localCondition_") xmlNode->getChild("conditionType");
      writeAttribute(xmlNode, "label", newLabel);
    }

    return TreeElementComplex::addChild(targetElement, newType, newLabel, xmlNode);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/
