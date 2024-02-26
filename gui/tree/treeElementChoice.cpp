/***********************************************/
/**
* @file treeElementChoice.cpp
*
* @brief Choice element with children.
*
* @author Torsten Mayer-Guerr
* @author Sebastian Strasser
* @date 2006-06-10
*/
/***********************************************/

#include <QtDebug>
#include <QComboBox>
#include "base/importGroops.h"
#include "tree/tree.h"
#include "tree/treeItem.h"
#include "tree/treeElementComplex.h"
#include "tree/treeElementChoice.h"

/***********************************************/

TreeElementChoice::TreeElementChoice(Tree *tree, TreeElementComplex *parentElement, XsdElementPtr xsdElement, const QString &defaultOverride,
                                     XmlNodePtr xmlNode, bool fillWithDefaults, bool recieveAutoComments)
  : TreeElementComplex(tree, parentElement, xsdElement, defaultOverride, xmlNode, recieveAutoComments)
{
  try
  {
    if(optional() && !unbounded())
    {
      addChoice("<none>", XsdElementPtr(nullptr), QJsonObject());
      schemaNameList.push_back("<none>");
      annotationList.push_back("");
    }

    // create children list from XSD schema
    for(auto &xsdElement : xsdElement->complex->elements)
    {
      addChoice(xsdElement->names.front(), xsdElement, defaultObject.value(xsdElement->names.front()).toObject());
      schemaNameList.push_back(xsdElement->names.front());
      annotationList.push_back(xsdElement->annotation);
    }

    int index = 0;
    if(isLinked())
      index = selectedIndex();
    else if(xmlNode) // values from XML
    {
      XmlNodePtr xmlNode2 = xmlNode->getNextChild();
      if(!xmlNode2)
        throw(Exception(xmlNode->getName().toStdString()+": xml node doesn't match with schema"));
      XsdElementPtr xsdChoice = xsdElement->complex->getXsdElement(xmlNode2->getName()); // find also renamed
      if(xsdChoice)
        index = findValueIndex(xsdChoice->names.front());
      else // unkown element
      {
        index = addChoice(xmlNode2->getName(), XsdElementPtr(nullptr), QJsonObject());
        schemaNameList.push_back("[unknown]");
        annotationList.push_back("");
      }
      _valueList[index] = xmlNode2->getName(); // if renamed
      createChildrenElements(index, xmlNode2);
    }
    else if(fillWithDefaults && !defaultObject.isEmpty()) // from default
      index = std::max(findValueIndex(defaultObject.begin().key()), 0);

    TreeElementComplex::setSelectedIndex(index);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

bool TreeElementChoice::isSelectionRenamedInSchema(int index) const
{
  return (index < schemaNameList.size()) && (_valueList[index] != schemaNameList[index]) && (schemaNameList[index] != "[unknown]");
}

/***********************************************/

bool TreeElementChoice::isSelectionUnknown(int index) const
{
  return (index < schemaNameList.size()) && (schemaNameList[index] == "[unknown]");
}

/***********************************************/

XmlNodePtr TreeElementChoice::createXmlTree(bool /*createRootEvenIfEmpty*/) const
{
  try
  {
    if(selectedValue()=="<none>")
      return XmlNodePtr(nullptr);
    XmlNodePtr xmlNode = createXmlBaseNode();
    if(isLinked())
      return xmlNode;
    XmlNodePtr xmlNodeChoice = XmlNode::create(selectedValue());
    xmlNode->addChild(xmlNodeChoice);
    createXmlChildren(xmlNodeChoice);
    return xmlNode;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

void TreeElementChoice::setSelectedIndex(int index)
{
  try
  {
    TreeElementComplex::setSelectedIndex(index);

    if(item())
    {
      if(item() == tree->selectedItem())
        item()->updateAnnotation((index < annotationList.size()) ? annotationList[index] : QString());
      else
        item()->updateAnnotation(annotation());
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

bool TreeElementChoice::overwrite(const QString &type, XmlNodePtr xmlNode, bool contentOnly)
{
  try
  {
    if(!canOverwrite(type))
      return false;

    if(!xmlNode && optional() && !unbounded()) // <none>
    {
      tree->undoStack->beginMacro("overwrite "+name());
      changeSelectedIndex(0);
      tree->undoStack->endMacro();
      return true;
    }

    if(!xmlNode || !xmlNode->hasChildren())
      return false;

    tree->undoStack->beginMacro("overwrite "+name());
    if(baseOverwrite(xmlNode, contentOnly))
    {
      XmlNodePtr xmlNode2 = xmlNode->getNextChild();
      int index = findValueIndex(xmlNode2->getName());
      if(index < 0) // possible renamed
      {
        XsdElementPtr xsdChoice = xsdElement->complex->getXsdElement(xmlNode2->getName()); // find also renamed
        if(xsdChoice)
          for(auto &name : xsdChoice->names)
          {
            index = findValueIndex(name);
            if(index >= 0)
              break;
          }
      }
      if(index < 0) // unknown element
      {
        index = addChoice(xmlNode2->getName(), XsdElementPtr(nullptr), QJsonObject());
        schemaNameList.push_back(xmlNode2->getName());
        annotationList.push_back("");
      }
      changeSelectedIndex(index);
      overwriteChildren(xmlNode2);
    }
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

bool TreeElementChoice::canUpdateName() const
{
  return TreeElementComplex::canUpdateName() || isSelectionRenamedInSchema(selectedIndex());
}

/***********************************************/

void TreeElementChoice::updateName()
{
  if(canUpdateName())
    tree->undoStack->push(new UndoCommandUpdateName(this, isSelectionRenamedInSchema(selectedIndex()) ? schemaNameList[selectedIndex()] : QString()));
}

/***********************************************/
/***********************************************/

QWidget *TreeElementChoice::createEditor()
{
  try
  {
    QComboBox *comboBox = createComboBox(false);
    if(!comboBox)
      return comboBox;
    for(int i=0; i<comboBox->count(); i++)
    {
      if(isSelectionUnknown(i))
        comboBox->setItemIcon(i, QIcon(":/icons/scalable/element-unknown.svg"));
      if(isSelectionRenamedInSchema(i))
        comboBox->setItemIcon(i, QIcon(":/icons/scalable/edit-rename.svg"));
    }
    connect(comboBox, SIGNAL(highlighted(int)), this, SLOT(comboBoxHighlighted(int)));

    return comboBox;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

void TreeElementChoice::comboBoxHighlighted(int index)
{
  if(item())
    item()->updateAnnotation((index < annotationList.size()) ? annotationList[index] : QString());
}

/***********************************************/

void TreeElementChoice::startSelected()
{
  if(item())
    item()->updateAnnotation((selectedIndex() < annotationList.size()) ? annotationList[selectedIndex()] : QString());
}

/***********************************************/


void TreeElementChoice::stopSelected()
{
  if(item())
    item()->updateAnnotation(annotation());
}

/***********************************************/
