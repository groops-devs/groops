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
                                     XmlNodePtr xmlNode, Bool fromFile, Bool recieveAutoComments)
  : TreeElementComplex(tree, parentElement, xsdElement, defaultOverride, xmlNode, recieveAutoComments)
{
  try
  {
    if((optional()) && (!unbounded()))
    {
      addChoice("<none>", XsdElementPtr(nullptr), QJsonObject());
      annotationList.push_back("");
    }

    // create children list from XSD schema
    XsdComplexPtr complex = xsdElement->complex;
    for(UInt i=0; i<complex->element.size(); i++)
    {
      addChoice(complex->element.at(i)->name, complex->element.at(i), defaultObject.value(complex->element.at(i)->name).toObject());
      annotationList.push_back(complex->element.at(i)->annotation);
      renamedList.push_back(QString());
    }

    // values from XML
    int index = 0;
    if(xmlNode && !isLinked())
    {
      XmlNodePtr xmlNode2 = xmlNode->getNextChild();
      if(!xmlNode2)
        throw(Exception("xml node doesn't match with schema"));
      else
      {
        index = findValueIndex(xmlNode2->getName());
        if(index < 0) // unknown element
        {
          XsdElementPtr xsdElement = complex->getXsdElement(xmlNode2->getName());
          index = addChoice(xmlNode2->getName(), xsdElement, QJsonObject());
          annotationList.push_back("");
          renamedList.push_back(xsdElement ? xsdElement->name : QString());
        }
      }
      createChildrenElements(index, xmlNode2);
    }

    if(isLinked())
      setSelectedIndex(selectedIndex());
    else if(!fromFile && !defaultObject.isEmpty())
    {
      int index = findValueIndex(defaultObject.begin().key());
      setSelectedIndex(index > 0 ? index : 0);
    }
    else
      setSelectedIndex(index);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

Bool TreeElementChoice::isSelectionRenamed(int index) const
{
  try
  {
    if(index >= 0 && index < renamedList.size())
      return !renamedList.at(index).isEmpty() && renamedList.at(index) != valueList().at(index);

    return false;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

QString TreeElementChoice::renamedSelection() const
{
  try
  {
    if(!isSelectionRenamed(selectedIndex()))
      return "";

    return renamedList.at(selectedIndex());
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

XmlNodePtr TreeElementChoice::getXML(Bool withEmptyNodes) const
{
  try
  {
    if(selectedValue()=="<none>")
      return XmlNodePtr(nullptr);
    XmlNodePtr xmlNode = getBaseXML();
    if((xmlNode==nullptr) || isLinked())
      return xmlNode;
    XmlNodePtr xmlNode2 = createXmlNode(xmlNode, selectedValue());
    getChildrenXML(xmlNode2, withEmptyNodes);
    return xmlNode;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

QString TreeElementChoice::annotationChild(int index) const
{
  if(index>=annotationList.size())
    return QString();
  return annotationList[index];
}

/***********************************************/

void TreeElementChoice::newSelectedIndex(int index)
{
  try
  {
    TreeElementComplex::newSelectedIndex(index);

    if(item())
    {
      if(item() == tree->currentItem())
        item()->updateAnnotation(annotationChild(index));
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
/***********************************************/

QWidget *TreeElementChoice::createEditor()
{
  try
  {
    QComboBox *comboBox = createComboBox(false);
    if(comboBox == nullptr)
      return comboBox;
    for(int i = 0; i < comboBox->count(); i++)
    {
      if(isSelectionRenamed(i))
        comboBox->setItemIcon(i, QIcon(":/icons/scalable/edit-rename.svg"));
      if(isSelectionUnknown(i))
        comboBox->setItemIcon(i, QIcon(":/icons/scalable/element-unknown.svg"));
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
    item()->updateAnnotation(annotationChild(index));
}

/***********************************************/

void TreeElementChoice::startSelected()
{
  if(item())
    item()->updateAnnotation(annotationChild(selectedIndex()));
}

/***********************************************/


void TreeElementChoice::stopSelected()
{
  if(item())
    item()->updateAnnotation(annotation());
}

/***********************************************/
