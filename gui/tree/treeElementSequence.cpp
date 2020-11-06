/***********************************************/
/**
* @file treeElementSequence.cpp
*
* @brief Element with children.
*
* @author Torsten Mayer-Guerr
* @author Sebastian Strasser
* @date 2006-06-10
*/
/***********************************************/

#include <QtDebug>
#include "base/importGroops.h"
#include "tree/tree.h"
#include "tree/treeItem.h"
#include "tree/treeElementComplex.h"
#include "tree/treeElementSequence.h"

/***********************************************/

TreeElementSequence::TreeElementSequence(Tree *tree, TreeElementComplex *parentElement, XsdElementPtr xsdElement,
                                         const QString &defaultOverride, XmlNodePtr xmlNode, Bool fromFile)
  : TreeElementComplex(tree, parentElement, xsdElement, defaultOverride, xmlNode)
{
  try
  {
    int index = 0;
    if((optional()) && (!unbounded()))
    {
      index = addChoice("<none>",    XsdElementPtr(nullptr), QJsonObject());
      index = addChoice("<enabled>", xsdElement,             defaultObject);
      if(xmlNode && (!isLinked()))
        createChildrenElements(index, xmlNode);
      else
        index = ((!defaultObject.isEmpty()) && !fromFile ? 1 : 0);
    }
    else
    {
      index = addChoice("", xsdElement, defaultObject);
      if(xmlNode && (!isLinked()))
        createChildrenElements(index, xmlNode);
    }

    setSelectedIndex(isLinked() ? selectedIndex() : index);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

XmlNodePtr TreeElementSequence::getXML(Bool withEmptyNodes) const
{
  try
  {
    if(selectedValue()=="<none>")
      return XmlNodePtr(nullptr);
    XmlNodePtr xmlNode = getBaseXML();
    if((xmlNode==nullptr) || isLinked())
      return xmlNode;
    getChildrenXML(xmlNode, withEmptyNodes);
    return xmlNode;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/
