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
                                         const QString &defaultOverride, XmlNodePtr xmlNode, bool fillWithDefaults)
  : TreeElementComplex(tree, parentElement, xsdElement, defaultOverride, xmlNode)
{
  try
  {
    int index = 0;
    if(optional() && !unbounded())
    {
      index = addChoice("<none>",    XsdElementPtr(nullptr), QJsonObject());
      index = addChoice("<enabled>", xsdElement,             defaultObject);
      if(!xmlNode || isLinked())
        index = (fillWithDefaults && !defaultObject.isEmpty()) ? 1 : 0;
    }
    else
      index = addChoice("", xsdElement, defaultObject);

    if(xmlNode && !isLinked())
      createChildrenElements(index, xmlNode);
    setSelectedIndex(isLinked() ? selectedIndex() : index);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

XmlNodePtr TreeElementSequence::createXmlTree(bool /*createRootEvenIfEmpty*/) const
{
  try
  {
    if(selectedValue()=="<none>")
      return XmlNodePtr(nullptr);
    XmlNodePtr xmlNode = createXmlBaseNode();
    if(!isLinked())
      createXmlChildren(xmlNode);
    return xmlNode;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

bool TreeElementSequence::overwrite(const QString &type, XmlNodePtr xmlNode, bool contentOnly)
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

    if(!xmlNode)
      return false;

    tree->undoStack->beginMacro("overwrite "+name());
    if(baseOverwrite(xmlNode, contentOnly))
    {
      changeSelectedIndex((optional() && !unbounded()) ? 1 : 0); // change to <enabled>
      overwriteChildren(xmlNode);
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
