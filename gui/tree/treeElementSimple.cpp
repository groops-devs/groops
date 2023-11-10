/***********************************************/
/**
* @file treeElementSimple.cpp
*
* @brief Element without children.
*
* @author Torsten Mayer-Guerr
* @author Sebastian Strasser
* @date 2006-06-10
*/
/***********************************************/

#include <QtDebug>
#include "base/importGroops.h"
#include "tree/tree.h"
#include "tree/treeElement.h"
#include "tree/treeElementComplex.h"
#include "tree/treeElementGlobal.h"
#include "tree/treeItem.h"
#include "tree/treeElementSimple.h"

/***********************************************/

TreeElementSimple::TreeElementSimple(Tree *tree, TreeElementComplex *parentElement, XsdElementPtr xsdElement,
                                     const QString &defaultOverride, XmlNodePtr xmlNode, bool fillWithDefaults)
  : TreeElement(tree, parentElement, xsdElement, defaultOverride, xmlNode)
{
  try
  {
    if(!isLinked())
    {
      if(xmlNode && xmlNode->hasChildren())
        throw(Exception("xml node doesn't match with schema"));
      else if(xmlNode && !xmlNode->getText().isEmpty())
        insertNewValue(xmlNode->getText(), false);
      else if(!fillWithDefaults || defaultValue.isEmpty())
        insertNewValue("", false);
    }
    if(!defaultValue.isEmpty())
      insertNewValue(defaultValue, false);

    setSelectedIndex(isLinked() ? selectedIndex() : 0);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

XmlNodePtr TreeElementSimple::createXmlTree(bool createRootEvenIfEmpty) const
{
  try
  {
    if(selectedValue().isEmpty() && !createRootEvenIfEmpty)
       return XmlNodePtr(nullptr);
    XmlNodePtr xmlNode = TreeElement::createXmlBaseNode();
    if(!isLinked())
      xmlNode->setText(selectedValue());
    return xmlNode;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

void TreeElementSimple::setSelectedIndex(int index)
{
  try
  {
    TreeElement::setSelectedIndex(index);

    if(tree->elementGlobal && tree->rootElement)
      updateParserResults(tree->elementGlobal->variableList(), false/*recursively*/);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

void TreeElementSimple::updateParserResults(const VariableList &varList, bool /*recursively*/)
{
  QString resultNew = parseExpression((isLinked()) ? "{"+selectedValue()+"}" : selectedValue(), varList);
  if(resultNew == result)
    return;
  result = resultNew;
  if(item())
    item()->updateValue();
}

/***********************************************/

QString TreeElementSimple::parseExpression(const QString &text, const VariableList &varList) const
{
  QString result = text;
  try
  {
    bool resolved = true;
    result = QString::fromStdString(StringParser::parse(name().toStdString(), text.toStdString(), varList, resolved));
    if(resolved && QStringList({"bool", "int", "uint", "double", "angle", "time", "expression"}).contains(type())) // only numerical values
    {
      Double d = ExpressionVariable::parse(result.toStdString(), varList);
      result.setNum(d, 'f', 7).remove(QRegularExpression("0+$")).remove(QRegularExpression("\\.$")); // %.7f with trailing zeros removed
    }
  }
  catch(std::exception &/*e*/)
  {
  }
  return result;
}


/***********************************************/

bool TreeElementSimple::overwrite(const QString &type, XmlNodePtr xmlNode, bool contentOnly)
{
  try
  {
    if(!canOverwrite(type) || !xmlNode)
      return false;

    tree->undoStack->beginMacro("overwrite "+name());
    if(baseOverwrite(xmlNode), contentOnly)
      changeSelectedValue(xmlNode->getText());
    tree->undoStack->endMacro();
    return true;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/
