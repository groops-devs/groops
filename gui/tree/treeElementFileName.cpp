/***********************************************/
/**
* @file treeElementFileName.cpp
*
* @brief Element with file selector.
*
* @author Torsten Mayer-Guerr
* @author Sebastian Strasser
* @date 2006-06-10
*/
/***********************************************/

#include <QtDebug>
#include <QHBoxLayout>
#include <QPushButton>
#include <QFileDialog>
#include "base/importGroops.h"
#include "tree/tree.h"
#include "tree/treeElement.h"
#include "tree/treeElementComplex.h"
#include "tree/treeElementGlobal.h"
#include "tree/treeElementFileName.h"

/***********************************************/

TreeElementFileName::TreeElementFileName(Tree *tree, TreeElementComplex *parentElement, XsdElementPtr xsdElement,
                                         const QString &defaultOverride, XmlNodePtr xmlNode, Bool fromFile)
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
      else if(fromFile || defaultValue().isEmpty())
        insertNewValue("", false);
    }
    if(!defaultValue().isEmpty())
      insertNewValue(defaultValue(), false);

    setSelectedIndex((isLinked() && selectedIndex() > 0) ? selectedIndex() : 0);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

XmlNodePtr TreeElementFileName::getXML(Bool withEmptyNodes) const
{
  try
  {
    if(selectedValue().isEmpty() && !withEmptyNodes)
       return XmlNodePtr(nullptr);
    XmlNodePtr xmlNode = TreeElement::getBaseXML();
    if((xmlNode==nullptr) || isLinked())
      return xmlNode;
    xmlNode->setText(selectedValue());
    return xmlNode;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

QWidget *TreeElementFileName::createEditor()
{
  try
  {
    // create layout
    QWidget *layoutWidget = new QWidget(tree);
    QHBoxLayout *layout   = new QHBoxLayout(layoutWidget);
    layout->setMargin(0);

    // create ComboBox
    QComboBox *comboBox = createComboBox(true);
    layout->addWidget(comboBox);
    layoutWidget->setFocusProxy(comboBox);

    // create FileSelector Button
    openFileButton = new QPushButton(layoutWidget);
    openFileButton->setIcon(QIcon(":/icons/scalable/document-open.svg"));
    layout->addWidget( openFileButton );
    // signals and slots connections
    connect(openFileButton, SIGNAL(clicked()), this, SLOT(openFileClicked()));

    return layoutWidget;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

void TreeElementFileName::interact()
{
    openFileClicked();
}

/***********************************************/

void TreeElementFileName::openFileClicked()
{
  try
  {
    QFileInfo startFile(tree->addXmlDirectory(selectedResult()));

    // if directory doesn't exist, repeatedly go up one directory until it exists
    while(!startFile.isFile() && !startFile.absoluteDir().exists() && !startFile.absoluteDir().isRoot())
    {
      QString dir = startFile.absolutePath();
      dir.truncate(dir.lastIndexOf("/")+1);
      startFile = QFileInfo(dir + startFile.fileName());
    }

    // lambda to replace beginning part of path with global variable
    auto replaceByVariable = [&](QString path, QString &variable, QString &parsedVariable)
    {
      for(int i = valueCount(); i < valueList().size(); i++)
      {
        QString parsed = parseExpression(valueList().at(i));
        if(path.startsWith(parsed) && parsed.size() > parsedVariable.size())
        {
          variable = "{"+valueList().at(i)+"}";
          parsedVariable = parsed;
        }
      }
    };

    // possible to add several elements?
    if(elementAdd())
    {
      QStringList files = QFileDialog::getOpenFileNames(openFileButton, name(), startFile.absoluteFilePath());
      if(files.size()==0)
        return;

      QString lastFile = tree->stripXmlDirectory(files.last());
      QString variable;
      QString parsedVariable;
      replaceByVariable(lastFile, variable, parsedVariable);

      for(int i=0; i<files.size()-1; i++)
      {
        XmlNodePtr xmlNode = getXML(true);
        if(!(xmlNode && parentElement))
          continue;
        xmlNode->setText(variable+tree->stripXmlDirectory(files[i]).mid(parsedVariable.size()));
        parentElement->addChild(this, type(), xmlNode);
      }
      changeSelectedValue(variable+lastFile.mid(parsedVariable.size()));
    }
    else
    {
      // File Selector Dialog
      QString selectedPath = QFileDialog::getSaveFileName(openFileButton, name(), startFile.absoluteFilePath(), "", nullptr, QFileDialog::DontConfirmOverwrite);
      if(!selectedPath.isEmpty())
      {
        selectedPath = tree->stripXmlDirectory(selectedPath);
        QString variable;
        QString parsedVariable;
        replaceByVariable(selectedPath, variable, parsedVariable);
        changeSelectedValue(variable+selectedPath.mid(parsedVariable.size()));
      }
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/
