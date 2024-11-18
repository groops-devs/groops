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
                                         const QString &defaultOverride, XmlNodePtr xmlNode, bool fillWithDefaults)
  : TreeElementSimple(tree, parentElement, xsdElement, defaultOverride, xmlNode, fillWithDefaults)
{
}

/***********************************************/

QWidget *TreeElementFileName::createEditor()
{
  try
  {
    // create layout
    QWidget *layoutWidget = new QWidget(tree);
    QHBoxLayout *layout   = new QHBoxLayout(layoutWidget);
    layout->setContentsMargins(0, 0, 0, 0);

    // create ComboBox
    QComboBox *comboBox = createComboBox(true);
    layout->addWidget(comboBox);
    layoutWidget->setFocusProxy(comboBox);

    // create FileSelector Button
    openFileButton = new QPushButton(layoutWidget);
    openFileButton->setIcon(QIcon(":/icons/scalable/document-open.svg"));
    layout->addWidget(openFileButton);
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
    // if directory doesn't exist, repeatedly go up one directory until it exists
    QFileInfo startFile(tree->addWorkingDirectory(selectedResult()));
    while(!startFile.isFile() && !startFile.absoluteDir().exists() && !startFile.absoluteDir().isRoot())
    {
      QString dir = startFile.absolutePath();
      dir.truncate(dir.lastIndexOf("/")+1);
      startFile = QFileInfo(dir + startFile.fileName());
    }

    // lambda to replace beginning part of path with global variable
    auto replaceByVariable = [&](QString path, QString &variable, QString &parsedVariable)
    {
      for(int i=_valueCount; i<_valueList.size(); i++) // all possible links
      {
        QString parsed = parseExpression("{"+_valueList.at(i)+"}", *varList);
        if(path.startsWith(parsed) && (parsed.size() > parsedVariable.size()))
        {
          variable = "{"+_valueList.at(i)+"}";
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

      QString lastFile = tree->stripWorkingDirectory(files.last());
      QString variable;
      QString parsedVariable;
      replaceByVariable(lastFile, variable, parsedVariable);

      for(int i=0; i<files.size()-1; i++)
      {
        XmlNodePtr xmlNode = createXmlTree(true);
        if(!(xmlNode && parentElement))
          continue;
        xmlNode->setText(variable+tree->stripWorkingDirectory(files[i]).mid(parsedVariable.size()));
        parentElement->addChild(this, type(), "", xmlNode);
      }
      changeSelectedValue(variable+lastFile.mid(parsedVariable.size()));
    }
    else
    {
      // File Selector Dialog
      QString selectedPath = QFileDialog::getSaveFileName(openFileButton, name(), startFile.absoluteFilePath(), "", nullptr, QFileDialog::DontConfirmOverwrite);
      if(!selectedPath.isEmpty())
      {
        selectedPath = tree->stripWorkingDirectory(selectedPath);
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
