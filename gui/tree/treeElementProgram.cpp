/***********************************************/
/**
* @file treeElementProgram.cpp
*
* @brief Program element with children.
*
* @author Torsten Mayer-Guerr
* @author Sebastian Strasser
* @date 2006-06-10
*/
/***********************************************/

#include "base/importGroops.h"
#include "tree/tree.h"
#include "tree/treeElement.h"
#include "tree/treeElementComplex.h"
#include "programDialog/programDialog.h"
#include "tree/treeElementProgram.h"

/***********************************************/

QWidget *TreeElementProgram::createEditor()
{
  try
  {
    // create program selection button
    openButton = new QPushButton(tree);
    openButton->setText(selectedValue());
    openButton->setStyleSheet("Text-align:left; padding: 2px; padding-left: 4px");
    openButton->setDefault(true);
    if(isSelectionRenamedInSchema(selectedIndex()))
      openButton->setIcon(QIcon(":/icons/scalable/edit-rename.svg"));
    else if(isSelectionUnknown(selectedIndex()))
      openButton->setIcon(QIcon(":/icons/scalable/element-unknown.svg"));

    // signals and slots connections
    connect(openButton, SIGNAL(clicked()), this, SLOT(openClicked()));

    return openButton;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

void TreeElementProgram::interact()
{
  openClicked();
}

/***********************************************/

void TreeElementProgram::setSelectedIndex(int index)
{
  TreeElementChoice::setSelectedIndex(index);
  if(openButton)
  {
    if(isSelectionRenamedInSchema(selectedIndex()))
      openButton->setIcon(QIcon(":/icons/scalable/edit-rename.svg"));
    else if(isSelectionUnknown(selectedIndex()))
      openButton->setIcon(QIcon(":/icons/scalable/element-unknown.svg"));
    else
      openButton->setIcon(QIcon());
    openButton->setText(selectedValue());
  }
}

/***********************************************/

void TreeElementProgram::openClicked()
{
  try
  {
    ProgramDialog dialog(this, tree);
    dialog.exec();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

