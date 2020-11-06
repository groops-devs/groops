/***********************************************/
/**
* @file treeElementProgram.h
*
* @brief Program element with children.
*
* @author Torsten Mayer-Guerr
* @author Sebastian Strasser
* @date 2006-06-10
*/
/***********************************************/

#ifndef __GROOPSGUI__TREEELEMENTPROGRAM__
#define __GROOPSGUI__TREEELEMENTPROGRAM__

#include "base/importGroops.h"
#include "tree/treeElementChoice.h"
#include <QPushButton>

/***** CLASS ***********************************/

class TreeElementProgram : public TreeElementChoice
{
  Q_OBJECT

  QPushButton *openButton;

public:
  TreeElementProgram(Tree *tree, TreeElementComplex *parentElement, XsdElementPtr xsdElement,
                       const QString &defaultOverride, XmlNodePtr xmlNode, Bool fromFile)
    : TreeElementChoice(tree, parentElement, xsdElement, defaultOverride, xmlNode, fromFile, true/*recieveAutoComments*/) {}
  virtual ~TreeElementProgram() override {}

  /** @brief creates an uneditable combo box with addtional selector. */
  virtual QWidget *createEditor() override;

  /** @brief Opens program selector dialog. */
  virtual void interact() override;

  /** @brief collects all program in a list. */
  virtual void getProgramList(QList<TreeElement*> &list) override {list.push_back(this);}

protected:
  /** @brief changes the current index.
  * if selected index is changed:
  * - update button text
  * - call TreeElementChoice::setSelectedIndex() */
  void setSelectedIndex(int index) override;

private slots:
  void openClicked();
};

/***********************************************/

#endif
