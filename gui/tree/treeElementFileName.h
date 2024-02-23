/***********************************************/
/**
* @file treeElementFileName.h
*
* @brief Element with file selector.
*
* @author Torsten Mayer-Guerr
* @author Sebastian Strasser
* @date 2006-06-10
*/
/***********************************************/

#ifndef __GROOPSGUI__TREEELEMENTFILENAME__
#define __GROOPSGUI__TREEELEMENTFILENAME__

#include "base/importGroops.h"
#include "tree/treeElementSimple.h"

/***** TYPES ***********************************/

class QPushButton;

/***** CLASS ***********************************/

class TreeElementFileName : public TreeElementSimple
{
  Q_OBJECT

  VariableList varList;

public:
  TreeElementFileName(Tree *tree, TreeElementComplex *parentElement, XsdElementPtr xsdElement,
                      const QString &defaultOverride, XmlNodePtr xmlNode, bool fillWithDefaults);

  /** @brief inform this element about changed variables.
  * recursively called for all children. */
  void updateParserResults(VariableList &varList) override {this->varList = varList; TreeElementSimple::updateParserResults(varList);}

  /** @brief creates an editable combo box with additional file selector. */
  QWidget *createEditor() override;

  /** @brief Opens file selector dialog. */
  void interact() override;

private:
  QPointer<QPushButton> openFileButton;

private slots:
  void openFileClicked();
};

/***********************************************/

#endif
