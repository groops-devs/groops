/***********************************************/
/**
* @file treeElementChoice.h
*
* @brief Choice element with children.
*
* @author Torsten Mayer-Guerr
* @author Sebastian Strasser
* @date 2006-06-10
*/
/***********************************************/

#ifndef __GROOPSGUI__TREEELEMENTCHOICE__
#define __GROOPSGUI__TREEELEMENTCHOICE__

#include "base/importGroops.h"
#include "tree/treeElementComplex.h"

/***** CLASS ***********************************/

class TreeElementChoice : public TreeElementComplex
{
  Q_OBJECT

  QStringList annotationList;
  QStringList renamedList;

public:
  TreeElementChoice(Tree *tree, TreeElementComplex *parentElement, XsdElementPtr xsdElement, const QString &defaultOverride,
                    XmlNodePtr xmlNode, Bool fromFile, Bool recieveAutoComments=false);
  virtual ~TreeElementChoice() override {}

/** @brief Is the selection renamed, i.e. has the selection been renamed in the schema? */
virtual Bool isSelectionRenamed(int index) const override;

/** @brief Returns the new name of the current selection in case the selection has been renamed in the schema. */
virtual QString renamedSelection() const;

/** @brief Generate XML-tree.
* recursively called for all children. */
virtual XmlNodePtr getXML(Bool withEmptyNodes=false) const override;

/** @brief event handler of current index change.
* This event handler is called by setSelectedValue whenever the value is changed. */
virtual void newSelectedIndex(int index) override;

/** @brief creates an uneditable combo box. */
virtual QWidget *createEditor() override;

/** @brief event handler called by item when it gets selected. */
void startSelected() override;

/** @brief event handler called by item when it losts the selection. */
void stopSelected() override;

QString annotationChild(int index) const;

protected slots:
  void comboBoxHighlighted(int index);
};

/***********************************************/

#endif
