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
#include "tree/treeElement.h"

/***** TYPES ***********************************/

class QPushButton;

/***** CLASS ***********************************/

class TreeElementFileName : public TreeElement
{
  Q_OBJECT

public:
  TreeElementFileName(Tree *tree, TreeElementComplex *parentElement, XsdElementPtr xsdElement,
                      const QString &defaultOverride, XmlNodePtr xmlNode, Bool fromFile);
  virtual ~TreeElementFileName() override {}

/** @brief Generate XML-tree. */
virtual XmlNodePtr getXML(Bool withEmptyNodes=false) const override;

/** @brief Values can be edited. */
virtual Bool isEditable() const override {return true;}

/** @brief creates an editable combo box with addtional file selector. */
virtual QWidget *createEditor() override;

/** @brief Opens file selector dialog. */
virtual void interact() override;

private:
  QPointer<QPushButton> openFileButton;

private slots:
  void openFileClicked();
};

/***********************************************/

#endif
