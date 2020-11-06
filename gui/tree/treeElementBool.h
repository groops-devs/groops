/***********************************************/
/**
* @file treeElementBool.h
*
* @brief Element with true or false.
*
* @author Torsten Mayer-Guerr
* @author Sebastian Strasser
* @date 2006-06-10
*/
/***********************************************/

#ifndef __GROOPSGUI__TREEELEMENTBOOL__
#define __GROOPSGUI__TREEELEMENTBOOL__

#include "base/importGroops.h"
#include "tree/treeElement.h"

/***** CLASS ***********************************/

class TreeElementBool : public TreeElement
{
public:
  TreeElementBool(Tree *tree, TreeElementComplex *parentElement, XsdElementPtr xsdElement,
                  const QString &defaultOverride, XmlNodePtr xmlNode);
  virtual ~TreeElementBool() override {}

/** @brief Generate XML-tree. */
virtual XmlNodePtr getXML(Bool withEmptyNodes=false) const override;

/** @brief creates an uneditable combo box. */
virtual QWidget *createEditor() override {return createComboBox(false);}
};

/***********************************************/

#endif
