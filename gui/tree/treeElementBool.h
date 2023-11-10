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
#include "tree/treeElementSimple.h"

/***** CLASS ***********************************/

class TreeElementBool : public TreeElement
{
public:
  TreeElementBool(Tree *tree, TreeElementComplex *parentElement, XsdElementPtr xsdElement,
                  const QString &defaultOverride, XmlNodePtr xmlNode, bool fillWithDefaults);

  /** @brief Values cannot be edited. */
  bool isEditable() const override {return false;}

  /** @brief Generate XML-tree. */
  XmlNodePtr createXmlTree(bool /*createRootEvenIfEmpty*/) const override;

  /** @brief Is it possible to overweite the element? */
  bool canOverwrite(const QString &type) override {return (this->type() == type);}

  /** @brief Copy the content of @a xmlNode into this.
  * Is undoable.
  * @return success */
  bool overwrite(const QString &type, XmlNodePtr xmlNode, bool contentOnly=false) override;

  /** @brief creates an uneditable combo box. */
  QWidget *createEditor() override {return createComboBox(false);}
};

/***********************************************/

#endif
