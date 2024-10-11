/***********************************************/
/**
* @file treeElementSimple.h
*
* @brief Element without children.
*
* @author Torsten Mayer-Guerr
* @author Sebastian Strasser
* @date 2006-06-10
*/
/***********************************************/

#ifndef __GROOPSGUI__TREEELEMENTSIMPLE__
#define __GROOPSGUI__TREEELEMENTSIMPLE__

#include "base/importGroops.h"
#include "tree/treeElement.h"

/***** CLASS ***********************************/

class TreeElementSimple : public TreeElement
{
  Q_OBJECT

  QString result; // parsed value

public:
  TreeElementSimple(Tree *tree, TreeElementComplex *parentElement, XsdElementPtr xsdElement,
                    const QString &defaultOverride, XmlNodePtr xmlNode, bool fillWithDefaults);

  /** @brief Generate XML-tree. */
  XmlNodePtr createXmlTree(bool createRootEvenIfEmpty) const override;

  /** @brief Values can be edited. */
  bool isEditable() const override {return true;}

  /** @brief the selected value as result of parser. */
  QString selectedResult() const override {return result;}

  /** @brief changes the current index.
  * calls TreeElement::selectIndex
  * Updates the parsed result. */
  void setSelectedIndex(int index) override;

  /** @brief inform this element about changed variables.
  * recursively called for all children.
  * If this element is a variable and @a addVariableInReturn an updated varList is returned. */
  VariableListPtr updateParserResults(VariableListPtr varList, Bool addVariableInReturn) override;

protected:
  virtual QString parseExpression(const QString &text, const VariableList &varList) const;

public:
  /** @brief Is it possible to overweite the element? */
  bool canOverwrite(const QString &type) override {return (this->type() == type);}

  /** @brief Copy the content of @a xmlNode into this.
  * Is undoable.
  * @return success */
  bool overwrite(const QString &type, XmlNodePtr xmlNode, bool contentOnly=false) override;

  /** @brief creates an editable combo box. */
  QWidget *createEditor() override {return createComboBox(true);}
};

/***********************************************/

#endif
