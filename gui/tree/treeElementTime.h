/***********************************************/
/**
* @file treeElementSimple.h
*
* @brief Element with date/time editor.
*
* @author Torsten Mayer-Guerr
* @author Sebastian Strasser
* @date 2006-06-10
*/
/***********************************************/

#ifndef __GROOPSGUI__TREEELEMENTTIME__
#define __GROOPSGUI__TREEELEMENTTIME__

#include "base/importGroops.h"
#include "tree/treeElement.h"
#include "tree/treeElementSimple.h"

/***** TYPES ***********************************/

class QComboBox;
class QDateTimeEdit;
class QDateTime;

/***** CLASS ***********************************/

class TreeElementTime : public TreeElementSimple
{
  Q_OBJECT

  VariableList varList;

public:
  TreeElementTime(Tree *tree, TreeElementComplex *parentElement, XsdElementPtr xsdElement,
                  const QString &defaultOverride, XmlNodePtr xmlNode, bool fillWithDefaults)
    : TreeElementSimple(tree, parentElement, xsdElement, defaultOverride, xmlNode, fillWithDefaults), mjd(0.) {}

  /** @brief changes the current index.
   * calls TreeElementSimple::selectIndex
  * Updates the timeEditor. */
  void setSelectedIndex(int index) override;

  /** @brief inform this element about changed variables.
  * recursively called for all children. */
  void updateParserResults(VariableList &varList) override {this->varList = varList; TreeElementSimple::updateParserResults(varList);}

  /** @brief creates an editable combo box + dateTimeEdit. */
  QWidget *createEditor() override;

  /** @brief Switches focus between comboBox and dateTimeEdit. */
  void interact() override;

private:
  mutable Double          mjd; // updated by parseExpression
  QPointer<QComboBox>     comboBox;
  QPointer<QDateTimeEdit> dateTimeEdit;
  bool changeNotComboBox, changeNotDateTime;

  QString   date2mjd(const QDateTime &dateTime) const;
  QDateTime mjd2date(Double mjd) const;

  QString parseExpression(const QString &text, const VariableList &varList) const override;

private slots:
  void comboBoxEditTextChanged(const QString &text);
  void dateTimeChanged(const QDateTime &dateTime);
};

/***********************************************/

#endif
