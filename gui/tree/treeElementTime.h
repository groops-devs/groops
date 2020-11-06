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

public:
  TreeElementTime(Tree *tree, TreeElementComplex *parentElement, XsdElementPtr xsdElement,
                  const QString &defaultOverride, XmlNodePtr xmlNode, Bool fromFile)
    : TreeElementSimple(tree, parentElement, xsdElement, defaultOverride, xmlNode, fromFile) {}
  virtual ~TreeElementTime() override {}

/** @brief Values can be edited. */
virtual Bool isEditable() const override {return true;}

/** @brief event handler of current index change.
* This event handler is called by setSelectedValue whenever the value is changed.
* Updates the timeEditor. */
virtual void newSelectedIndex(int index) override;

/** @brief creates an editable combo box + dateTimeEdit. */
virtual QWidget *createEditor() override;

/** @brief Switches focus between comboBox and dateTimeEdit. */
virtual void interact() override;

private:
  QPointer<QComboBox>     comboBox;
  QPointer<QDateTimeEdit> dateTimeEdit;
  Bool changeNotComboBox, changeNotDateTime;

  QString   date2mjd(const QDateTime &dateTime) const;
  QDateTime mjd2date(const QString &text) const;

  virtual QString parseExpression(const QString &value) const override;

private slots:
  void comboBoxEditTextChanged(const QString &text);
  void dateTimeChanged(const QDateTime &dateTime);
};

/***********************************************/

#endif
