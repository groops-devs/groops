/***********************************************/
/**
* @file addVariableDialog.h
*
* @brief Dialog to add variables.
*
* @author Sebastian Strasser
* @date 2016-10-17
*/
/***********************************************/

#ifndef __GROOPSGUI__ADDVARIABLELDIALOG__
#define __GROOPSGUI__ADDVARIABLELDIALOG__

#include "base/importGroops.h"
#include "tree/treeElementGlobal.h"
#include <QDialog>
#include <QRegularExpression>

/***** TYPES ***********************************/

namespace Ui
{
  class AddVariableDialog;
}

class Tree;
class TreeElement;
class TreeElementGlobal;

/***** CLASS ***********************************/

class AddVariableDialog : public QDialog
{
  Q_OBJECT

  QString sourceType;
  bool    disableCreateLink;

public:
  AddVariableDialog(Tree *tree, const QString &type, const QString &label, bool disablePlace, bool disableCreateLink, QWidget *parent=nullptr);
 ~AddVariableDialog();

  QString type()     const;
  QString label()    const;
  Bool    inGlobal() const;
  Bool    setLink()  const;

private slots:
  void nameEdited(const QString &label);
  void typeChanged(int index);

private:
  Ui::AddVariableDialog *ui;
};

/***********************************************/

#endif
