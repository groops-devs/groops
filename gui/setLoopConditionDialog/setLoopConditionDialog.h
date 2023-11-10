/***********************************************/
/**
* @file setLoopConditionDialog.h
*
* @brief Set loop/condition for an unbounded element.
*
* @author Sebastian Strasser
* @date 2017-02-09
*/
/***********************************************/

#ifndef __GROOPSGUI__SETLOOPCONDITIONDIALOG__
#define __GROOPSGUI__SETLOOPCONDITIONDIALOG__

#include "base/importGroops.h"
#include "tree/treeElementGlobal.h"
#include <QDialog>

/***** TYPES ***********************************/

namespace Ui
{
  class SetLoopConditionDialog;
}

class Tree;
class TreeElement;
class TreeElementGlobal;

/***** CLASS ***********************************/

class SetLoopConditionDialog : public QDialog
{
  Q_OBJECT

  QStringList globalElementNames;

public:
  SetLoopConditionDialog(TreeElementGlobal *globalRoot, const QString &type, QWidget *parent=nullptr);
 ~SetLoopConditionDialog();

  QString name();

private slots:
  void radioButtonToggled(bool checked);
  void nameEdited(const QString &name);

private:
  Ui::SetLoopConditionDialog *ui;
};

/***********************************************/

#endif
