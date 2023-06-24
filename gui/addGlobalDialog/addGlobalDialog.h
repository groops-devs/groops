/***********************************************/
/**
* @file addGlobalDialog.h
*
* @brief Dialog to add global variable.
*
* @author Sebastian Strasser
* @date 2016-10-17
*/
/***********************************************/

#ifndef __GROOPSGUI__ADDGLOBALDIALOG__
#define __GROOPSGUI__ADDGLOBALDIALOG__

#include "base/importGroops.h"
#include "tree/treeElementGlobal.h"
#include <QDialog>
#include <QRegularExpression>

/***** TYPES ***********************************/

namespace Ui
{
  class AddGlobalDialog;
}

class Tree;
class TreeElement;
class TreeElementGlobal;

/***** CLASS ***********************************/

class AddGlobalDialog : public QDialog
{
  Q_OBJECT

  std::vector<XsdElementPtr> xsdElements;
  QStringList                globalElementNames;

public:
  AddGlobalDialog(TreeElementGlobal *globalRoot, QWidget *parent=nullptr);
  ~AddGlobalDialog();

  QString       elementName();
  XsdElementPtr elementType();

private slots:
  void nameEdited(const QString &name);

private:
  Ui::AddGlobalDialog *ui;
};

/***********************************************/

#endif
