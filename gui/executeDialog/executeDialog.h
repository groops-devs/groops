/***********************************************/
/**
* @file executeDialog.h
*
* @brief Execute dialog.
*
* @author Torsten Mayer-Guerr
* @author Sebastian Strasser
* @date 2006-10-13
*/
/***********************************************/

#ifndef __GROOPSGUI__EXECUTEDIALOG__
#define __GROOPSGUI__EXECUTEDIALOG__

#include <QDialog>
#include <QList>
#include "base/importGroops.h"

/***** TYPES ***********************************/

namespace Ui
{
  class ExecuteDialog;
}

class QSettings;
class QTreeWidgetItem;
class QAbstractButton;
class Tree;
class TreeElement;

/***** CLASS ***********************************/

class ExecuteDialog : public QDialog
{
  Q_OBJECT

  Ui::ExecuteDialog      *ui;
  QSettings              *settings;
  Tree                   *tree;
  QList<TreeElement*>     programList;
  QList<QTreeWidgetItem*> itemList;

public:
  ExecuteDialog(Tree *tree, QWidget *parent=nullptr);
 ~ExecuteDialog();

public slots:
  void editCommand();
  void toggleLogfile();
  void browseLogfile();
  void clickedAll();
  void clickedNone();
  void clickedApply();
  void clickedStdBox(QAbstractButton *button);

protected:
  void changeEvent(QEvent *e);

private:
  void addChildren(QTreeWidgetItem *parentItem, TreeElement *parentElement);
};

/***********************************************/

#endif
