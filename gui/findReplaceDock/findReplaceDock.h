/***********************************************/
/**
* @file findReplaceDock.h
*
* @brief Find/replace dock.
*
* @author Torsten Mayer-Guerr
* @author Sebastian Strasser
* @date 2013-10-30
*/
/***********************************************/

#ifndef __GROOPSGUI__FINDREPLACEDOCK__
#define __GROOPSGUI__FINDREPLACEDOCK__

#include <QDockWidget>
#include "base/importGroops.h"
#include "tree/tree.h"

/***** TYPES ***********************************/

namespace Ui
{
  class FindReplaceDock;
}

class MainWindow;

/***** CLASS ***********************************/

class FindReplaceDock : public QDockWidget
{
  Q_OBJECT

  Ui::FindReplaceDock *ui;
  MainWindow          *mainWindow;
  Tree                *tree;
  QSettings           *settings;

  void    messageNotFound();
  QRegularExpression getExpression();
  void    findNext(Bool forwards);

public:
  FindReplaceDock(MainWindow *parent);
 ~FindReplaceDock();

  void setFocusOnFind();
  void setFindText(QString text);
  QString getFindText() const;

public slots:
  void clickedNext();
  void clickedPrevious();
  void clickedReplace();
  void clickedReplaceAll();
  void editFindTextChanged(const QString &text);
  void searchTypeChanged();

protected:
  void closeEvent(QCloseEvent *e);
  void changeEvent(QEvent *e);
};

/***********************************************/

#endif
