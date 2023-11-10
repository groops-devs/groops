/***********************************************/
/**
* @file settingsCommandDialog.h
*
* @brief Dialog to set up commands.
*
* @author Torsten Mayer-Guerr
* @author Sebastian Strasser
* @date 2006-10-13
*/
/***********************************************/

#ifndef __GROOPSGUI__SETTINGSCOMMANDDIALOG__
#define __GROOPSGUI__SETTINGSCOMMANDDIALOG__

#include <QDialog>
#include <QSettings>
#include "base/importGroops.h"

/***** TYPES ***********************************/

namespace Ui
{
  class SettingsCommandDialog;
}

class QAbstractButton;

/***** CLASS ***********************************/

class SettingsCommandDialog : public QDialog
{
  Q_OBJECT

  Ui::SettingsCommandDialog *ui;
  QSettings                  settings;

  void insertItem(int row, const QString &label, const QString &command);

public:
  SettingsCommandDialog(QWidget *parent=nullptr);
 ~SettingsCommandDialog();

 static void readCommandList(QStringList &lableList, QStringList &commandList);

public slots:
  void clickedUp();
  void clickedDown();
  void clickedApply();
  void clickedStdBox(QAbstractButton *button);

private slots:
  void tableCellChanged(int row, int column);

protected:
  void changeEvent(QEvent *e);
};

/***********************************************/

#endif
