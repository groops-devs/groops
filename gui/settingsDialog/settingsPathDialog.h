/***********************************************/
/**
* @file settingsPathDialog.h
*
* @brief Dialog to set up default paths and files.
*
* @author Torsten Mayer-Guerr
* @author Sebastian Strasser
* @date 2006-10-13
*/
/***********************************************/

#ifndef __GROOPSGUI__SETTINGSPATHDIALOG__
#define __GROOPSGUI__SETTINGSPATHDIALOG__

#include <QDialog>
#include "base/importGroops.h"

/***** TYPES ***********************************/

namespace Ui
{
  class SettingsPathDialog;
}

class QSettings;

/***** CLASS ***********************************/

class SettingsPathDialog : public QDialog
{
  Q_OBJECT

  Ui::SettingsPathDialog *ui;
  QSettings              *settings;

public:
  SettingsPathDialog(QWidget *parent=nullptr);
 ~SettingsPathDialog();

public slots:
  void clickedSchemaAdd();
  void clickedSchemaRemove();
  void clickedTemplateFile();
  void clickedWorkingDir();
  void clickedDocumentationDir();
  void clickedOk();
  void schemaListChanged();

protected:
  void changeEvent(QEvent *e);
};

/***********************************************/

#endif
