/***********************************************/
/**
* @file mainWindow.h
*
* @brief The main window.
*
* @author Torsten Mayer-Guerr
* @author Sebastian Strasser
* @date 2011-09-02
*/
/***********************************************/

#ifndef __GROOPSGUI__MAINWINDOW_
#define __GROOPSGUI__MAINWINDOW_

#include <QMainWindow>
#include <QSettings>

/***** TYPES ***********************************/

namespace Ui
{
  class MainWindow;
}

class QMenu;
class QUndoGroup;
class QUndoView;
class TabEnvironment;
class FindReplaceDock;
class Tree;
class SideBar;
class ProgramListWidget;
class SchemaSelector;

/***** CLASS ***********************************/

class MainWindow : public QMainWindow
{
  Q_OBJECT

  Ui::MainWindow    *ui;
  QSettings          settings;
  QMenu             *menuFileLastOpened;
  QUndoGroup        *undoGroup;
  QUndoView         *undoView;
  TabEnvironment    *tabEnvironment;
  FindReplaceDock   *findReplaceDock;
  SideBar           *sideBar;
  ProgramListWidget *programListWidget;
  SchemaSelector    *schemaSelector;

public:
  MainWindow(QWidget *parent=nullptr);
 ~MainWindow();

  Tree *getCurrentTree() const;

public slots:
  bool fileNew();
  bool fileOpen(const QString &fileName=QString());
  bool fileReOpen();
  bool fileClose();
  void fileCloseOther();
  void fileExit();
  void fileLastOpened(QAction *whichAction);
  void editFindReplace();
  void settingsCommand();
  void settingsPath();
  void settingsFont();
  void helpAbout();
  void fileChanged(const QString &caption, const QString &fileName, bool isClean);

protected:
  void closeEvent(QCloseEvent *e);
  void changeEvent(QEvent *e);
};

/***** CLASS ***********************************/

class ActionList
{
public:
  QAction *fileNewAction;
  QAction *fileOpenAction;
  QAction *fileReOpenAction;
  QAction *fileOpenRecentAction;
  QAction *fileShowInManagerAction;
  QAction *fileSaveAction;
  QAction *fileSaveAsAction;
  QAction *fileRunAction;
  QAction *fileCloseAction;
  QAction *fileCloseOtherAction;
  QAction *editCutAction;
  QAction *editCopyAction;
  QAction *editPasteAction;
  QAction *editPasteOverwriteAction;
  QAction *editAddAction;
  QAction *editRemoveAction;
  QAction *editSetGlobalAction;
  QAction *editSetLoopAction;
  QAction *editRemoveLoopAction;
  QAction *editSetConditionAction;
  QAction *editRemoveConditionAction;
  QAction *editEnabledAction;
  QAction *editEnableAllAction;
  QAction *editDisableAllAction;
  QAction *editRenameAction;
  QAction *editUpdateNameAction;
  QAction *editAddCommentAction;
  QAction *editCollapseAllAction;
  QAction *editOpenExternallyAction;
  QAction *settingsPathAction;
  QAction *helpShowDescriptionsAction;
  QAction *helpShowResultsAction;
  QAction *helpOpenDocumentationAction;
};

/***********************************************/

#endif
