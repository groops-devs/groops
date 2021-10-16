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

/***** TYPES ***********************************/

namespace Ui
{
  class MainWindow;
}

class QSettings;
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

  Ui::MainWindow      *ui;
  QSettings           *settings;
  QMenu               *menuFileLastOpened;
  QUndoGroup          *undoGroup;
  QUndoView           *undoView;
  TabEnvironment      *workspace;
  FindReplaceDock     *findReplaceDock;
  SideBar             *sideBar;
  ProgramListWidget   *programListWidget;
  SchemaSelector      *schemaSelector;

public:
  MainWindow(QWidget *parent=nullptr);
 ~MainWindow();

  Tree *getCurrentTree() const;

public slots:
  void fileNew();
  void fileOpen(const QString &fileName = "");
  void fileOpenInitial();
  void fileReOpen();
  void fileShowInManager();
  void fileSave();
  void fileSaveAs();
  void fileRun();
  void fileClose();
  void fileCloseOther();
  void fileExit();
  void fileLastOpened(QAction *whichAction);
  void editFindReplace();
  void settingsCommand();
  void settingsPath();
  void settingsFont();
  void showDescriptionsToggle(bool state);
  void showResultsToggle(bool state);
  void helpAbout();
  void fileChanged(const QString &fileName, bool changed);
  void addToRecentFiles(const QString &fileName);
  void openDocumentationExternally();

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
  QAction *fileCloseAction;
  QAction *fileCloseOtherAction;
  QAction *fileSaveAction;
  QAction *fileSaveAsAction;
  QAction *fileRunAction;
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
  QAction *editCommentAction;
  QAction *editCollapseAllAction;
  QAction *openExternallyAction;
};

/***********************************************/

#endif
