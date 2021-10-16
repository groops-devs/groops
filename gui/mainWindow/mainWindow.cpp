/***********************************************/
/**
* @file mainWindow.cpp
*
* @brief The main window.
*
* @author Torsten Mayer-Guerr
* @author Sebastian Strasser
* @date 2011-09-02
*
*/
/***********************************************/

#include "ui_mainWindow.h"
#include <QtDebug>
#include <QDesktopServices>
#include <QtWidgets>
#include <QFontDialog>
#include <QUndoGroup>
#include <QUndoView>
#include "base/importGroops.h"
#include "tree/tree.h"
#include "findReplaceDock/findReplaceDock.h"
#include "programDialog/programDialog.h"
#include "settingsDialog/settingsCommandDialog.h"
#include "settingsDialog/settingsPathDialog.h"
#include "mainWindow/schemaSelector.h"
#include "mainWindow/sideBar.h"
#include "mainWindow/tabs.h"
#include "mainWindow.h"

/***********************************************/

MainWindow::MainWindow(QWidget *parent) : QMainWindow(parent), ui(new Ui::MainWindow)
{
  try
  {
    ui->setupUi(this);
    ui->statusBar->setHidden(true);
    workspace = nullptr;
    undoView  = nullptr;
    settings  = new QSettings(this);
    undoGroup = new QUndoGroup(this);

    // Undo/Redo actions
    // -----------------
    QAction *editUndoAction = undoGroup->createUndoAction(this, tr("&Undo"));
    editUndoAction->setShortcuts(QKeySequence::Undo);
    editUndoAction->setIcon(QIcon(":/icons/scalable/edit-undo.svg"));
    ui->menuEdit->insertAction(ui->editCutAction, editUndoAction);
    ui->mainToolBar->addSeparator();
    ui->mainToolBar->addAction(editUndoAction);

    QAction *editRedoAction = undoGroup->createRedoAction(this, tr("&Redo"));
    editRedoAction->setShortcuts(QKeySequence::Redo);
    editRedoAction->setIcon(QIcon(":/icons/scalable/edit-redo.svg"));
    ui->menuEdit->insertAction(ui->editCutAction, editRedoAction);
    ui->mainToolBar->addAction(editRedoAction);
    ui->menuEdit->insertSeparator(ui->editCutAction);

    schemaSelector = new SchemaSelector(ui->mainToolBar);
    schemaSelector->setAction(ui->mainToolBar->addWidget(schemaSelector));

    // signal slot connections
    // -----------------------
    connect(ui->fileNewAction,         SIGNAL(triggered(bool)), this, SLOT(fileNew()));
    connect(ui->fileOpenAction,        SIGNAL(triggered(bool)), this, SLOT(fileOpen()));
    connect(ui->fileReOpenAction,      SIGNAL(triggered(bool)), this, SLOT(fileReOpen()));
    connect(ui->fileShowInManagerAction, SIGNAL(triggered(bool)), this, SLOT(fileShowInManager()));
    connect(ui->fileCloseAction,       SIGNAL(triggered(bool)), this, SLOT(fileClose()));
    connect(ui->fileCloseOtherAction,  SIGNAL(triggered(bool)), this, SLOT(fileCloseOther()));
    connect(ui->fileSaveAction,        SIGNAL(triggered(bool)), this, SLOT(fileSave()));
    connect(ui->fileSaveAsAction,      SIGNAL(triggered(bool)), this, SLOT(fileSaveAs()));
    connect(ui->fileExitAction,        SIGNAL(triggered(bool)), this, SLOT(fileExit()));
    connect(ui->fileRunAction,         SIGNAL(triggered(bool)), this, SLOT(fileRun()));
    connect(ui->editFindReplaceAction, SIGNAL(triggered(bool)), this, SLOT(editFindReplace()));
    connect(ui->settingsCommandAction, SIGNAL(triggered(bool)), this, SLOT(settingsCommand()));
    connect(ui->settingsPathAction,    SIGNAL(triggered(bool)), this, SLOT(settingsPath()));
    connect(ui->settingsFontAction,    SIGNAL(triggered(bool)), this, SLOT(settingsFont()));
    connect(ui->helpShowDescriptionsAction, SIGNAL(triggered(bool)), this, SLOT(showDescriptionsToggle(bool)));
    connect(ui->helpShowResultsAction, SIGNAL(triggered(bool)), this, SLOT(showResultsToggle(bool)));
    connect(ui->helpOpenDocumentation, SIGNAL(triggered(bool)), this, SLOT(openDocumentationExternally()));
    connect(ui->helpAboutAction,       SIGNAL(triggered(bool)), this, SLOT(helpAbout()));

    // additional keyboard shortcuts
    // -----------------------------
    ui->editFindReplaceAction->setShortcuts({tr("Ctrl+F"), tr("Ctrl+H")});

    // restore window
    // --------------
    setMinimumSize(QGuiApplication::primaryScreen()->size()/3);
    if(settings->contains("mainWindow/geometry"))
      restoreGeometry(settings->value("mainWindow/geometry").toByteArray());
    else {
      resize(minimumSizeHint());
      move(QStyle::alignedRect(Qt::LeftToRight, Qt::AlignCenter, size(), QGuiApplication::primaryScreen()->geometry()).topLeft());
    }
    restoreState(settings->value("mainWindow/state").toByteArray());
/*    QApplication::setFont(settings->value("misc/font", QApplication::font()).toString());*/
// qWarning()<<settings->value("misc/font", QApplication::font()).toString();

    // fill 'Open Recent' menu
    // ------------------------------------
    QStringList recentFileList = settings->value("recentFiles").toStringList();
    menuFileLastOpened = new QMenu(this);
    ui->fileOpenRecentAction->setEnabled(recentFileList.size()>0);
    connect(menuFileLastOpened,   SIGNAL(triggered(QAction*)), this, SLOT(fileLastOpened(QAction *)));
    for(int i=0; i<recentFileList.size(); i++)
      menuFileLastOpened->addAction(recentFileList[i])->setData(recentFileList[i]);
    ui->fileOpenRecentAction->setMenu(menuFileLastOpened);

    // ShowDescriptions status
    // -----------------------
    ui->helpShowDescriptionsAction->setCheckable(true);
    ui->helpShowDescriptionsAction->setChecked( settings->value("misc/showDescriptions", true).toBool() );
    ui->helpShowResultsAction->setCheckable(true);
    ui->helpShowResultsAction->setChecked( settings->value("misc/showResults", true).toBool() );

    // ActionList
    // ----------
    ActionList actionList;
    actionList.fileNewAction             = ui->fileNewAction;
    actionList.fileOpenAction            = ui->fileOpenAction;
    actionList.fileReOpenAction          = ui->fileReOpenAction;
    actionList.fileShowInManagerAction   = ui->fileShowInManagerAction;
    actionList.fileCloseAction           = ui->fileCloseAction;
    actionList.fileCloseOtherAction      = ui->fileCloseOtherAction;
    actionList.fileSaveAction            = ui->fileSaveAction;
    actionList.fileSaveAsAction          = ui->fileSaveAsAction;
    actionList.fileRunAction             = ui->fileRunAction;
    actionList.editCutAction             = ui->editCutAction;
    actionList.editCopyAction            = ui->editCopyAction;
    actionList.editPasteAction           = ui->editPasteAction;
    actionList.editPasteOverwriteAction  = ui->editPasteOverwriteAction;
    actionList.editAddAction             = ui->editAddAction;
    actionList.editRemoveAction          = ui->editRemoveAction;
    actionList.editSetGlobalAction       = ui->editSetGlobalAction;
    actionList.editSetLoopAction         = ui->editSetLoopAction;
    actionList.editRemoveLoopAction      = ui->editRemoveLoopAction;
    actionList.editSetConditionAction    = ui->editSetConditionAction;
    actionList.editRemoveConditionAction = ui->editRemoveConditionAction;
    actionList.editEnabledAction         = ui->editEnabledAction;
    actionList.editEnableAllAction       = ui->editEnableAllAction;
    actionList.editDisableAllAction      = ui->editDisableAllAction;
    actionList.editRenameAction          = ui->editRenameAction;
    actionList.editUpdateNameAction      = ui->editUpdateNameAction;
    actionList.editCommentAction         = ui->editCommentAction;
    actionList.editCollapseAllAction     = ui->editCollapseAllAction;
    actionList.openExternallyAction      = ui->openExternallyAction;

    // TabEnvironment
    // --------------
    show();
    workspace = new TabEnvironment(this, &actionList, undoGroup);
    workspace->setShowDescriptions(ui->helpShowDescriptionsAction->isChecked());
    workspace->setShowResults(ui->helpShowResultsAction->isChecked());
    connect(workspace, SIGNAL(fileChanged(const QString &, bool)), this,           SLOT(fileChanged(const QString &, bool)));
    connect(workspace, SIGNAL(fileOpened(const QString &)),        this,           SLOT(addToRecentFiles(const QString &)));
    connect(workspace, SIGNAL(schemaChanged()),                    schemaSelector, SLOT(updateList()));

    // Side bar and side bar widgets
    // -----------------------------

    // Open Files widget
    OpenFilesTreeWidget *openFilesTreeWidget = new OpenFilesTreeWidget();
    openFilesTreeWidget->init(workspace, &actionList);

    // Undo Stack widget
    undoView = new QUndoView(undoGroup);
    undoView->setAlternatingRowColors(true);

    // Program List widget
    programListWidget = new ProgramListWidget();

    // File System Browser widget
//    QTreeWidget *fileBrowserWidget = new QTreeWidget;
//    QTreeWidgetItem *treeWidgetItem = new QTreeWidgetItem;
//    treeWidgetItem->setText(0, "Test");
//    fileBrowserWidget->addTopLevelItem(treeWidgetItem);

    // Add widgets to side bar
    sideBar = new SideBar(this);
    sideBar->setSideBarWidgetsHidden(settings->value("sideBar/isHidden", false).toBool());
    sideBar->addSideBarWidget("Open Files", openFilesTreeWidget);
//    sideBar->addSideBarWidget("File System Browser", fileBrowserWidget);
    sideBar->addSideBarWidget("Program List", programListWidget);
    sideBar->addSideBarWidget("Undo Stack", undoView);

    // Basic layout
    // ------------

    // Horizontal splitter (side bar widgets | tab environment)
    QSplitter *splitter = new QSplitter;
    splitter->setOrientation(Qt::Horizontal);
    splitter->addWidget(sideBar->stackedWidget());
    splitter->addWidget(workspace);
    splitter->setChildrenCollapsible(false);
    splitter->setStretchFactor(0, 0);
    splitter->setStretchFactor(1, 1);
    QList<int> widths;
    widths.push_back(sideBar->stackedWidget()->lastWidth());
    splitter->setSizes(widths);
    connect(splitter, SIGNAL(splitterMoved(int, int)), sideBar->stackedWidget(), SLOT(widthChanged(int, int)));

    // Set layout (side bar | splitter)
    QHBoxLayout *layout = new QHBoxLayout;
    layout->addWidget(sideBar);
    layout->addWidget(splitter);
    layout->setAlignment(sideBar, Qt::AlignTop);
    layout->setContentsMargins(0,0,0,0);
    layout->setSpacing(0);

    // Set layout as central widget
    QWidget *window = new QWidget;
    window->setLayout(layout);
    setCentralWidget(window);
    workspace->setFocus();

    // Find/Replace Dock
    // -----------------
    findReplaceDock = new FindReplaceDock(this);
    findReplaceDock->setAllowedAreas(Qt::TopDockWidgetArea | Qt::BottomDockWidgetArea);
    addDockWidget(Qt::BottomDockWidgetArea, findReplaceDock);
    findReplaceDock->hide();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

MainWindow::~MainWindow()
{
  delete ui;
  delete undoView;
}

/***********************************************/

Tree *MainWindow::getCurrentTree() const
{
  return workspace->currentTree();
}

/***********************************************/

void MainWindow::fileNew()                           {workspace->newFile();}
void MainWindow::fileOpen(const QString &fileName)   {workspace->openFile(fileName);}
void MainWindow::fileOpenInitial()                   {workspace->openInitialFiles();}
void MainWindow::fileReOpen()                        {workspace->reopenFile();}
void MainWindow::fileShowInManager()                 {workspace->showFileInManager();}
void MainWindow::fileSave()                          {workspace->saveFile();}
void MainWindow::fileSaveAs()                        {workspace->saveAsFile();}
void MainWindow::fileRun()                           {workspace->execFile();}
void MainWindow::fileClose()                         {workspace->closeFile();}
void MainWindow::fileCloseOther()                    {workspace->closeOtherFiles();}
void MainWindow::fileExit()                          {close();}

/***********************************************/

void MainWindow::fileLastOpened(QAction *whichAction)
{
  workspace->openFile(whichAction->data().toString());
}

/***********************************************/

void MainWindow::editFindReplace()
{
  QString oldFindText = findReplaceDock->getFindText();
  QComboBox* comboBox = dynamic_cast<QComboBox*>(QApplication::focusWidget());
  if(comboBox && comboBox->lineEdit() && !comboBox->lineEdit()->selectedText().isEmpty())
    findReplaceDock->setFindText(comboBox->lineEdit()->selectedText());
  else
  {
    QLineEdit* lineEdit = dynamic_cast<QLineEdit*>(QApplication::focusWidget());
    if(lineEdit && !lineEdit->selectedText().isEmpty())
      findReplaceDock->setFindText(lineEdit->selectedText());
  }

  if(!findReplaceDock->isVisible())
    findReplaceDock->setVisible(true);
  else if(findReplaceDock->getFindText() == oldFindText)
    findReplaceDock->setVisible(false);

  if(findReplaceDock->isVisible())
    findReplaceDock->setFocusOnFind();
}

/***********************************************/

void MainWindow::settingsCommand()
{
  SettingsCommandDialog dialog(this);
  dialog.exec();
}

/***********************************************/

void MainWindow::settingsPath()
{
  SettingsPathDialog dialog(this);
  if(dialog.exec() == QDialog::Accepted)
    schemaSelector->updateList();
}

/***********************************************/

void MainWindow::settingsFont()
{
  bool ok;
  QFont font = QFontDialog::getFont(&ok, QApplication::font());
  if(!ok)
    return;
  settings->setValue("misc/font", font);
  QApplication::setFont(font);
}

/***********************************************/

void MainWindow::showDescriptionsToggle(bool state)
{
  settings->setValue("misc/showDescriptions", state);
  workspace->setShowDescriptions(state);
}


/***********************************************/

void MainWindow::showResultsToggle(bool state)
{
  settings->setValue("misc/showResults", state);
  workspace->setShowResults(state);
}

/***********************************************/

void MainWindow::helpAbout()
{
  QMessageBox::about(this, tr("About GROOPS GUI"),
                     tr("<h2>GROOPS Graphical User Interface</h2>"
                        "<p>Gravity Recovery Object Oriented Programming System"
                        "<p>Torsten Mayer-G&uuml;rr"
                        "<br>Wolfgang Mayer-G&uuml;rr"
                        "<br>Sebastian Strasser"
                        "<p>GitHub repository: <a href=\"https://github.com/groops-devs/groops\">https://github.com/groops-devs/groops</a></p>"
                        "<p>This GUI application uses the <a href=\"https://www.qt.io\">Qt framework</a> and Google's <a href=\"https://github.com/google/material-design-icons\">Material design icons</a>.</p>"));
}

/***********************************************/
/***********************************************/

void MainWindow::fileChanged(const QString &fileName, bool changed)
{
  try
  {
    // adapt window title
    // ------------------
    QString capt = QFileInfo(fileName).fileName();
    if(fileName.isEmpty())
      capt = workspace->currentTabText();
    if(changed)
      capt += tr(" [changed]");
    capt += tr(" - GROOPS");
    this->setWindowTitle(capt);

    // update Program List widget with programs from current schema
    // ----------------------------------------------------------------
    Tree *tree = workspace->currentTree();
    if(tree == nullptr)
      return;
    if(programListWidget)
    {
      disconnect(programListWidget, SIGNAL(programSelected(int)), nullptr, nullptr);
      programListWidget->init(tree);
      connect(programListWidget, SIGNAL(programSelected(int)), tree, SLOT(addProgram(int)));
    }

    schemaSelector->setCurrentTreeSchema(tree->schemaFile());

    ui->fileShowInManagerAction->setEnabled(!tree->fileName().isEmpty());
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

void MainWindow::addToRecentFiles(const QString &fileName)
{
  try
  {
    if(!fileName.startsWith("/")) // exclude "new*" tabs
      return;

    // update recentFiles list
    // -----------------------
    QStringList recentFileList = settings->value("recentFiles").toStringList();
    recentFileList.prepend( QFileInfo(fileName).absoluteFilePath() );
    // is the name in the list? -> remove
    for(int i=1; i<recentFileList.size(); i++)
      if(QFileInfo(fileName).absoluteFilePath() == recentFileList[i])
        recentFileList.removeAt(i--);
    // Maximum of 10 entriesD
    while(recentFileList.size()>10)
      recentFileList.removeLast();
    settings->setValue("recentFiles", recentFileList);

    // update 'recently openend files' menu
    // ------------------------------------
    menuFileLastOpened->clear();
    for(int i=0; i<recentFileList.size(); i++)
      menuFileLastOpened->addAction(recentFileList[i])->setData(recentFileList[i]);
    menuFileLastOpened->setEnabled(recentFileList.size()>0);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

void MainWindow::openDocumentationExternally()
{
  QUrl url = workspace->documentationUrl();
  if(url.isEmpty())
    return;

  if(!QFileInfo(url.toLocalFile()).exists())
  {
    url = QUrl::fromLocalFile(url.path().replace(url.fileName(), "index.html"));
    if(!QFileInfo(url.toLocalFile()).exists())
      return;
  }

  QDesktopServices::openUrl(url);
}

/***********************************************/

void MainWindow::closeEvent(QCloseEvent *e)
{
  try
  {
    if((!workspace) || workspace->okToAbandon())
    {
      settings->setValue("mainWindow/geometry", saveGeometry());
      settings->setValue("mainWindow/state",    saveState());
      settings->setValue("sideBar/width",       sideBar->stackedWidget()->lastWidth());
      settings->setValue("sideBar/isHidden",    sideBar->stackedWidget()->isHidden());
      e->accept();
    }
    else
      e->ignore();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

void MainWindow::changeEvent(QEvent *e)
{
  QMainWindow::changeEvent(e);
  switch(e->type())
  {
    case QEvent::LanguageChange:
      ui->retranslateUi(this);
      break;
    default:
      break;
  }
}

/***********************************************/
