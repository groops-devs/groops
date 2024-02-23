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
    tabEnvironment = nullptr;
    undoView  = nullptr;
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
    connect(ui->fileCloseAction,       SIGNAL(triggered(bool)), this, SLOT(fileClose()));
    connect(ui->fileCloseOtherAction,  SIGNAL(triggered(bool)), this, SLOT(fileCloseOther()));
    connect(ui->fileExitAction,        SIGNAL(triggered(bool)), this, SLOT(fileExit()));
    connect(ui->editFindReplaceAction, SIGNAL(triggered(bool)), this, SLOT(editFindReplace()));
    connect(ui->settingsCommandAction, SIGNAL(triggered(bool)), this, SLOT(settingsCommand()));
    connect(ui->settingsPathAction,    SIGNAL(triggered(bool)), this, SLOT(settingsPath()));
    connect(ui->settingsFontAction,    SIGNAL(triggered(bool)), this, SLOT(settingsFont()));
    connect(ui->helpAboutAction,       SIGNAL(triggered(bool)), this, SLOT(helpAbout()));

    // additional keyboard shortcuts
    // -----------------------------
    ui->editFindReplaceAction->setShortcuts({tr("Ctrl+F"), tr("Ctrl+H")});
    ui->editRemoveAction->setShortcuts({tr("Del"), tr("Ctrl+-")});

    // restore window
    // --------------
    setMinimumSize(QGuiApplication::primaryScreen()->size()/3);
    if(settings.contains("mainWindow/geometry"))
      restoreGeometry(settings.value("mainWindow/geometry").toByteArray());
    else
    {
      resize(minimumSizeHint());
      move(QStyle::alignedRect(Qt::LeftToRight, Qt::AlignCenter, size(), QGuiApplication::primaryScreen()->geometry()).topLeft());
    }
    restoreState(settings.value("mainWindow/state").toByteArray());
/*    QApplication::setFont(settings.value("misc/font", QApplication::font()).toString());*/
// qWarning()<<settings.value("misc/font", QApplication::font()).toString();

    // fill 'Open Recent' menu
    // ------------------------------------
    QStringList recentFileList = settings.value("recentFiles").toStringList();
    menuFileLastOpened = new QMenu(this);
    ui->fileOpenRecentAction->setEnabled(recentFileList.size()>0);
    connect(menuFileLastOpened,   SIGNAL(triggered(QAction*)), this, SLOT(fileLastOpened(QAction *)));
    for(int i=0; i<recentFileList.size(); i++)
      menuFileLastOpened->addAction(recentFileList[i])->setData(recentFileList[i]);
    ui->fileOpenRecentAction->setMenu(menuFileLastOpened);

    // ShowDescriptions status
    // -----------------------
    ui->helpShowDescriptionsAction->setCheckable(true);
    ui->helpShowResultsAction->setCheckable(true);

    // ActionList
    // ----------
    ActionList actionList;
    actionList.fileNewAction               = ui->fileNewAction;
    actionList.fileOpenAction              = ui->fileOpenAction;
    actionList.fileReOpenAction            = ui->fileReOpenAction;
    actionList.fileShowInManagerAction     = ui->fileShowInManagerAction;
    actionList.fileCloseAction             = ui->fileCloseAction;
    actionList.fileCloseOtherAction        = ui->fileCloseOtherAction;
    actionList.fileSaveAction              = ui->fileSaveAction;
    actionList.fileSaveAsAction            = ui->fileSaveAsAction;
    actionList.fileRunAction               = ui->fileRunAction;
    actionList.editCutAction               = ui->editCutAction;
    actionList.editCopyAction              = ui->editCopyAction;
    actionList.editPasteAction             = ui->editPasteAction;
    actionList.editPasteOverwriteAction    = ui->editPasteOverwriteAction;
    actionList.editAddAction               = ui->editAddAction;
    actionList.editRemoveAction            = ui->editRemoveAction;
    actionList.editAddVariableAction       = ui->editAddVariableAction;
    actionList.editSetLoopAction           = ui->editSetLoopAction;
    actionList.editSetConditionAction      = ui->editSetConditionAction;
    actionList.editEnabledAction           = ui->editEnabledAction;
    actionList.editEnableAllAction         = ui->editEnableAllAction;
    actionList.editDisableAllAction        = ui->editDisableAllAction;
    actionList.editRenameAction            = ui->editRenameAction;
    actionList.editUpdateNameAction        = ui->editUpdateNameAction;
    actionList.editAddCommentAction        = ui->editAddCommentAction;
    actionList.editCollapseAllAction       = ui->editCollapseAllAction;
    actionList.editOpenExternallyAction    = ui->editOpenExternallyAction;
    actionList.settingsPathAction          = ui->settingsPathAction;
    actionList.helpShowDescriptionsAction  = ui->helpShowDescriptionsAction;
    actionList.helpShowResultsAction       = ui->helpShowResultsAction;
    actionList.helpOpenDocumentationAction = ui->helpOpenDocumentationAction;

    // TabEnvironment
    // --------------
    show();
    tabEnvironment = new TabEnvironment(this, &actionList, undoGroup);
    connect(tabEnvironment, SIGNAL(fileChanged(const QString&, const QString&, bool)), this, SLOT(fileChanged(const QString&, const QString&, bool)));

    // Side bar and side bar widgets
    // -----------------------------
    // Open Files widget
    OpenFilesTreeWidget *openFilesTreeWidget = new OpenFilesTreeWidget(tabEnvironment, actionList);

    // Undo Stack widget
    undoView = new QUndoView(undoGroup);
    undoView->setAlternatingRowColors(true);

    // Program List widget
    programListWidget = new ProgramListWidget();

    // Add widgets to side bar
    sideBar = new SideBar(this);
    sideBar->addSideBarWidget("Open Files",   openFilesTreeWidget);
    sideBar->addSideBarWidget("Program List", programListWidget);
    sideBar->addSideBarWidget("Undo Stack",   undoView);

    // Basic layout
    // ------------
    // Horizontal splitter (side bar widgets | tab environment)
    QSplitter *splitter = new QSplitter;
    splitter->setOrientation(Qt::Horizontal);
    splitter->addWidget(sideBar->stackedWidget());
    splitter->addWidget(tabEnvironment);
    splitter->setChildrenCollapsible(false);
    splitter->setStretchFactor(0, 0);
    splitter->setStretchFactor(1, 1);
    splitter->setSizes({sideBar->lastWidth()});
    connect(splitter, SIGNAL(splitterMoved(int, int)), sideBar, SLOT(widthChanged(int, int)));

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
    tabEnvironment->setFocus();

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
  return tabEnvironment->currentTree();
}

/***********************************************/

bool MainWindow::fileNew()                            {return tabEnvironment->fileNew();}
bool MainWindow::fileOpen(const QString &fileName)    {return tabEnvironment->fileOpen(fileName);}
bool MainWindow::fileReOpen()                         {return tabEnvironment->fileReOpen();}
void MainWindow::fileLastOpened(QAction *whichAction) {tabEnvironment->fileOpen(whichAction->data().toString());}
bool MainWindow::fileClose()                          {return tabEnvironment->fileClose();}
void MainWindow::fileCloseOther()                     {tabEnvironment->fileCloseOther();}
void MainWindow::fileExit()                           {close();}

/***********************************************/

void MainWindow::editFindReplace()
{
  QString oldFindText = findReplaceDock->getFindText();

  auto widget = QApplication::focusWidget();
  QLineEdit *lineEdit = dynamic_cast<QLineEdit*>(widget);
  if(!lineEdit && dynamic_cast<QComboBox*>(widget))
    lineEdit = dynamic_cast<QComboBox*>(widget)->lineEdit();
  if(lineEdit && !lineEdit->selectedText().isEmpty())
    findReplaceDock->setFindText(lineEdit->selectedText());

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
  settings.setValue("misc/font", font);
  QApplication::setFont(font);
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

void MainWindow::fileChanged(const QString &caption, const QString &fileName, bool isClean)
{
  try
  {
    Tree *tree = tabEnvironment->currentTree();
    if(tree == nullptr)
      return;

    // adapt window title
    // ------------------
    QString capt = caption;
    if(!isClean)
      capt += tr(" [changed]");
    capt += tr(" - GROOPS");
    this->setWindowTitle(capt);

    // update Program List widget with programs from current schema
    // ----------------------------------------------------------------
    disconnect(programListWidget, SIGNAL(addProgram(const QString &)), nullptr, nullptr);
    programListWidget->init(tree);
    connect(programListWidget, SIGNAL(addProgram(const QString &)), tree, SLOT(addProgram(const QString &)));

    schemaSelector->setCurrentTreeSchema(tree->fileNameSchema());

    // update recentFiles list
    // -----------------------
    if(!fileName.isEmpty() && isClean) // !isClean cannot be a new file
    {
      QStringList fileNames = settings.value("recentFiles").toStringList();
      fileNames.removeOne(fileName); // is the name in the list? -> remove
      fileNames.prepend(fileName);
      while(fileNames.size() > 10)   // Maximum of 10 entries
        fileNames.removeLast();
      settings.setValue("recentFiles", fileNames);

      // update 'recently openend files' menu
      // ------------------------------------
      menuFileLastOpened->clear();
      for(auto &fileName : fileNames)
        menuFileLastOpened->addAction(fileName)->setData(fileName);
      menuFileLastOpened->setEnabled(fileNames.size());
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

void MainWindow::closeEvent(QCloseEvent *e)
{
  try
  {
    if(!tabEnvironment || tabEnvironment->okToAbandon())
    {
      settings.setValue("mainWindow/geometry", saveGeometry());
      settings.setValue("mainWindow/state",    saveState());
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
