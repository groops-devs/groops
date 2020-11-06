/***********************************************/
/**
* @file tabs.cpp
*
* @brief Tab environment.
*
* @author Torsten Mayer-Guerr
* @author Sebastian Strasser
* @date 2011-07-16
*/
/***********************************************/

#include <QApplication>
#include <QSettings>
#include <QTabWidget>
#include <QMenu>
#include <QContextMenuEvent>
#include <QHeaderView>
#include <QPushButton>
#include <QLabel>
#include <QVBoxLayout>
#include <QUrl>
#include <QUndoGroup>
#include <QDrag>
#include <QProcess>
#include <QDebug>
#include <QFileDialog>
#include <QMimeData>
#include <QDesktopServices>
#include "base/importGroops.h"
#include "tree/tree.h"
#include "tree/treeItem.h"
#include "tree/treeElement.h"
#include "mainWindow/mainWindow.h"
#include "tabs.h"

/***********************************************/

TabEnvironment::TabEnvironment(QWidget *parent, ActionList *actionList, QUndoGroup *undoGroup) : QTabWidget(parent)
{
  try
  {
    TabBar *tabBar = new TabBar(this, this);
    setTabBar(tabBar);
    setMovable(true);
    setTabsClosable(true);
    setDocumentMode(true);
    this->actionList = *actionList;
    this->undoGroup  = undoGroup;
    this->settings   = new QSettings(this);
    this->newFileCounter = 0;
    this->isFullyLoaded = false;

    QPushButton *addButton = new QPushButton(QIcon(":/icons/scalable/document-new.svg"), "");
    addButton->setFlat(true);
    addButton->resize(tabBar->contentsRect().height(), tabBar->contentsRect().height());
    addButton->setIconSize(addButton->contentsRect().size()*0.9);
    setCornerWidget(addButton, Qt::TopRightCorner);
    connect(addButton, SIGNAL(clicked()),              this, SLOT(newFile()));

    connect(this,      SIGNAL(currentChanged(int)),    this, SLOT(tabChanged(int)));
    connect(this,      SIGNAL(tabCloseRequested(int)), this, SLOT(tabClose(int)));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

TabEnvironment::~TabEnvironment()
{
  Tree *tree = currentTree();
  if(tree)
  {
    const Double width = this->width();
    settings->setValue("tree/relativeColumnWidth0", tree->columnWidth(0)/width);
    settings->setValue("tree/relativeColumnWidth1", tree->columnWidth(1)/width);
    settings->setValue("tree/relativeColumnWidth2", tree->columnWidth(2)/width);
  }
  delete settings;
}

/***********************************************/
/***********************************************/

void TabEnvironment::saveOpenFiles()
{
  try
  {
    QStringList fileList;
    for(int i=0; i<count(); i++)
    {
      Tree *tree = treeAt(i);
      if(tree)
      {
        QString fileName = tree->fileName();
        if(!fileName.isEmpty())
          fileList<<fileName;
      }
    }
    settings->setValue("openFiles", fileList);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

Tree *TabEnvironment::currentTree() const
{
  TabPage *tabPage = dynamic_cast<TabPage*>(currentWidget());
  return (tabPage ? tabPage->tree() : nullptr);
}

/***********************************************/

Tree *TabEnvironment::treeAt(int index) const
{
  TabPage *tabPage = dynamic_cast<TabPage*>(widget(index));
  return (tabPage ? tabPage->tree() : nullptr);
}

/***********************************************/

QStringList TabEnvironment::allFileNames() const
{
  QStringList fileNames;
  for(int i=0; i<count(); i++)
  {
    Tree *tree = treeAt(i);
    if(tree && !tree->fileName().isEmpty())
      fileNames.push_back(tree->fileName());
    else if(tree && tree->fileName().isEmpty())
      fileNames.push_back(tabText(i));
  }

  return fileNames;
}

/***********************************************/

QStringList TabEnvironment::openedFileNames() const
{
  QStringList fileNames;
  for(int i=0; i<count(); i++)
  {
    Tree *tree = treeAt(i);
    if(tree && !tree->fileName().isEmpty())
      fileNames.push_back(tree->fileName());
  }

  return fileNames;
}

/***********************************************/

QStringList TabEnvironment::newFileNames() const
{
  QStringList fileNames;
  for(int i=0; i<count(); i++)
  {
    Tree *tree = treeAt(i);
    if(tree && tree->fileName().isEmpty())
      fileNames.push_back(tabText(i));
  }

  return fileNames;
}

/***********************************************/

QString TabEnvironment::currentTabText() const
{
  return tabText(currentIndex());
}

/***********************************************/

void TabEnvironment::resizeTreeColumns(const std::vector<int> &columnWidths)
{
  for(int i=0; i<count(); i++)
  {
    Tree *tree = treeAt(i);
    if(tree)
      for(size_t i = 0; i < columnWidths.size(); i++)
        tree->setColumnWidth(i, columnWidths.at(i));
  }
}

/***********************************************/

void TabEnvironment::setShowDescriptions(Bool state)
{
  for(int i=0; i<count(); i++)
  {
    Tree *tree = treeAt(i);
    if(tree)
      tree->setShowDescriptions(state);
  }
}

/***********************************************/

void TabEnvironment::setShowResults(Bool state)
{
  for(int i=0; i<count(); i++)
  {
    Tree *tree = treeAt(i);
    if(tree)
      tree->setShowResults(state);
  }
}

/***********************************************/

Bool TabEnvironment::newFile()
{
  try
  {
    Tree *tree = new Tree(nullptr, &actionList, this);
    undoGroup->addStack(tree->undoStack());
    connect(tree, SIGNAL(fileChanged(const QString &, bool)), this, SLOT(treeChanged(const QString &, bool)));
    TabPage *page = new TabPage(tree, this);
    insertTab(currentIndex()+1, page, QString("new%1").arg(++newFileCounter).toUtf8());
    tree->header()->blockSignals(true);
    Tree *currentTree = this->currentTree();
    if(currentTree)
    {
      tree->setColumnWidth(0, currentTree->columnWidth(0));
      tree->setColumnWidth(1, currentTree->columnWidth(1));
      tree->setColumnWidth(2, currentTree->columnWidth(2));
    }
    setCurrentWidget(page);
    tree->header()->blockSignals(false);
    emit fileChanged(currentTabText(), false);
    return true;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

Bool TabEnvironment::closeFile()
{
  try
  {
    Tree *tree = currentTree();
    int index = currentIndex();
    if(!tree || !tree->okToAbandon())
      return false;
    if(count()==1)
      newFile();
    undoGroup->removeStack(tree->undoStack());
    tree->deleteLater();
    removeTab(index);
    return true;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

Bool TabEnvironment::closeOtherFiles()
{
  try
  {
    bool allOthersClosed = true;
    Tree* remainingTree = currentTree();
    for(int i = count(); i --> 0; )
    {
      setCurrentIndex(i);
      Tree* tree = currentTree();
      if(!tree || tree == remainingTree)
        continue;

      if(tree->okToAbandon())
      {
        undoGroup->removeStack(tree->undoStack());
        tree->deleteLater();
        removeTab(i);
      }
      else
        allOthersClosed = false;
    }

    remainingTree = currentTree();
    if(remainingTree)
      emit fileChanged(remainingTree->fileName(), false);
    return allOthersClosed;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

Bool TabEnvironment::reopenFile()
{
  Tree *tree = currentTree();
  bool success = tree ? tree->reopenFile() : false;
  if(success)
    emit fileOpened(tree->fileName());
  return success;
}

/***********************************************/

void TabEnvironment::showFileInManager() const
{
  Tree *tree = currentTree();
  if(tree)
    QDesktopServices::openUrl(QUrl::fromLocalFile(QFileInfo(tree->fileName()).absoluteDir().path()));
}

/***********************************************/

Bool TabEnvironment::openFile(const QString &fileName, bool emitFileOpened)
{
  if(!fileName.isEmpty() && openedFileNames().contains(fileName))
  {
    setCurrentIndex(allFileNames().indexOf(fileName));
    Tree *tree = currentTree();
    bool success = tree ? tree->openFile(fileName) : false;
    if(success && emitFileOpened)
      emit fileOpened(fileName);
    return success;
  }
  else if(!fileName.isEmpty() && QFileInfo(fileName).isFile())
  {
    if(newFile())
    {
      Tree *tree = currentTree();
      bool success = tree ? tree->openFile(fileName) : false;
      if(success && emitFileOpened)
        emit fileOpened(fileName);
      else if(!success)
        closeFile();
      return success;
    }
  }
  else
  {
    QString name = fileName;
    if(name.isEmpty())
    {
      Tree *tree = currentTree();
      if(tree)
        name = QFileInfo(tree->fileName()).absolutePath();
    }
    if(name.isEmpty())
      name = settings->value("files/workingDirectory", QString("../scenario")).toString();

    name = QFileDialog::getOpenFileName(this, tr("Open File - GROOPS"), name, tr("XML files (*.xml)"));
    if(!name.isEmpty())
    {
      bool success = openFile(name);
      if(success && emitFileOpened)
        emit fileOpened(fileName);
      return success;
    }
  }

  return false;
}

void TabEnvironment::openInitialFiles()
{
  QStringList fileList = settings->value("openFiles").toStringList();
  if(fileList.size())
  {
    for(int i=0; i<fileList.size(); i++)
      if(QFileInfo(fileList[i]).isFile())
        openFile(fileList[i], false/*emitFileOpened*/);
  }
  else
    newFile();

  if(fileList.size())
    this->newFileCounter = 0;
  else
    this->newFileCounter = 1;
}

/***********************************************/

Bool TabEnvironment::saveFile()
{
  Tree *tree = currentTree();
  return tree ? tree->saveFile() : false;
}

/***********************************************/

Bool TabEnvironment::saveAsFile(const QString &fileName)
{
  Tree *tree = currentTree();
  return tree ? tree->saveAsFile(fileName) : false;
}

/***********************************************/

Bool TabEnvironment::execFile()
{
  Tree *tree = currentTree();
  return tree ? tree->execFile() : false;
}

/***********************************************/

Bool TabEnvironment::okToAbandon()
{
  for(int i=0; i<count(); i++)
  {
    setCurrentIndex(i);
    Tree *tree = treeAt(i);
    if(tree && !tree->okToAbandon())
      return false;
  }
  saveOpenFiles();
  return true;
}

/***********************************************/

QUrl TabEnvironment::documentationUrl()
{
  Tree *tree = currentTree();
  if(tree)
  {
    TreeItem *item = tree->getSelectedItem();
    if(item)
    {
      TreeElement *element = item->treeElement();
      QString name = (element->isProgram() ? element->selectedValue() : element->type())+".html";
      return QUrl::fromLocalFile(settings->value("files/documentationDirectory").toString()+"/"+name);
    }
  }
  return QUrl::fromLocalFile(settings->value("files/documentationDirectory").toString()+"/index.html");
}

/***********************************************/
/***********************************************/

void TabEnvironment::showEvent(QShowEvent *event)
{
  if(!isFullyLoaded)
  {
    for(int i = 0; i < count(); i++)
    {
      Tree *tree = treeAt(i);
      if(tree)
      {
        tree->setColumnWidth(0, static_cast<int>(std::round(settings->value("tree/relativeColumnWidth0", 0.2).toDouble()*width())));
        tree->setColumnWidth(1, static_cast<int>(std::round(settings->value("tree/relativeColumnWidth1", 0.4).toDouble()*width())));
        tree->setColumnWidth(2, static_cast<int>(std::round(settings->value("tree/relativeColumnWidth2", 0.2).toDouble()*width())));
      }
    }
    isFullyLoaded = true;
  }

  QTabWidget::showEvent(event);
}

/***********************************************/

void TabEnvironment::tabClose(int index)
{
  // deselect all other tabs
  for(int i=0; i<count(); i++)
  {
    Tree *tree = treeAt(i);
    if(tree && i != index)
      tree->setSelectedItem(nullptr);
  }
  setCurrentIndex(index);
  closeFile();
}

/***********************************************/

void TabEnvironment::tabChanged(int index)
{
  if(tabBar()->currentIndex() != index)
    tabBar()->setCurrentIndex(index);

  // deselect all other tabs
  for(int i=0; i<count(); i++)
  {
    Tree *tree = treeAt(i);
    if(tree && i != index)
      tree->setSelectedItem(nullptr);
  }

  Tree *tree = currentTree();
  if(tree)
  {
    undoGroup->setActiveStack(tree->undoStack());
    emit fileChanged(tree->fileName(), tree->isChanged());
  }
}

/***********************************************/

void TabEnvironment::treeChanged(const QString &/*fileName*/, bool /*changed*/)
{
  for(int i=0; i<count(); i++)
  {
    Tree *tree = treeAt(i);
    if(tree)
    {
      QString fileName = tree->fileName();
      setTabToolTip(i, fileName);
      setTabText   (i, (fileName.isEmpty()) ? tabText(i) : QFileInfo(fileName).fileName());
      setTabIcon   (i, ((tree->isChanged()) ? QIcon(":/icons/scalable/document-save.svg") : QIcon()));
    }
  }

  Tree *tree = currentTree();
  if(tree)
  {
    undoGroup->setActiveStack(tree->undoStack());
    emit fileChanged(tree->fileName(), tree->isChanged());
  }
}

/***********************************************/
/***********************************************/

TabBar::TabBar(TabEnvironment *tabEnvironment, QWidget *parent) : QTabBar(parent)
{
  tabs = tabEnvironment;
  setAcceptDrops(true);
}

/***********************************************/

TabBar::~TabBar()
{
}

/***********************************************/

void TabBar::mouseDoubleClickEvent(QMouseEvent */*e*/)
{
  tabs->newFile();
}

/***********************************************/

void TabBar::contextMenuEvent(QContextMenuEvent *e)
{
  int index = tabAt(e->pos());
  if(index < 0)
    return;
  tabs->setCurrentIndex(index);

  QMenu *contextMenu = new QMenu(this);
  contextMenu->addAction(tabs->actionList.fileNewAction);
  contextMenu->addAction(tabs->actionList.fileOpenAction);
  contextMenu->addAction(tabs->actionList.fileReOpenAction);
  contextMenu->addAction(tabs->actionList.fileShowInManagerAction);
  contextMenu->addSeparator();
  contextMenu->addAction(tabs->actionList.fileSaveAction);
  contextMenu->addAction(tabs->actionList.fileSaveAsAction);
  contextMenu->addSeparator();
  contextMenu->addAction(tabs->actionList.fileRunAction);
  contextMenu->addSeparator();
  contextMenu->addAction(tabs->actionList.editEnableAllAction);
  contextMenu->addAction(tabs->actionList.editDisableAllAction);
  contextMenu->addAction(tabs->actionList.editCollapseAllAction);
  contextMenu->addSeparator();
  contextMenu->addAction(tabs->actionList.fileCloseAction);
  contextMenu->addAction(tabs->actionList.fileCloseOtherAction);
  contextMenu->exec(e->globalPos());
  delete contextMenu;
}

/***********************************************/

void TabBar::dragEnterEvent(QDragEnterEvent *event)
{
  try
  {
    // file names?
    if(event->mimeData()->hasUrls())
      event->acceptProposedAction();
    else
      QTabBar::dragEnterEvent(event);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

void TabBar::dragMoveEvent(QDragMoveEvent *event)
{
  try
  {
    if(this->geometry().contains(event->pos()))
    {
      // hacky way to abort QDrag (may not work on Windows)
      QKeyEvent esc(QEvent::KeyPress, Qt::Key_Escape, Qt::NoModifier);
      QApplication::sendEvent(this, &esc);
      return;
    }

    // file names?
    if(event->mimeData()->hasUrls())
      event->acceptProposedAction();
    else
      QTabBar::dragMoveEvent(event);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

void TabBar::dropEvent(QDropEvent *event)
{
  try
  {
    // file names?
    if(!event->mimeData()->hasUrls())
    {
      QTabBar::dropEvent(event);
      return;
    }

    QList<QUrl> urls = event->mimeData()->urls();
    for(int i=0; i<urls.size(); i++)
    {
      tabs->openFile(urls.at(i).toLocalFile());
    }
    event->acceptProposedAction();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

void TabBar::mousePressEvent(QMouseEvent *event)
{
    if (event->button() == Qt::LeftButton)
        dragStartPosition = event->pos();

    QTabBar::mousePressEvent(event);
}

/***********************************************/

void TabBar::mouseMoveEvent(QMouseEvent *event)
{
  if(!(event->buttons() & Qt::LeftButton))
    return;

  // if cursor is within TabBar, call original event handler
  if(this->geometry().contains(event->pos()))
  {
    QTabBar::mouseMoveEvent(event);
    return;
  }

  if((event->pos() - dragStartPosition).manhattanLength() < QApplication::startDragDistance())
    return;

  int index = tabAt(dragStartPosition);
  if(index < 0)
    return;

  // Tab is dragged outside of TabBar
  // --------------------------------
  tabs->setCurrentIndex(index);
  Tree *tree = tabs->currentTree();
  if(!tree)
    return;

  QDrag *drag = new QDrag(this);
  QMimeData *mimeData = new QMimeData;

  // add file to MIME data
  QList<QUrl> urls;
  if(!tree->fileName().isEmpty())
    urls.push_back(QUrl::fromLocalFile(tree->fileName()));
  mimeData->setUrls(urls);
  drag->setMimeData(mimeData);

  // add small picture to mouse cursor
  QPixmap pixmap(tree->size());
  tree->render(&pixmap);
  drag->setPixmap(pixmap.scaled(200,200, Qt::KeepAspectRatio));
  drag->setHotSpot(QPoint(-20,-20));

  Qt::DropAction dropAction = drag->exec(Qt::CopyAction | Qt::MoveAction);

  if(dropAction != Qt::IgnoreAction)
  {
    QWidget *target = qobject_cast<QWidget*>(drag->target());

    if(target && !tree->fileName().isEmpty())
    {
      QProcess *process = new QProcess(this);
      QString program = QFileInfo(QCoreApplication::applicationFilePath()).fileName();
      process->startDetached(program, {"-c", "-f", tree->fileName()});
    }

    // TabBar refresh after drag requires mouse button release
    QMouseEvent click(QEvent::MouseButtonRelease, dragStartPosition, Qt::LeftButton, Qt::NoButton, Qt::NoModifier);
    QApplication::sendEvent(this, &click);

    // close tab if file is moved (Shift+Drag)
    Qt::KeyboardModifiers modifiers = QApplication::keyboardModifiers();
    if(modifiers == Qt::ShiftModifier)
      tabs->closeFile();
  }
}

/***********************************************/

TabPage::TabPage(Tree *tree, TabEnvironment *parent) : QWidget(parent)
{
  this->_tree = tree;

  // Bar handling external file changes
  {
    QPushButton *buttonReopen = new QPushButton(QIcon(":/icons/scalable/view-refresh.svg"), "Reopen", this);
    QPushButton *buttonIgnore = new QPushButton(QIcon(":/icons/scalable/ignore.svg"), "Ignore", this);
    QLabel *iconLabel = new QLabel(this);
    iconLabel->setPixmap(QIcon(":/icons/scalable/warning.svg").pixmap(24,24));
    QHBoxLayout *layoutBar = new QHBoxLayout(this);
    layoutBar->addWidget(iconLabel);
    layoutBar->addWidget(new QLabel("File was modified externally. Reopen?", this), 1);
    layoutBar->addWidget(buttonReopen);
    layoutBar->addWidget(buttonIgnore);
    layoutBar->setMargin(3);
    barFileChanged = new QFrame(this);
    barFileChanged->setFrameStyle(QFrame::Box);
    barFileChanged->setLayout(layoutBar);
    barFileChanged->setVisible(false);
    const QString highlightColor = barFileChanged->palette().highlight().color().name().right(6);
    barFileChanged->setStyleSheet(".QFrame { color: #"+highlightColor+"; background-color: #4d"+highlightColor+" }");

    connect(buttonReopen, &QPushButton::clicked,  this, &TabPage::fileChangeClickedReopen);
    connect(buttonIgnore, &QPushButton::clicked,  this, &TabPage::barClickedIgnore);
    connect(tree,         &Tree::treeFileChanged, this, &TabPage::treeFileChanged);
  }

  // Bar handling unknown elements
  {
    QPushButton *buttonShowAll = new QPushButton(QIcon(":/icons/scalable/edit-find-replace.svg"), "Show all", this);
    QPushButton *buttonRemoveAll = new QPushButton(QIcon(":/icons/scalable/edit-delete.svg"), "Remove all", this);
    QPushButton *buttonIgnore = new QPushButton(QIcon(":/icons/scalable/ignore.svg"), "Ignore", this);
    buttonRemoveAll->setMinimumWidth(95);
    QLabel *iconLabel = new QLabel(this);
    iconLabel->setPixmap(QIcon(":/icons/scalable/help-about.svg").pixmap(24,24));
    QHBoxLayout *layoutBar = new QHBoxLayout(this);
    labelUnknownElements = new QLabel(this);
    layoutBar->addWidget(iconLabel);
    layoutBar->addWidget(labelUnknownElements, 1);
    layoutBar->addWidget(buttonShowAll);
    layoutBar->addWidget(buttonRemoveAll);
    layoutBar->addWidget(buttonIgnore);
    layoutBar->setMargin(3);
    barUnknownElements = new QFrame(this);
    barUnknownElements->setFrameStyle(QFrame::Box);
    barUnknownElements->setLayout(layoutBar);
    barUnknownElements->setVisible(false);
    const QString highlightColor = barUnknownElements->palette().highlight().color().name().right(6);
    barUnknownElements->setStyleSheet(".QFrame { color: #"+highlightColor+"; background-color: #4d"+highlightColor+" }");

    connect(buttonShowAll,   &QPushButton::clicked,         tree, &Tree::expandUnknownElements);
    connect(buttonRemoveAll, &QPushButton::clicked,         tree, &Tree::removeAllUnknownElements);
    connect(buttonIgnore,    &QPushButton::clicked,         this, &TabPage::barClickedIgnore);
    connect(tree,            &Tree::unknownElementsChanged, this, &TabPage::unknownElementsChanged);
  }

  // Bar handling renamed elements
  {
    QPushButton *buttonShowAll = new QPushButton(QIcon(":/icons/scalable/edit-find-replace.svg"), "Show all", this);
    QPushButton *buttonUpdateAll = new QPushButton(QIcon(":/icons/scalable/edit-rename.svg"), "Update all", this);
    QPushButton *buttonIgnore = new QPushButton(QIcon(":/icons/scalable/ignore.svg"), "Ignore", this);
    buttonUpdateAll->setMinimumWidth(95);
    QLabel *iconLabel = new QLabel(this);
    iconLabel->setPixmap(QIcon(":/icons/scalable/help-about.svg").pixmap(24,24));
    QHBoxLayout *layoutBar = new QHBoxLayout(this);
    labelRenamedElements = new QLabel(this);
    layoutBar->addWidget(iconLabel);
    layoutBar->addWidget(labelRenamedElements, 1);
    layoutBar->addWidget(buttonShowAll);
    layoutBar->addWidget(buttonUpdateAll);
    layoutBar->addWidget(buttonIgnore);
    layoutBar->setMargin(3);
    barRenamedElements = new QFrame(this);
    barRenamedElements->setFrameStyle(QFrame::Box);
    barRenamedElements->setLayout(layoutBar);
    barRenamedElements->setVisible(false);
    const QString highlightColor = barRenamedElements->palette().highlight().color().name().right(6);
    barRenamedElements->setStyleSheet(".QFrame { color: #"+highlightColor+"; background-color: #4d"+highlightColor+" }");

    connect(buttonShowAll,   &QPushButton::clicked,         tree, &Tree::expandRenamedElements);
    connect(buttonUpdateAll, &QPushButton::clicked,         tree, &Tree::updateAllRenamedElements);
    connect(buttonIgnore,    &QPushButton::clicked,         this, &TabPage::barClickedIgnore);
    connect(tree,            &Tree::renamedElementsChanged, this, &TabPage::renamedElementsChanged);
  }

  QVBoxLayout *layout = new QVBoxLayout();
  layout->setSpacing(3);
  layout->setContentsMargins(0, 3, 0, 0);
  layout->addWidget(barFileChanged);
  layout->addWidget(barUnknownElements);
  layout->addWidget(barRenamedElements);
  layout->addWidget(tree, 1);
  setLayout(layout);
}

/***********************************************/

void TabPage::treeFileChanged(const QString &/*fileName*/, bool changed)
{
  if(barFileChanged)
    barFileChanged->setVisible(changed);
}

/***********************************************/

void TabPage::unknownElementsChanged(int count)
{
  if(barUnknownElements)
    barUnknownElements->setVisible(count>0);
  if(labelUnknownElements)
    labelUnknownElements->setText(QString("File contains %1 unknown elements.").arg(count));
}

/***********************************************/

void TabPage::renamedElementsChanged(int count)
{
  if(barRenamedElements)
    barRenamedElements->setVisible(count>0);
  if(labelRenamedElements)
    labelRenamedElements->setText(QString("File contains %1 elements that were renamed in the schema.").arg(count));
}

/***********************************************/

void TabPage::fileChangeClickedReopen()
{
  if(tree() && barFileChanged)
    if(tree()->reopenFile())
      barFileChanged->setHidden(true);
}

/***********************************************/

void TabPage::barClickedIgnore()
{
  QFrame *bar = dynamic_cast<QFrame*>(sender()->parent());
  if(bar)
    bar->setHidden(true);
}

/***********************************************/
