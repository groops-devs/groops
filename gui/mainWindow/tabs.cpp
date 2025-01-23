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
#include "tree/treeElementProgram.h"
#include "mainWindow/mainWindow.h"
#include "tabs.h"

/***********************************************/

TabEnvironment::TabEnvironment(QWidget *parent, ActionList *actionList, QUndoGroup *undoGroup) : QTabWidget(parent)
{
  try
  {
    TabBar *tabBar = new TabBar(this);
    setTabBar(tabBar);
    setMovable(true);
    setTabsClosable(true);
    setDocumentMode(true);
    this->actionList     = *actionList;
    this->undoGroup      = undoGroup;
    this->newFileCounter = 0;
    this->isFullyLoaded  = false;

    QPushButton *addButton = new QPushButton(QIcon(":/icons/scalable/document-new.svg"), "");
    addButton->setFlat(true);
    addButton->resize(tabBar->contentsRect().height(), tabBar->contentsRect().height());
    addButton->setIconSize(addButton->contentsRect().size()*0.9);
    setCornerWidget(addButton, Qt::TopRightCorner);

    connect(addButton, SIGNAL(clicked()),              this, SLOT(fileNew()));
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
    settings.setValue("tree/relativeColumnWidth0", tree->columnWidth(0)/width);
    settings.setValue("tree/relativeColumnWidth1", tree->columnWidth(1)/width);
    settings.setValue("tree/relativeColumnWidth2", tree->columnWidth(2)/width);
  }
}

/***********************************************/
/***********************************************/

Tree *TabEnvironment::currentTree() const
{
  return dynamic_cast<Tree*>(currentWidget());
}

/***********************************************/

void TabEnvironment::setCurrentTree(Tree *tree)
{
  setCurrentWidget(tree);
}

/***********************************************/

QVector<Tree*> TabEnvironment::trees() const
{
  QVector<Tree*> trees(count());
  for(int i=0; i<count(); i++)
    trees[i] = dynamic_cast<Tree*>(widget(i));
  return trees;
}

/***********************************************/

Tree *TabEnvironment::treeAt(int index) const
{
  return dynamic_cast<Tree*>(widget(index));
}

/***********************************************/
/***********************************************/

bool TabEnvironment::newTab(const QString &fileName)
{
  try
  {
    Tree *tree = new Tree(nullptr, &actionList, this);

    bool ok = false;
    if(fileName.isEmpty())
      ok = tree->fileNew(QString("new%1").arg(++newFileCounter));
    else
      ok = tree->fileOpen(fileName);
    if(!ok)
    {
      delete tree;
      return false;
    }

    // copy sections sizes from current tree
    if(currentTree())
    {
      tree->setColumnWidth(0, currentTree()->columnWidth(0));
      tree->setColumnWidth(1, currentTree()->columnWidth(1));
      tree->setColumnWidth(2, currentTree()->columnWidth(2));
    }

    undoGroup->addStack(tree->undoStack);
    connect(tree, SIGNAL(fileChanged(const QString&, const QString&, bool)), this, SLOT(treeFileChanged(const QString&, const QString&, bool)));
    connect(tree, SIGNAL(sectionResized(int, int, int)), this, SLOT(treeSectionResized(int, int, int)));
    setCurrentIndex(insertTab(currentIndex()+1, tree, tree->caption()));
    return true;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

bool TabEnvironment::fileOpen(const QString &fileNameConst)
{
  try
  {
    QStringList fileNames;
    QString     fileName = fileNameConst;

    // no file name -> open file selector
    if(fileName.isEmpty())
    {
      if(currentTree())
        fileName = QFileInfo(currentTree()->fileName()).absolutePath();
      if(fileName.isEmpty())
        fileName = settings.value("files/workingDirectory", QString("../scenario")).toString();
      fileNames = QFileDialog::getOpenFileNames(this, tr("Open File - GROOPS"), fileName, tr("XML files (*.xml)"));
      if(fileNames.isEmpty())
        return false;
    }
    else
      fileNames << fileName;

    for(auto fileName : fileNames)
    {
      // file is already opened? -> select tab and reOpen
      bool found = false;
      for(int i=0; i<count(); i++)
      {
        if(treeAt(i)->fileName() == QFileInfo(fileName).absoluteFilePath())
        {
          found = true;
          setCurrentIndex(i);
          if(!fileReOpen())
            return false;
          break;
        }
      }
      if(found)
        continue;

      // current tab is new and unchanged -> overwrite
      if(currentTree() && currentTree()->fileName().isEmpty() && currentTree()->isClean())
      {
        if(!currentTree()->fileOpen(fileName))
          return false;
      }
      else if(!newTab(fileName))
          return false;

    }
    return true;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

bool TabEnvironment::fileReOpen()
{
  return currentTree() ? currentTree()->fileOpen(currentTree()->fileName()) : false;
}

/***********************************************/

bool TabEnvironment::fileClose()
{
  try
  {
    Tree *tree = currentTree();
    if(!tree || !tree->okToAbandon())
      return false;
    int index = currentIndex();
    if(count()==1)
      fileNew();
    removeTab(index);
    undoGroup->removeStack(tree->undoStack);
    tree->deleteLater();
    return true;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

void TabEnvironment::fileCloseOther()
{
  try
  {
    Tree *remainingTree = currentTree();
    if(remainingTree)
      for(int i=count(); i-->0;)
      {
        setCurrentIndex(i);
        if(currentTree() != remainingTree)
          fileClose();
      }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

bool TabEnvironment::okToAbandon()
{
  QStringList fileList;
  for(int i=0; i<count(); i++)
  {
    setCurrentIndex(i);
    if(!treeAt(i)->fileName().isEmpty())
      fileList<<treeAt(i)->fileName();
    if(!treeAt(i)->okToAbandon())
      return false;
  }

  settings.setValue("openFiles", fileList);
  return true;
}

/***********************************************/
/***********************************************/

void TabEnvironment::showEvent(QShowEvent *event)
{
  if(!isFullyLoaded)
    for(Tree *tree : trees())
    {
      tree->setColumnWidth(0, static_cast<int>(std::round(settings.value("tree/relativeColumnWidth0", 0.2).toDouble()*width())));
      tree->setColumnWidth(1, static_cast<int>(std::round(settings.value("tree/relativeColumnWidth1", 0.4).toDouble()*width())));
      tree->setColumnWidth(2, static_cast<int>(std::round(settings.value("tree/relativeColumnWidth2", 0.2).toDouble()*width())));
    }
  isFullyLoaded = true;
  QTabWidget::showEvent(event);
}

/***********************************************/

void TabEnvironment::resizeEvent(QResizeEvent *event)
{
  Tree *tree = currentTree();
  if(tree && (event->oldSize().width() > 0))
  {
    const Double scale = 1.*event->size().width()/event->oldSize().width();
    tree->setColumnWidth(0, static_cast<int>(std::round(scale*tree->columnWidth(0))));
    tree->setColumnWidth(1, static_cast<int>(std::round(scale*tree->columnWidth(1))));
    tree->setColumnWidth(2, static_cast<int>(std::round(scale*tree->columnWidth(2))));
  }
  QTabWidget::resizeEvent(event);
}

/***********************************************/
/***********************************************/

void TabEnvironment::treeSectionResized(int logicalIndex, int /*oldSize*/, int newSize)
{
  if(logicalIndex > 2)
    return;
  for(Tree *tree : trees())
  {
    tree->blockSignals(true);
    tree->setColumnWidth(logicalIndex, newSize);
    tree->blockSignals(false);
  }
}

/***********************************************/

void TabEnvironment::treeFileChanged(const QString &/*caption*/, const QString &/*fileName*/, bool /*isClean*/)
{
  for(int i=0; i<count(); i++)
  {
    setTabToolTip(i, treeAt(i)->fileName());
    setTabText   (i, treeAt(i)->caption());
    setTabIcon   (i, treeAt(i)->isClean() ? QIcon() : QIcon(":/icons/scalable/document-save.svg"));
  }

  if(currentTree())
    emit fileChanged(currentTree()->caption(), currentTree()->fileName(), currentTree()->isClean());
}

/***********************************************/

void TabEnvironment::tabClose(int index)
{
  // deselect all other tabs
  for(int i=0; i<count(); i++)
    if(i != index)
      treeAt(i)->setSelectedItem(nullptr);
  setCurrentIndex(index);
  fileClose();
}

/***********************************************/

void TabEnvironment::tabChanged(int index)
{
  if(tabBar()->currentIndex() != index)
    tabBar()->setCurrentIndex(index);

  // deselect all other tabs
  for(int i=0; i<count(); i++)
    if(i != index)
      treeAt(i)->setCurrent(false);

  if(currentTree())
  {
    currentTree()->setCurrent(true);
    undoGroup->setActiveStack(currentTree()->undoStack);
    emit fileChanged(currentTree()->caption(), currentTree()->fileName(), currentTree()->isClean());
  }
}

/***********************************************/
/***********************************************/

TabBar::TabBar(TabEnvironment *tabEnvironment) : QTabBar(tabEnvironment), tabEnvironment(tabEnvironment)
{
  setAcceptDrops(true);
}

/***********************************************/

void TabBar::mouseDoubleClickEvent(QMouseEvent */*e*/)
{
  tabEnvironment->fileNew();
}

/***********************************************/

void TabBar::contextMenuEvent(QContextMenuEvent *e)
{
  int index = tabAt(e->pos());
  if(index < 0)
    return;
  tabEnvironment->setCurrentIndex(index);

  QMenu *contextMenu = new QMenu(this);
  contextMenu->addAction(tabEnvironment->actionList.fileNewAction);
  contextMenu->addAction(tabEnvironment->actionList.fileOpenAction);
  contextMenu->addAction(tabEnvironment->actionList.fileReOpenAction);
  contextMenu->addAction(tabEnvironment->actionList.fileShowInManagerAction);
  contextMenu->addSeparator();
  contextMenu->addAction(tabEnvironment->actionList.fileSaveAction);
  contextMenu->addAction(tabEnvironment->actionList.fileSaveAsAction);
  contextMenu->addSeparator();
  contextMenu->addAction(tabEnvironment->actionList.fileRunAction);
  contextMenu->addSeparator();
  contextMenu->addAction(tabEnvironment->actionList.editEnableAllAction);
  contextMenu->addAction(tabEnvironment->actionList.editDisableAllAction);
  contextMenu->addAction(tabEnvironment->actionList.editCollapseAllAction);
  contextMenu->addSeparator();
  contextMenu->addAction(tabEnvironment->actionList.fileCloseAction);
  contextMenu->addAction(tabEnvironment->actionList.fileCloseOtherAction);
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
  if(event->mimeData()->hasUrls()) // file names?
  {
    for(auto &url : event->mimeData()->urls())
      tabEnvironment->fileOpen(url.toLocalFile());
    event->acceptProposedAction();
  }
  else
    QTabBar::dropEvent(event);
}

/***********************************************/

void TabBar::mousePressEvent(QMouseEvent *event)
{
  if(event->button() == Qt::LeftButton)
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
  tabEnvironment->setCurrentIndex(index);
  Tree *tree = tabEnvironment->currentTree();
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
      tabEnvironment->fileClose();
  }
}

/***********************************************/
