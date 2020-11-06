/***********************************************/
/**
* @file programDialog.cpp
*
* @brief Program List Widget and Dialog.
*
* @author Torsten Mayer-Guerr
* @author Sebastian Strasser
* @date 2015-05-15
*/
/***********************************************/

#include "ui_programDialog.h"
#include <QSettings>
#include <QTreeWidget>
#include <QTreeWidgetItem>
#include <QList>
#include <QLabel>
#include <QDragEnterEvent>
#include <QMimeData>
#include <QDrag>
#include <QDialogButtonBox>
#include "base/importGroops.h"
#include "tree/tree.h"
#include "tree/treeElementProgram.h"
#include "programDialog.h"

/***********************************************/

ProgramListWidget::ProgramListWidget(QWidget *parent) : QWidget(parent)
{
  try
  {
    treeWidget               = new TreeWidget(this);
    lineEdit                 = new QLineEdit(this);
    QLabel      *label       = new QLabel("Filter:", this);
    QGridLayout *gridLayout  = new QGridLayout(this);

    gridLayout->addWidget(treeWidget, 0, 0, 1, 2);
    gridLayout->addWidget(label,      1, 0);
    gridLayout->addWidget(lineEdit,   1, 1);
    gridLayout->setContentsMargins(0, 0, 0, gridLayout->spacing());

    this->settings = new QSettings(this);

    // init table
    // ----------
    QStringList headerLabel;
    treeWidget->setHeaderLabels(headerLabel << "Program" << "Comment" << "Tags");
    treeWidget->hideColumn(2);
    treeWidget->setColumnWidth(0, settings->value("programListWidget/columnWidth", 200).toInt());
    treeWidget->setAlternatingRowColors(true);
    treeWidget->setRootIsDecorated(false);
    treeWidget->setItemsExpandable(false);
    treeWidget->setSelectionMode(QAbstractItemView::SingleSelection);
    treeWidget->setDragEnabled(true);
    treeWidget->setAcceptDrops(false);
    treeWidget->setDragDropMode(QAbstractItemView::DragOnly);

    connect(treeWidget, SIGNAL(itemSelectionChanged()),               this, SLOT(itemSelectionChanged()));
    connect(treeWidget, SIGNAL(itemActivated(QTreeWidgetItem*, int)), this, SLOT(itemActivated(QTreeWidgetItem*, int)));
    connect(lineEdit,   SIGNAL(textEdited(QString)),                  this, SLOT(filterChanged(QString)));

    lineEdit->setFocus();

    lineEdit->installEventFilter(this);
    treeWidget->installEventFilter(this);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

ProgramListWidget::ProgramListWidget(Tree *tree, int indexSelected, QWidget *parent) : ProgramListWidget(parent)
{
  init(tree, indexSelected);
}

/***********************************************/

ProgramListWidget::~ProgramListWidget()
{
  settings->setValue("programListWidget/columnWidth", columnWidth());

  delete settings;
}

/***********************************************/

void ProgramListWidget::init(Tree *tree, int indexSelected)
{
  try
  {
    if(tree)
      this->tree = tree;
    else
      throw(Exception("tree is null pointer"));

    // clear list
    treeWidget->clear();

    // fill list
    programList = tree->programListFromSchema();
    QTreeWidgetItem *selectedItem = nullptr;
    for(size_t i = 0; i < programList.size(); i++)
    {
      QTreeWidgetItem *item = new QTreeWidgetItem();
      item->setFlags(Qt::ItemIsEnabled | Qt::ItemIsSelectable);
      item->setIcon(0, QIcon(":/icons/scalable/program.svg"));
      item->setText(0, programList.at(i)->name);
      item->setText(1, programList.at(i)->annotation);
      item->setText(2, programList.at(i)->tags.join(", "));

      if(!programList.at(i)->tags.size()) // old schema without tags
        programList.at(i)->tags.push_back("Untagged");
      QList<QTreeWidgetItem*> tagItems = treeWidget->findItems(programList.at(i)->tags.at(0), Qt::MatchExactly, 0);
      QTreeWidgetItem *tagItem;
      if(tagItems.size())
        tagItem = tagItems.at(0);
      else
      {
        tagItem = new QTreeWidgetItem();
        tagItem->setText(0, programList.at(i)->tags.at(0));
        treeWidget->addTopLevelItem(tagItem);
      }
      tagItem->addChild(item);

      if(static_cast<int>(i) == indexSelected)
        selectedItem = item;
    }

    if(selectedItem)
    {
      selectedItem->setSelected(true);
      treeWidget->setCurrentItem(selectedItem);
    }

    treeWidget->expandAll();

    filterChanged(lineEdit->text());
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

void ProgramListWidget::setColumnWidth(int width)
{
  try
  {
    if(width < 0 || !treeWidget)
      return;

    treeWidget->setColumnWidth(0, width);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

int ProgramListWidget::selectedIndex()
{
  try
  {
    return treeWidget->getItemIndex(treeWidget->currentItem());
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

void ProgramListWidget::setUnknownProgram(const QString &name)
{
  try
  {
    QTreeWidgetItem *tagItem = new QTreeWidgetItem();
    tagItem->setText(0, "Unknown");
    treeWidget->addTopLevelItem(tagItem);

    QTreeWidgetItem *item = new QTreeWidgetItem();
    item->setFlags(Qt::ItemIsEnabled | Qt::ItemIsSelectable);
    item->setIcon(0, QIcon(":/icons/scalable/element-unknown.svg"));
    item->setText(0, name);
    tagItem->addChild(item);
    treeWidget->expandAll();

    if(selectedIndex() < 0)
    {
      item->setSelected(true);
      treeWidget->setCurrentItem(item);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}


/***********************************************/

bool ProgramListWidget::eventFilter(QObject *obj, QEvent *event)
{
  try
  {
    if(event->type() == QEvent::KeyPress)
    {
      QKeyEvent *keyEvent = static_cast<QKeyEvent*>(event);
      const bool isFilterKey = (!keyEvent->text().isEmpty() || keyEvent->key() == Qt::Key_Left   || keyEvent->key() == Qt::Key_Right
                                                            || keyEvent->key() == Qt::Key_Delete || keyEvent->key() == Qt::Key_Backspace);
      if(obj != treeWidget && !isFilterKey)
      {
        treeWidget->setFocus();
        return QApplication::sendEvent(treeWidget, keyEvent);
      }
      else if(obj != lineEdit && isFilterKey)
      {
        lineEdit->setFocus();
        return QApplication::sendEvent(lineEdit, keyEvent);
      }
    }

    return QWidget::eventFilter(obj, event);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

void ProgramListWidget::itemActivated(QTreeWidgetItem *item, int /*column*/)
{
  try
  {
    if(!item || item->childCount() || item->isHidden())
      return;

    emit programSelected(treeWidget->getItemIndex(item));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

void ProgramListWidget::filterChanged(QString text)
{
  try
  {
    QStringList filterStrings = text.split(" ");
    if(!filterStrings.size())
      return;

    // find items containing the first filter string
    QList<QTreeWidgetItem*> itemsFound;
    itemsFound.append(treeWidget->findItems(filterStrings.at(0), Qt::MatchContains | Qt::MatchRecursive, 0));
    itemsFound.append(treeWidget->findItems(filterStrings.at(0), Qt::MatchContains | Qt::MatchRecursive, 1));
    itemsFound.append(treeWidget->findItems(filterStrings.at(0), Qt::MatchContains | Qt::MatchRecursive, 2));
    itemsFound = itemsFound.toSet().values();

    // remove items that don't contain further filter strings
    for(int i = 1; i < filterStrings.size(); i++)
      for(int j = itemsFound.size()-1; j >= 0; j--)
        if(!itemsFound.at(j)->text(0).contains(filterStrings.at(i), Qt::CaseInsensitive) && !itemsFound.at(j)->text(1).contains(filterStrings.at(i), Qt::CaseInsensitive))
          itemsFound.removeAt(j);

    // hide all items
    QTreeWidgetItemIterator it(treeWidget);
    while(*it) {
      (*it)->setHidden(true);
      ++it;
    }

    // show found items
    for(int i = 0; i < itemsFound.size(); i++)
    {
      itemsFound.at(i)->setHidden(false);
      if(itemsFound.at(i)->parent() != nullptr)
        itemsFound.at(i)->parent()->setHidden(false);
    }

    itemSelectionChanged();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

void ProgramListWidget::itemSelectionChanged()
{
  try
  {
    QTreeWidgetItem *item = treeWidget->currentItem();
    emit selectionChanged(item && !item->childCount() && !item->isHidden());
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/
/***********************************************/

ProgramListWidget::TreeWidget::TreeWidget(ProgramListWidget *parent) : QTreeWidget(parent)
{
  try
  {
    if(parent)
      programListWidget = parent;
    else
      throw(Exception("parent is null pointer"));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}



/***********************************************/

int ProgramListWidget::TreeWidget::getItemIndex(QTreeWidgetItem *item)
{
  try
  {
    if(item == nullptr)
      return -1;

    QTreeWidgetItemIterator it(this, QTreeWidgetItemIterator::NoChildren);
    int index = 0;
    while(*it) {
      if((*it) == item)
        return index;
      ++it;
      ++index;
    }

    return -1;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

void ProgramListWidget::TreeWidget::mousePressEvent(QMouseEvent *event)
{
  try
  {
    if(event->button() == Qt::LeftButton)
      dragStartPosition = event->pos();
    QTreeWidget::mousePressEvent(event);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

void ProgramListWidget::TreeWidget::mouseMoveEvent(QMouseEvent *event)
{
  try
  {
    // use original event handler if it is not a movement with pressed left mouse button starting over an item
    QTreeWidgetItem *item = itemAt(event->pos());
    if((!(event->buttons() & Qt::LeftButton)) || ((event->pos()-dragStartPosition).manhattanLength() < QApplication::startDragDistance()) ||
       !item || (item && item->childCount()))
    {
      QTreeWidget::mouseMoveEvent(event);
      return;
    }

    // create xml data
    QString xmlData;
    int index = getItemIndex(item);
    if(index >= 0 && index < static_cast<int>(programListWidget->programList.size()))
    {
      TreeElement *element = TreeElementProgram::newTreeElement(programListWidget->tree, nullptr, programListWidget->programList.at(index), "", XmlNodePtr(nullptr), false);
      if(element->name().startsWith("==="))
        return;
      QString type = element->tree->programType();
      if(type.isEmpty())
        return;
      QString name = type.left(type.lastIndexOf("Type"));
      XmlNodePtr xmlNode = XmlNodePtr(new XmlNode(name));
      writeAttribute(xmlNode, "xsdType", type);
      xmlNode->addChild(element->getXML());

      QTextStream stream(&xmlData, QIODevice::WriteOnly);
      XmlNode::write(stream, xmlNode);
    }

    // create mime data
    QMimeData *mimeData = new QMimeData;
    mimeData->setData("application/x-groops", xmlData.toUtf8());
    if(mimeData == nullptr)
    {
      QTreeWidget::mouseMoveEvent(event);
      return;
    }

    // perform drag
    QDrag *drag = new QDrag(viewport());
    drag->setMimeData(mimeData);
    drag->exec(Qt::MoveAction|Qt::CopyAction);

    event->accept();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/
/***********************************************/

ProgramDialog::ProgramDialog(TreeElementProgram *treeElement, QWidget *parent) : QDialog(parent), ui(new Ui::ProgramDialog)
{
  try
  {
    ui->setupUi(this);

    this->treeElement = treeElement;
    this->settings = new QSettings(this);

    // restore size of window
    // ----------------------
    QRect parentRect(parentWidget()->mapToGlobal(QPoint(0, 0)), parentWidget()->size());
    resize(settings->value("programDialog/size", size()).toSize());
    move(settings->value("programDialog/position", QStyle::alignedRect(Qt::LeftToRight, Qt::AlignCenter, size(), parentRect).topLeft()).toPoint());

    programListWidget = new ProgramListWidget(treeElement->tree, treeElement->selectedIndex(), this);
    programListWidget->setColumnWidth(settings->value("programDialog/columnWidth", 350).toInt());
    ui->verticalLayout->addWidget(programListWidget);

    if(programListWidget->programCount() < treeElement->valueList().size())
      programListWidget->setUnknownProgram(treeElement->selectedValue());

    QDialogButtonBox *buttonBox = new QDialogButtonBox(QDialogButtonBox::Ok |QDialogButtonBox::Cancel);
    ui->verticalLayout->addWidget(buttonBox);

    connect(programListWidget, SIGNAL(selectionChanged(bool)), buttonBox->button(QDialogButtonBox::Ok), SLOT(setEnabled(bool)));
    connect(buttonBox, &QDialogButtonBox::accepted, this, &QDialog::accept);
    connect(buttonBox, &QDialogButtonBox::rejected, this, &QDialog::reject);
    connect(programListWidget, SIGNAL(programSelected(int)), this, SLOT(programSelected(int)));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

ProgramDialog::~ProgramDialog()
{
  settings->setValue("programDialog/position",    pos());
  settings->setValue("programDialog/size",        size());
  settings->setValue("programDialog/columnWidth", programListWidget->columnWidth());

  delete ui;
  delete settings;
}

/***********************************************/

void ProgramDialog::programSelected(int index)
{
  treeElement->changeSelectedIndex(index);
  QDialog::accept();
}

/***********************************************/

void ProgramDialog::accept()
{
  int index = programListWidget->selectedIndex();
  if(index > 0)
    programSelected(index);
}


/***********************************************/

void ProgramDialog::changeEvent(QEvent *e)
{
  QDialog::changeEvent(e);
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
