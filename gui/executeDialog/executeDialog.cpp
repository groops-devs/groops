/***********************************************/
/**
* @file executeDialog.cpp
*
* @brief Execute dialog.
*
* @author Torsten Mayer-Guerr
* @author Sebastian Strasser
* @date 2006-10-13
*
*/
/***********************************************/

#include "ui_executeDialog.h"
#include <QFileDialog>
#include <QSettings>
#include <QList>
#include <QCheckBox>
#include <QScreen>
#include "base/importGroops.h"
#include "tree/tree.h"
#include "tree/treeElement.h"
#include "tree/treeElementAdd.h"
#include "tree/treeElementProgram.h"
#include "settingsDialog/settingsCommandDialog.h"
#include "executeDialog.h"

/***********************************************/

ExecuteDialog::ExecuteDialog(Tree *tree, QWidget *parent)
  : QDialog(parent), ui(new Ui::ExecuteDialog)
{
  try
  {
    ui->setupUi(this);
    this->tree = tree;

    // restore size of window
    // ----------------------
    setMinimumSize(QGuiApplication::primaryScreen()->size()/4);
    QRect parentRect(parentWidget()->mapToGlobal(QPoint(0, 0)), parentWidget()->size());
    resize(settings.value("executeDialog/size", minimumSizeHint()).toSize());
    move(settings.value("executeDialog/position", QStyle::alignedRect(Qt::LeftToRight, Qt::AlignCenter, size(), parentRect).topLeft()).toPoint());

    // init table
    // ----------
    QStringList headerLable;
    ui->tableWidget->setHeaderLabels(headerLable << "Program" << "Comment");
    ui->tableWidget->setColumnWidth(0, settings.value("executeDialog/columnWidth", 350).toInt());
    ui->tableWidget->setAlternatingRowColors(true);
    ui->tableWidget->setSelectionMode(QAbstractItemView::SingleSelection);

    // set rectangular icon size
    // -------------------------
    QCheckBox *checkBox = new QCheckBox(this);
    int height = checkBox->sizeHint().height();
    ui->tableWidget->setIconSize(QSize(2*height, height));
    delete checkBox;

    // fill list
    // ---------
    addChildren(ui->tableWidget->invisibleRootItem(), tree->rootElement);

    QStringList labelList, commandList;
    SettingsCommandDialog::readCommandList(labelList, commandList);

    ui->comboBoxCommand->addItems( labelList );
    ui->comboBoxCommand->setCurrentIndex( settings.value("execute/commandIndex", int(0)).toInt() );
    ui->checkBoxLogFile->setChecked( settings.value("execute/useLogFile", bool(false)).toBool() );
    ui->editorLogfile->setText( settings.value("execute/logFile", QString("groops.log")).toString() );
    ui->editorLogfile->setEnabled( ui->checkBoxLogFile->isChecked() );
    ui->buttonBrowseLogfile->setEnabled( ui->checkBoxLogFile->isChecked() );
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

ExecuteDialog::~ExecuteDialog()
{
  delete ui;
}

/***********************************************/

void ExecuteDialog::editCommand()
{
  try
  {
    settings.setValue("execute/commandIndex", ui->comboBoxCommand->currentIndex());

    SettingsCommandDialog dialog(this);
    dialog.exec();

    ui->comboBoxCommand->clear();
    ui->comboBoxCommand->addItems( settings.value("execute/commandLabels").toStringList() );
    ui->comboBoxCommand->setCurrentIndex( settings.value("execute/commandIndex", int(0)).toInt() );
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

void ExecuteDialog::toggleLogfile()
{
  try
  {
    ui->editorLogfile->setEnabled( ui->checkBoxLogFile->isChecked() );
    ui->buttonBrowseLogfile->setEnabled( ui->checkBoxLogFile->isChecked() );
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

void ExecuteDialog::browseLogfile()
{
  try
  {
    QString name = tree->addWorkingDirectory(ui->editorLogfile->text());
    name = QFileDialog::getSaveFileName(this, tr("Choose logfile - GROOPS"),
                                        name, "", nullptr, QFileDialog::DontConfirmOverwrite);
    if(!name.isEmpty())
      ui->editorLogfile->setText(tree->stripWorkingDirectory(name));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

void ExecuteDialog::clickedAll()
{
  try
  {
    for(int i=0; i<itemList.size(); i++)
      if(!itemList[i]->parent()) // only top level items
        itemList[i]->setCheckState(0, Qt::Checked);
    ui->tableWidget->update();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/*******************************************************************************/

void ExecuteDialog::clickedNone()
{
  try
  {
    for(int i=0; i<itemList.size(); i++)
      if(!itemList[i]->parent()) // only top level items
        itemList[i]->setCheckState(0, Qt::Unchecked);
    ui->tableWidget->update();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

void ExecuteDialog::clickedApply()
{
  try
  {
    for(int i=0; i<programList.size(); i++)
      programList[i]->setDisabled( (itemList[i]->checkState(0)==Qt::Unchecked) );
    tree->fileSave();

    settings.setValue("execute/useLogFile",   ui->checkBoxLogFile->isChecked());
    settings.setValue("execute/logFile",      ui->editorLogfile->text());
    settings.setValue("execute/commandIndex", ui->comboBoxCommand->currentIndex());

    settings.setValue("executeDialog/position",    pos());
    settings.setValue("executeDialog/size",        size());
    settings.setValue("executeDialog/columnWidth", ui->tableWidget->columnWidth(0));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

void ExecuteDialog::clickedStdBox(QAbstractButton *button)
{
  try
  {
    switch(ui->stdButtonBox->standardButton(button))
    {
      case QDialogButtonBox::Ok:
        clickedApply();
        accept();
        break;
      case QDialogButtonBox::Apply:
        clickedApply();
        break;
      case QDialogButtonBox::Cancel:
        reject();
        break;
      default:
        qWarning("Unknown button clicked");
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

void ExecuteDialog::changeEvent(QEvent *e)
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

void ExecuteDialog::addChildren(QTreeWidgetItem *parentItem, TreeElement *parentElement)
{
  try
  {
    TreeElementComplex *parent = dynamic_cast<TreeElementComplex*>(parentElement);

    for(auto &element : parent->children())
      if(dynamic_cast<TreeElementProgram*>(element))
      {
        QTreeWidgetItem *item = new QTreeWidgetItem();
        item->setFlags(Qt::ItemIsEnabled | Qt::ItemIsSelectable | Qt::ItemIsUserCheckable);
        item->setIcon(0, QIcon(":/icons/scalable/program.svg"));
        item->setText(0, element->selectedValue());
        item->setText(1, element->comment());
        item->setCheckState(0, (element->disabled()) ? Qt::Unchecked : Qt::Checked);
        itemList.push_back(item);
        programList.push_back(element);
        parentItem->addChild(item);

        // recursively add children
        addChildren(item, element);
      }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/
