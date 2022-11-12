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

#include "ui_settingsCommandDialog.h"
#include <QtDebug>
#include <QSettings>
#include <QList>
#include <QAbstractButton>
#include <QDropEvent>
#include <QScreen>
#include "base/importGroops.h"
#include "settingsCommandDialog.h"

/***********************************************/

SettingsCommandDialog::SettingsCommandDialog(QWidget *parent) : QDialog(parent), ui(new Ui::SettingsCommandDialog)
{
  try
  {
    ui->setupUi(this);
    this->settings = new QSettings(this);

    // restore size of window
    // ----------------------
    setMinimumSize(QGuiApplication::primaryScreen()->size()/4);
    QRect parentRect(parentWidget()->mapToGlobal(QPoint(0, 0)), parentWidget()->size());
    resize(settings->value("settingsCommandDialog/size", minimumSizeHint()).toSize());
    move(settings->value("settingsCommandDialog/position", QStyle::alignedRect(Qt::LeftToRight, Qt::AlignCenter, size(), parentRect).topLeft()).toPoint());

    // commands
    // --------
    QStringList labelList, commandList;
    readCommandList(settings, labelList, commandList);
    int commandIndex = settings->value("execute/commandIndex", int(0)).toInt();

    // init table
    // ----------
    QStringList headerLable;
    ui->table->setColumnCount(2);
    ui->table->setRowCount(labelList.size()+1);
    ui->table->horizontalHeader()->setStretchLastSection(true);
    ui->table->setColumnWidth(0, settings->value("settingsCommandDialog/columnWidth", int(200)).toInt());
    ui->table->setHorizontalHeaderLabels(headerLable << "Label" << "Command");
    ui->table->setSelectionMode(QAbstractItemView::SingleSelection);
    ui->table->setSelectionBehavior(QAbstractItemView::SelectRows);

    for(int i=0; i<labelList.size();i++)
      insertItem(i, labelList[i], commandList[i]);
    insertItem(labelList.size(), "", "");
    ui->table->setCurrentCell(commandIndex, 1);

    connect(ui->table,        SIGNAL(cellChanged(int, int)),     this, SLOT(tableCellChanged(int, int)));
    connect(ui->stdButtonBox, SIGNAL(clicked(QAbstractButton*)), this, SLOT(clickedStdBox(QAbstractButton*)));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

SettingsCommandDialog::~SettingsCommandDialog()
{
  delete ui;
  delete settings;
}

/***********************************************/

void SettingsCommandDialog::readCommandList(QSettings *settings, QStringList &labelList, QStringList &commandList)
{
  try
  {
    labelList    = settings->value("execute/commandLabels").toStringList();
    commandList  = settings->value("execute/commands").toStringList();
    if(labelList.size() != commandList.size())
      throw(Exception("size mismatch between command labels and commands"));
    if(labelList.isEmpty())
    {
      labelList<<"groops (Windows)"<<"groops (KDE)"<<"groops (GNOME)";
      labelList<<"groopsMPI (Windows, 4 processes)"<<"groopsMPI (KDE, 4 processes)"<<"groopsMPI (GNOME, 4 processes)";
      commandList<<"cd /d %w && groops.exe %f";
      commandList<<"konsole --workdir %w -e bash -ic \"groops %f; bash\"";
      commandList<<"gnome-terminal --working-directory=%w -x bash -ic \"groops %f; bash\"";
      commandList<<"cd /d %w && mpiexec /genv OMP_NUM_THREADS 1 -n 4 groopsMPI.exe %f";
      commandList<<"konsole --workdir %w -e bash -ic \"OMP_NUM_THREADS=1 mpiexec -x OMP_NUM_THREADS -n 4 groopsMPI %f; bash\"";
      commandList<<"gnome-terminal --working-directory=%w -x bash -ic \"OMP_NUM_THREADS=1 mpiexec -x OMP_NUM_THREADS -n 4 groopsMPI %f; bash\"";
      settings->setValue("execute/commandLabels", labelList);
      settings->setValue("execute/commands",      commandList);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

void SettingsCommandDialog::insertItem(int row, const QString &label, const QString &command)
{
  try
  {
    QTableWidgetItem *item1 = new QTableWidgetItem(label);
    QTableWidgetItem *item2 = new QTableWidgetItem(command);
    item1->setFlags(Qt::ItemIsEnabled | Qt::ItemIsSelectable | Qt::ItemIsEditable);
    item2->setFlags(Qt::ItemIsEnabled | Qt::ItemIsSelectable | Qt::ItemIsEditable);
    ui->table->setItem(row, 0, item1);
    ui->table->setItem(row, 1, item2);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

void SettingsCommandDialog::clickedUp()
{
  try
  {
    int index = ui->table->currentRow();
    if((index<1) || (index >= ui->table->rowCount()-1))
      return;

    QString label, command;
    if(ui->table->item(index, 0)) label   = ui->table->item(index, 0)->text();
    if(ui->table->item(index, 1)) command = ui->table->item(index, 1)->text();

    ui->table->removeRow(index);
    ui->table->insertRow(index-1);
    insertItem(index-1, label, command);
    ui->table->setCurrentCell(index-1, 1);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

void SettingsCommandDialog::clickedDown()
{
  try
  {
    int index = ui->table->currentRow();
    if((index<0) || (index >= ui->table->rowCount()-2))
      return;

    QString label, command;
    if(ui->table->item(index, 0)) label   = ui->table->item(index, 0)->text();
    if(ui->table->item(index, 1)) command = ui->table->item(index, 1)->text();

    ui->table->removeRow(index);
    ui->table->insertRow(index+1);
    insertItem(index+1, label, command);
    ui->table->setCurrentCell(index+1, 1);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

void SettingsCommandDialog::clickedApply()
{
  try
  {
    QStringList labelList, commandList;
    for(int i=0; i<ui->table->rowCount()-1; i++)
    {
      QString label, command;
      if(ui->table->item(i, 0)) label   = ui->table->item(i, 0)->text();
      if(ui->table->item(i, 1)) command = ui->table->item(i, 1)->text();

      if((!label.isEmpty()) || (!command.isEmpty()))
      {
        labelList.append(label);
        commandList.append(command);
      }
      else
        ui->table->removeRow(i--);
    }
    int commandIndex = ui->table->currentRow();

    settings->setValue("execute/commandLabels",             labelList);
    settings->setValue("execute/commands",                  commandList);
    settings->setValue("execute/commandIndex",              commandIndex);
    settings->setValue("settingsCommandDialog/position",    pos());
    settings->setValue("settingsCommandDialog/size",        size());
    settings->setValue("settingsCommandDialog/columnWidth", ui->table->columnWidth(0));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

void SettingsCommandDialog::clickedStdBox(QAbstractButton *button)
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

void SettingsCommandDialog::tableCellChanged(int /*row*/, int /*column*/)
{
  try
  {
    QString label, command;
    int rows  = ui->table->rowCount();
    if(ui->table->item(rows-1, 0)) label   = ui->table->item(rows-1, 0)->text();
    if(ui->table->item(rows-1, 1)) command = ui->table->item(rows-1, 1)->text();

    if((!label.isEmpty()) || (!command.isEmpty()))
    {
      ui->table->insertRow(rows);
      insertItem(rows, "", "");
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

void SettingsCommandDialog::changeEvent(QEvent *e)
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
