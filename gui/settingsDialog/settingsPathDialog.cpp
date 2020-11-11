/***********************************************/
/**
* @file settingsPathDialog.cpp
*
* @brief Dialog to set up default paths and files.
*
* @author Torsten Mayer-Guerr
* @author Sebastian Strasser
* @date 2006-10-13
*/
/***********************************************/

#include "ui_settingsPathDialog.h"
#include <QtDebug>
#include <QSettings>
#include <QFileDialog>
#include <QMessageBox>
#include <QSysInfo>
#include <QScreen>
#include "base/importGroops.h"
#include "base/schema.h"
#include "settingsPathDialog.h"

/***********************************************/

SettingsPathDialog::SettingsPathDialog(QWidget *parent) : QDialog(parent), ui(new Ui::SettingsPathDialog)
{
  try
  {
    ui->setupUi(this);
    this->settings = new QSettings(this);

    // restore size of window
    // ----------------------
    setMinimumSize(QGuiApplication::primaryScreen()->size()/4);
    QRect parentRect(parentWidget()->mapToGlobal(QPoint(0, 0)), parentWidget()->size());
    resize(settings->value("settingsPathDialog/size", minimumSizeHint()).toSize());
    move(settings->value("settingsPathDialog/position", QStyle::alignedRect(Qt::LeftToRight, Qt::AlignCenter, size(), parentRect).topLeft()).toPoint());

    connect(ui->listSchema->model(), SIGNAL(rowsInserted(QModelIndex,int,int)), this, SLOT(schemaListChanged()));
    connect(ui->listSchema->model(), SIGNAL(rowsRemoved (QModelIndex,int,int)), this, SLOT(schemaListChanged()));

    const QString baseDir(QSysInfo::productType() == "windows" ? "C:" : QDir::homePath());
    const QString schemaFile = settings->value("files/schemaFile", QFileInfo(baseDir+"/groops/groops.xsd").absoluteFilePath()).toString();
    ui->listSchema->addItems(settings->value("files/schemaFiles", QStringList(schemaFile)).toStringList());
    ui->editTemplateFile->setText(settings->value("files/templateFile").toString());
    ui->editWorkingDir->setText(settings->value("files/workingDirectory", QDir(baseDir+"/groops/scenario").absolutePath()).toString());
    ui->editDocumentationDir->setText(settings->value("files/documentationDirectory", QDir(baseDir+"/groops/docs/html").absolutePath()).toString());
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

SettingsPathDialog::~SettingsPathDialog()
{
  delete ui;
  delete settings;
}

void SettingsPathDialog::schemaListChanged()
{
   ui->buttonSchemaChange->setEnabled(ui->listSchema->count() > 0);
   ui->buttonSchemaRemove->setEnabled(ui->listSchema->count() > 1);
}

/***********************************************/

void SettingsPathDialog::clickedSchemaChange()
{
  QListWidgetItem *item = ui->listSchema->currentItem();
  if(!item)
    return;

  QString name = QFileDialog::getOpenFileName(this, tr("Change Schema File - GROOPS"),
                                              item->text(), tr("XML Schema files (*.xsd)"));
  if(!name.isEmpty() && !ui->listSchema->findItems(name, Qt::MatchFixedString).size())
  {
    // validate XSD schema file
    if(!Schema::validateSchema(name))
    {
      QMessageBox::StandardButton button =
        QMessageBox::critical(this , tr("Change Schema File - GROOPS"),
                              tr("File '%1' seems not to be a valid XSD schema. Try another file?").arg(name),
                              QMessageBox::Ok | QMessageBox::Cancel);
      if(button == QMessageBox::Ok)
        clickedSchemaChange();
      return;
    }

    item->setText(name);
    ui->listSchema->sortItems();
    ui->listSchema->setItemSelected(item, true);
  }
}


/***********************************************/

void SettingsPathDialog::clickedSchemaAdd()
{
  QString name = QFileDialog::getOpenFileName(this, tr("Add Schema File - GROOPS"),
                                              ui->listSchema->currentItem() ? ui->listSchema->currentItem()->text() : "groops.xsd", tr("XML Schema files (*.xsd)"));
  if(!name.isEmpty() && !ui->listSchema->findItems(name, Qt::MatchFixedString).size())
  {
    // validate XSD schema file
    if(!Schema::validateSchema(name))
    {
      QMessageBox::StandardButton button =
        QMessageBox::critical(this , tr("Add Schema File - GROOPS"),
                              tr("File '%1' seems not to be a valid XSD schema. Try another file?").arg(name),
                              QMessageBox::Ok | QMessageBox::Cancel);
      if(button == QMessageBox::Ok)
        clickedSchemaAdd();
      return;
    }

    ui->listSchema->addItem(name);
    ui->listSchema->sortItems();
    ui->listSchema->setItemSelected(ui->listSchema->findItems(name, Qt::MatchFixedString).at(0), true);
  }
}

/***********************************************/

void SettingsPathDialog::clickedSchemaRemove()
{
  qDeleteAll(ui->listSchema->selectedItems());
}

/***********************************************/

void SettingsPathDialog::clickedTemplateFile()
{
  QString name = QFileDialog::getOpenFileName(this, tr("Choose System Directory - GROOPS"), ui->editTemplateFile->text(), tr("XML files (*.xml)"));
  if(!name.isEmpty())
    ui->editTemplateFile->setText(name);
}

/***********************************************/

void SettingsPathDialog::clickedWorkingDir()
{
  QString name = QFileDialog::getExistingDirectory(this, tr("Choose Working Directory - GROOPS"), ui->editWorkingDir->text());
  if(!name.isEmpty())
    ui->editWorkingDir->setText(name);
}

/***********************************************/

void SettingsPathDialog::clickedDocumentationDir()
{
  QString name = QFileDialog::getExistingDirectory(this, tr("Choose Documentation Directory - GROOPS"), ui->editDocumentationDir->text());
  if(!name.isEmpty())
    ui->editDocumentationDir->setText(name);
}


/***********************************************/

void SettingsPathDialog::clickedOk()
{
  try
  {
    QStringList schemaFiles;
    for(int i = 0; i < ui->listSchema->count(); i++)
      schemaFiles.append(ui->listSchema->item(i)->text());
    if(!schemaFiles.size())
    {
      QMessageBox::warning(this , tr("Path Settings - GROOPS"), tr("No XML Schema Files specified. At least one schema is required."), QMessageBox::Ok);
      return;
    }

    settings->setValue("files/schemaFile",       schemaFiles.at(0));
    settings->setValue("files/schemaFiles",      schemaFiles);
    settings->setValue("files/templateFile",     ui->editTemplateFile->text());
    settings->setValue("files/workingDirectory", ui->editWorkingDir->text());
    settings->setValue("files/documentationDirectory", ui->editDocumentationDir->text());
    settings->setValue("settingsPathDialog/position",  pos());
    settings->setValue("settingsPathDialog/size",      size());
    accept();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

void SettingsPathDialog::changeEvent(QEvent *e)
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
