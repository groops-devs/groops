/***********************************************/
/**
* @file schemaSelector.h
*
* @brief Schema selector.
*
* @author Sebastian Strasser
* @date 2017-06-17
*/
/***********************************************/

#include "schemaSelector.h"
#include "ui_schemaSelector.h"

/***********************************************/

SchemaSelector::SchemaSelector(QWidget *parent) : QWidget(parent), ui(new Ui::SchemaSelector)
{
  ui->setupUi(this);
  settings = new QSettings(this);
  action   = nullptr;

  ui->icon->hide();
  int iconSize = parentWidget()->minimumSizeHint().height()*7/10;
  ui->icon->setMaximumSize(iconSize, iconSize);
  ui->icon->resize(ui->icon->maximumSize());

  updateList();

  connect(ui->comboBox, SIGNAL(currentIndexChanged(int)), this, SLOT(selectionChanged()));
}

/***********************************************/

SchemaSelector::~SchemaSelector()
{
  delete ui;
}

/***********************************************/

void SchemaSelector::setAction(QAction *action)
{
  this->action = action;

  if(action != nullptr)
    action->setVisible(ui->comboBox->count() > 1);
}

/***********************************************/

void SchemaSelector::setCurrentTreeSchema(QString schemaFile)
{
  currentTreeSchema = schemaFile;
  ui->icon->setVisible(schemaFile != ui->comboBox->currentText());
}

/***********************************************/

void SchemaSelector::updateList()
{
  QString     schemaFile  = settings->value("files/schemaFile").toString();
  QStringList schemaFiles = settings->value("files/schemaFiles").toStringList();
  ui->comboBox->blockSignals(true);
  ui->comboBox->clear();
  if(schemaFiles.size())
    ui->comboBox->addItems(schemaFiles);
  else
    ui->comboBox->addItem(schemaFile);
  int index = ui->comboBox->findText(schemaFile);
  ui->comboBox->setCurrentIndex(index >= 0 ? index : 0);
  ui->comboBox->blockSignals(false);

  if(schemaFile != ui->comboBox->currentText())
    selectionChanged();

  ui->icon->setVisible(currentTreeSchema != ui->comboBox->currentText());

  // setVisible() has to be called on QAction when QWidget is added to QMenuBar
  if(action != nullptr)
    action->setVisible(ui->comboBox->count() > 1);
}

/***********************************************/

void SchemaSelector::selectionChanged()
{
  settings->setValue("files/schemaFile", ui->comboBox->currentText());
  ui->icon->setVisible(currentTreeSchema != ui->comboBox->currentText());
}

/***********************************************/
