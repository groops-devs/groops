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

#ifndef __GROOPSGUI__SCHEMASELECTOR_H__
#define __GROOPSGUI__SCHEMASELECTOR_H__

#include <QAction>
#include <QSettings>
#include <QWidget>

/***** CLASS ***********************************/

namespace Ui {
class SchemaSelector;
}

class SchemaSelector : public QWidget
{
  Q_OBJECT

public:
  explicit SchemaSelector(QWidget *parent);
 ~SchemaSelector();
  void setAction(QAction *action);
  void setCurrentTreeSchema(QString schemaFile);

public slots:
  void updateList();
  void selectionChanged();

private:
  Ui::SchemaSelector *ui;
  QSettings           settings;
  QAction            *action;
  QString             currentTreeSchema;
};

/***********************************************/

#endif // __GROOPSGUI__SCHEMASELECTOR_H__
