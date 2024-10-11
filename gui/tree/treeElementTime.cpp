/***********************************************/
/**
* @file treeElementTime.cpp
*
* @brief Element with date/time editor.
*
* @author Torsten Mayer-Guerr
* @author Sebastian Strasser
* @date 2006-06-10
*/
/***********************************************/

#include <QtDebug>
#include <QHBoxLayout>
#include <QComboBox>
#include <QDateTimeEdit>
#include "base/importGroops.h"
#include "tree/tree.h"
#include "tree/treeElement.h"
#include "tree/treeElementSimple.h"
#include "tree/treeElementGlobal.h"
#include "tree/treeElementTime.h"

/***********************************************/

void TreeElementTime::setSelectedIndex(int index)
{
  try
  {
    TreeElementSimple::setSelectedIndex(index); // updates also mjd

    if(dateTimeEdit && !changeNotDateTime) // avoid infinite recursion
    {
      changeNotComboBox = true;
      dateTimeEdit->setDateTime(mjd2date(mjd));
      changeNotComboBox = false;
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

QWidget *TreeElementTime::createEditor()
{
  try
  {
    changeNotComboBox = changeNotDateTime = false;

    // create layout
    QWidget *layoutWidget = new QWidget(tree);
    QHBoxLayout *layout   = new QHBoxLayout(layoutWidget);
    layout->setContentsMargins(0, 0, 0, 0);

    // create ComboBox
    comboBox = createComboBox(true);
    layout->addWidget(comboBox);
    layoutWidget->setFocusProxy(comboBox);
    connect(comboBox, SIGNAL(editTextChanged(const QString &)), this, SLOT(comboBoxEditTextChanged(const QString &)));

    // create DateTime-Editor
    dateTimeEdit = new QDateTimeEdit();
    dateTimeEdit-> setDisplayFormat("yyyy-MM-dd hh:mm:ss");
    dateTimeEdit->setMinimumDateTime(QDate(1858, 11, 17).startOfDay()); // (mjd = 0)
    dateTimeEdit->setDateTime(mjd2date(mjd));
    layout->addWidget(dateTimeEdit);
    connect(dateTimeEdit, SIGNAL(dateTimeChanged(const QDateTime &)), this, SLOT(dateTimeChanged(const QDateTime &)));

    return layoutWidget;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

void TreeElementTime::interact()
{
  if(!dateTimeEdit->hasFocus())
    dateTimeEdit->setFocus();
  else
    comboBox->setFocus();
}

/***********************************************/

void TreeElementTime::comboBoxEditTextChanged(const QString &text)
{
  try
  {
    if(dateTimeEdit && !changeNotDateTime) // avoid infinite recursion
    {
      mjd = 0;
      if(tree->elementGlobal)
      {
        try
        {
          bool resolved = true;
          auto result = StringParser::parse(name().toStdString(), text.toStdString(), *varList, resolved);
          if(resolved)
            mjd = ExpressionVariable::parse(result, *varList);
        }
        catch(std::exception &/*e*/)
        {} // if not a number use mjd = 0
      }

      changeNotComboBox = true;
      dateTimeEdit->setDateTime(mjd2date(mjd));
      changeNotComboBox = false;
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

void TreeElementTime::dateTimeChanged(const QDateTime &dateTime)
{
  try
  {
    if(comboBox && !changeNotComboBox) // avoid infinite recursion
    {
      changeNotDateTime = true;
      QString text = date2mjd(dateTime);
      comboBox->setEditText(text);
      comboBoxTextEdited(text);
      changeNotDateTime = false;
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/
/***********************************************/

QString TreeElementTime::date2mjd(const QDateTime &dateTime) const
{
  try
  {
    int Y   = dateTime.date().year();
    int M   = dateTime.date().month();
    int D   = dateTime.date().day();
    int mjd = (1461 * (Y + 4800 + (M - 14)/12))/4 +(367 * (M - 2 - 12 * ((M - 14)/12)))/12 - (3 * ((Y + 4900 + (M - 14)/12)/100))/4 + D - 32075-2400001;

    int    hour = dateTime.time().hour();
    int    min  = dateTime.time().minute();
    double sec  = dateTime.time().second(); //+dateTime.time().msec()/1000.;

    QString     string;
    QTextStream stream(&string);
    bool        set=false;
    if(mjd!=0)  {stream<<mjd; set=true;}
    if(hour!=0) {stream<<((set)?"+":"")<<hour<<"/24";   set=true;}
    if(min!=0)  {stream<<((set)?"+":"")<<min<<"/1440";  set=true;}
    if(sec!=0.) {stream<<((set)?"+":"")<<sec<<"/86400"; set=true;}
    return string;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

QDateTime TreeElementTime::mjd2date(Double mjd) const
{
  try
  {
    Time time = mjd2time(mjd);
    UInt year, month, day, hour, minute;
    Double second;
    time.date(year, month, day, hour, minute, second);

    QDateTime dateTime(QDate(year, month, day), QTime(hour, minute, second));
    if(std::round(second - dateTime.time().second()) > 0)
      dateTime = dateTime.addSecs(1);
    return dateTime;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e);
  }
}

/***********************************************/

QString TreeElementTime::parseExpression(const QString &text, const VariableList &varList) const
{
  mjd = 0;
  QString result = text;
  try
  {
    bool resolved = true;
    result = QString::fromStdString(StringParser::parse(name().toStdString(), text.toStdString(), varList, resolved));
    if(result.isEmpty() || !resolved)
      return result;
    mjd = ExpressionVariable::parse(result.toStdString(), varList);
    result.setNum(mjd, 'f', 7).remove(QRegularExpression("0+$")).remove(QRegularExpression("\\.$")); // %.7f with trailing zeros removed
    if(mjd == 0)
      return result;

    // is a date?
    if(mjd > 1000)
    {
      QString dateString = mjd2date(mjd).toString("yyyy-MM-dd hh:mm:ss");
      if(dateString.endsWith(" 00:00:00"))
        dateString.truncate(dateString.lastIndexOf(" 00:00:00"));
      return dateString.trimmed();
    }

    // time span
    QDateTime dateTime = mjd2date(std::fabs(mjd));
    if(dateTime.time() == QTime(0,0,0))
      return result;
    QString dateString;
    if(mjd < 0)
      dateString += "-";
    if(std::floor(std::abs(mjd)))
      dateString += QString::number(static_cast<int>(std::floor(std::abs(mjd))))+"d "; // days
    dateString += dateTime.toString("hh:mm:ss");
    if(result != text)
      dateString += " ("+result+")";
    return dateString;
  }
  catch(std::exception &/*e*/)
  {
  }

  return result;
}

/***********************************************/
/***********************************************/
