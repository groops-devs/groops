#-------------------------------------------------
#
# Project created by QtCreator 2011-02-19T14:42:42
#
#-------------------------------------------------

!defined(GROOPS_DIR, var) {
  GROOPS_DIR = ../
}

!defined(OBJ_DIR, var) {
  OBJ_DIR=build
}

!defined(DESTDIR, var) {
  DESTDIR=$$GROOPS_DIR/bin
}

!defined(TARGET, var) {
  TARGET=groopsGui
}

QT += xml
QT += widgets

CONFIG += c++14

win32:RC_ICONS += resources/groops.ico

TEMPLATE = app

INCLUDEPATH += $$GROOPS_DIR/source

OBJECTS_DIR=$$OBJ_DIR
MOC_DIR=$$OBJ_DIR
UI_DIR=$$OBJ_DIR
RCC_DIR=$$OBJ_DIR

SOURCES += main.cpp \
           base/xml.cpp \
           base/schema.cpp \
           tree/tree.cpp \
           tree/treeItem.cpp \
           tree/treeElement.cpp \
           tree/treeElementSimple.cpp \
           tree/treeElementBool.cpp \
           tree/treeElementFileName.cpp \
           tree/treeElementTime.cpp \
           tree/treeElementComplex.cpp \
           tree/treeElementSequence.cpp \
           tree/treeElementChoice.cpp \
           tree/treeElementAdd.cpp \
           tree/treeElementGlobal.cpp \
           tree/treeElementProgram.cpp \
           tree/treeElementUnknown.cpp \
           executeDialog/executeDialog.cpp \
           findReplaceDock/findReplaceDock.cpp \
           programDialog/programDialog.cpp \
           settingsDialog/settingsCommandDialog.cpp \
           settingsDialog/settingsPathDialog.cpp \
           mainWindow/tabs.cpp \
           mainWindow/mainWindow.cpp \
           mainWindow/sideBar.cpp \
           addGlobalDialog/addGlobalDialog.cpp \
           setLoopConditionDialog/setLoopConditionDialog.cpp \
           mainWindow/schemaSelector.cpp \
           $$GROOPS_DIR/source/parser/expressionParser.cpp \
           $$GROOPS_DIR/source/parser/stringParser.cpp \
           $$GROOPS_DIR/source/base/format.cpp \
           $$GROOPS_DIR/source/base/time.cpp \
           $$GROOPS_DIR/source/base/constants.cpp \

HEADERS  += \
            base/xml.h \
            base/schema.h \
            tree/tree.h \
            tree/treeItem.h \
            tree/treeElement.h \
            tree/treeElementSimple.h \
            tree/treeElementBool.h \
            tree/treeElementFileName.h \
            tree/treeElementTime.h \
            tree/treeElementComplex.h \
            tree/treeElementSequence.h \
            tree/treeElementChoice.h \
            tree/treeElementAdd.h \
            tree/treeElementGlobal.h \
            tree/treeElementProgram.h \
            tree/treeElementUnknown.h \
            executeDialog/executeDialog.h \
            findReplaceDock/findReplaceDock.h \
            programDialog/programDialog.h \
            settingsDialog/settingsCommandDialog.h \
            settingsDialog/settingsPathDialog.h \
            mainWindow/tabs.h \
            mainWindow/mainWindow.h \
            mainWindow/sideBar.h \
            addGlobalDialog/addGlobalDialog.h \
            setLoopConditionDialog/setLoopConditionDialog.h \
            mainWindow/schemaSelector.h \
            base/importGroops.h \
            $$GROOPS_DIR/source/parser/expressionParser.h \
            $$GROOPS_DIR/source/parser/stringParser.h \
            $$GROOPS_DIR/source/base/format.h \
            $$GROOPS_DIR/source/base/exception.h \
            $$GROOPS_DIR/source/base/importStd.h \
            $$GROOPS_DIR/source/base/time.h \
            $$GROOPS_DIR/source/base/portable.h \
            $$GROOPS_DIR/source/base/constants.h \

FORMS += mainWindow/mainWindow.ui \
         executeDialog/executeDialog.ui \
         settingsDialog/settingsCommandDialog.ui \
         findReplaceDock/findReplaceDock.ui \
         programDialog/programDialog.ui \
         settingsDialog/settingsPathDialog.ui \
         addGlobalDialog/addGlobalDialog.ui \
         setLoopConditionDialog/setLoopConditionDialog.ui \
         mainWindow/schemaSelector.ui \

RESOURCES += resources/icons.qrc
