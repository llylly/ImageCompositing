#-------------------------------------------------
#
# Project created by QtCreator 2018-01-13T13:40:29
#
#-------------------------------------------------

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = ImageCompositing
TEMPLATE = app

INCLUDEPATH += "/usr/local/include"

LIBS        += -L "/usr/local/lib" -lgmp
LIBS        += -L "/usr/local/lib" -lCGAL

# The following define makes your compiler emit warnings if you use
# any feature of Qt which as been marked as deprecated (the exact warnings
# depend on your compiler). Please consult the documentation of the
# deprecated API in order to know how to port your code away from it.
DEFINES += QT_DEPRECATED_WARNINGS

# You can also make your code fail to compile if you use deprecated APIs.
# In order to do so, uncomment the following line.
# You can also select to disable deprecated APIs only up to a certain version of Qt.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0


SOURCES += main.cpp \
    MainWindow.cpp \
    Context.cpp \
    ImageView.cpp \
    Document.cpp \
    DocumentLayer.cpp \
    NewDialog.cpp \
    LayerListWidget.cpp \
    AdjustSizeDialog.cpp \
    ImageLabel.cpp \
    QuadTreeCompositing.cpp \
    QuadTreeNode.cpp \
    StatusDialog.cpp \
    AboutDialog.cpp \
    MVCCompositingThread.cpp \
    MVCHierarchyList.cpp

HEADERS  += \
    MainWindow.h \
    Context.h \
    ImageView.h \
    Document.h \
    DocumentLayer.h \
    NewDialog.h \
    LayerListWidget.h \
    AdjustSizeDialog.h \
    ImageLabel.h \
    QuadTreeCompositing.h \
    QuadTreeNode.h \
    StatusDialog.h \
    AboutDialog.h \
    MVCCompositingThread.h \
    MVCHierarchyList.h
