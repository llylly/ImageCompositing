#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QApplication>
#include <QDesktopWidget>
#include <QMainWindow>
#include <QtWidgets/QHBoxLayout>
#include <QtWidgets/QVBoxLayout>
#include <QtWidgets/QLabel>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QScrollArea>
#include <QtWidgets/QListView>

#include "Context.h"
#include "ImageView.h"
#include "LayerListWidget.h"
#include "Document.h"
#include "NewDialog.h"
#include "AdjustSizeDialog.h"
#include "StatusDialog.h"
#include "AboutDialog.h"

#include "QuadTreeCompositing.h"

class MainWindow : public QMainWindow {
    Q_OBJECT

private:
    QWidget *centralWidget;
    QHBoxLayout *mainLayout;
    ImageView *imageView;
    QGridLayout *rightLayout;
    QPushButton *newButton;
    QPushButton *adjustSizeButton;
    QPushButton *autoAdjustSizeButton;
    QPushButton *addLayerButton;
    QPushButton *compositeButton;
    QPushButton *saveButton;
    QPushButton *aboutButton;
    LayerListWidget *layerView;

    NewDialog *newDialog;
    AdjustSizeDialog *adjustSizeDialog;
    StatusDialog *statusDialog;
    AboutDialog *aboutDialog;

    QuadTreeCompositing *quadTreeCompositingThread;

private slots:
    void newClick();
    void adjustSizeClick();
    void autoAdjustClick();
    void addLayerClick();
    void viewSaveClick();
    void viewCancelClick();
    void viewMaskCropClick();
    void saveClick();
    void aboutClick();
    void dragHandler(QPoint delta);
    void clickhandler(QPoint pos);
    void quadTreeCompositingToggle();

public slots:
    void viewUpdate();

    void quadTreeCompositingFinish();

public:
    MainWindow(QWidget *parent = 0);
    ~MainWindow();

    Document *curDocment;
};

#endif // MAINWINDOW_H
