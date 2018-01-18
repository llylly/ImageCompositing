#ifndef ABOUTDIALOG_H
#define ABOUTDIALOG_H

#include <QObject>
#include <QDialog>
#include <QtWidgets>

class AboutDialog: public QDialog {
    Q_OBJECT

public:
    AboutDialog(QWidget *widget);
};

#endif // ABOUTDIALOG_H
