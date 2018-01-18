#ifndef NEWDIALOG_H
#define NEWDIALOG_H

#include <QObject>
#include <QDialog>
#include <QtWidgets>

class NewDialog : public QDialog {
    Q_OBJECT

public:
    NewDialog(QWidget *parent);

    int width, height;
    bool available;

    QVBoxLayout *mainLayout;
    QHBoxLayout *widthLayout, *heightLayout, *buttonLayout;
    QLabel *widthLabel, *heightLabel;
    QLineEdit *widthEdit, *heightEdit;
    QPushButton *okButton, *cancelButton;

private slots:
    void ok();
};

#endif // NEWDIALOG_H
