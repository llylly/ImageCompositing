#ifndef ADJUSTSIZEDIALOG_H
#define ADJUSTSIZEDIALOG_H

#include <iostream>
#include <QObject>
#include <QtWidgets>

using namespace std;

class AdjustSizeDialog : public QDialog {
    Q_OBJECT

public:
    AdjustSizeDialog(QWidget *parent);

    int width, height;
    bool available;

    QVBoxLayout *mainLayout;
    QHBoxLayout *widthLayout, *heightLayout, *buttonLayout;
    QLabel *widthLabel, *heightLabel;
    QLineEdit *widthEdit, *heightEdit;
    QPushButton *okButton, *cancelButton;

    void initRefresh(int width, int height);

private slots:
    void ok();
};

#endif // ADJUSTSIZEDIALOG_H
