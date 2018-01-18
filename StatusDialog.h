#ifndef STATUSDIALOG_H
#define STATUSDIALOG_H

#include <QObject>
#include <QtWidgets>

using namespace std;

class StatusDialog: public QDialog {
    Q_OBJECT

public:
    StatusDialog(QWidget *widget);

public slots:
    void cleanStatus();
    void addStatus(QString s);

private:
    QVBoxLayout *mainLayout;
    QListWidget *listWidget;


};

#endif // STATUSDIALOG_H
