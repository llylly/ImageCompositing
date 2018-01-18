#ifndef LAYERLISTWIDGET_H
#define LAYERLISTWIDGET_H

#include <QObject>
#include <QtWidgets>
#include "Document.h"

class LayerListWidget: public QGroupBox {

    Q_OBJECT

public:
    LayerListWidget(QWidget *parent = Q_NULLPTR);
    virtual ~LayerListWidget();
    void updateView();
    void setNewDocument(Document *document);

    Document *document;

private slots:
    void buttonDisable();
    void clickButtonUpdate(QListWidgetItem *curItem);
    void doubleClicked(QModelIndex index);
    void upClicked();
    void downClicked();
    void delClicked();

private:
    QVBoxLayout *mainLayout;
    QHBoxLayout *buttonLayout;
    QPushButton *upButton, *downButton, *delButton;
    QListWidget *listWidget;

};

#endif // LAYERLISTWIDGET_H
