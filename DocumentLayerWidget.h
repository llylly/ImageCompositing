#ifndef DOCUMENTLAYERWIDGET_H
#define DOCUMENTLAYERWIDGET_H

#include <QObject>
#include <QtWidgets>

class LayerListWidget;

class DocumentLayerWidget: public QListWidgetItem {

public:
    DocumentLayerWidget(int index, LayerListWidget *parent = Q_NULLPTR, int type = Type);
    virtual ~DocumentLayerWidget();

    QWidget *mainWidget;

private:
    int index;

    QHBoxLayout *mainLayout;
    QLabel *indexLabel;
    QSpacerItem *spacer;
    QPushButton *upButton, *downButton, *delButton;

};

#endif // DOCUMENTLAYERWIDGET_H
