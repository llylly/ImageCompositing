#include "DocumentLayerWidget.h"
#include "LayerListWidget.h"

DocumentLayerWidget::DocumentLayerWidget(int index, LayerListWidget *parent, int type): QListWidgetItem(parent, type) {
    this->index = index;
    this->mainWidget = new QWidget(parent);
    this->mainLayout = new QHBoxLayout(this->mainWidget);
    this->mainWidget->setMinimumHeight(50);
    this->indexLabel = new QLabel(QString("Layer ") + QString::number(index + 1));
    this->mainLayout->addWidget(this->indexLabel);
    //this->spacer = new QSpacerItem(20, 20, QSizePolicy::Expanding);
    //this->mainLayout->addSpacerItem(this->spacer);
    this->upButton = new QPushButton("↑");
    this->mainLayout->addWidget(this->upButton);
    this->downButton = new QPushButton("↓");
    this->mainLayout->addWidget(this->downButton);
    this->delButton = new QPushButton("X");
    this->mainLayout->addWidget(this->delButton);
}

DocumentLayerWidget::~DocumentLayerWidget() {
//    this->mainLayout->removeWidget(this->indexLabel);
//    this->mainLayout->removeWidget(this->upButton);
//    this->mainLayout->removeWidget(this->downButton);
//    this->mainLayout->removeWidget(this->delButton);
//    delete this->indexLabel;
//    delete this->spacer;
//    delete this->upButton;
//    delete this->downButton;
//    delete this->delButton;
//    delete this->mainLayout;
//    delete this->mainWidget;
}

