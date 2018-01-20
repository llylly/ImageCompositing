#include "LayerListWidget.h"

LayerListWidget::LayerListWidget(QWidget *parent): QGroupBox(parent) {
    this->document = NULL;

    this->setTitle(QString("Layers"));
    this->mainLayout = new QVBoxLayout(this);
    this->mainLayout->setSpacing(1);
    this->buttonLayout = new QHBoxLayout();
    this->buttonLayout->setSpacing(1);
    this->buttonLayout->setMargin(0);
    this->buttonLayout->addSpacerItem(new QSpacerItem(10, 10, QSizePolicy::Expanding));
    this->upButton = new QPushButton("↑");
    this->upButton->setSizePolicy(QSizePolicy::Minimum, QSizePolicy::Minimum);
    this->downButton = new QPushButton("↓");
    this->downButton->setSizePolicy(QSizePolicy::Minimum, QSizePolicy::Minimum);
    this->delButton = new QPushButton("X");
    this->delButton->setSizePolicy(QSizePolicy::Minimum, QSizePolicy::Minimum);
    this->buttonLayout->addWidget(this->upButton);
    this->buttonLayout->addWidget(this->downButton);
    this->buttonLayout->addWidget(this->delButton);
    this->mainLayout->addLayout(this->buttonLayout);
    this->listWidget = new QListWidget();
    this->listWidget->setSelectionMode(QAbstractItemView::SingleSelection);
    this->mainLayout->addWidget(this->listWidget);

    this->buttonDisable();

    connect(this->listWidget, SIGNAL(itemClicked(QListWidgetItem*)), this, SLOT(clickButtonUpdate(QListWidgetItem*)));
    connect(this->listWidget, SIGNAL(doubleClicked(QModelIndex)), this, SLOT(doubleClicked(QModelIndex)));
    connect(this->upButton, SIGNAL(clicked(bool)), this, SLOT(upClicked()));
    connect(this->downButton, SIGNAL(clicked(bool)), this, SLOT(downClicked()));
    connect(this->delButton, SIGNAL(clicked(bool)), this, SLOT(delClicked()));
}

LayerListWidget::~LayerListWidget() {

}

void LayerListWidget::updateView() {
    this->listWidget->clear();
    if (this->document != NULL) {
        int totLayers = this->document->layers.size();
        for (int i=0; i<totLayers; ++i) {
            this->listWidget->addItem(QString(this->document->layers[i]->name.c_str()));
        }
    }
    this->buttonDisable();
}

void LayerListWidget::setNewDocument(Document *document) {
    this->document = document;
}

void LayerListWidget::buttonDisable() {
    this->upButton->setEnabled(false);
    this->downButton->setEnabled(false);
    this->delButton->setEnabled(false);
}

void LayerListWidget::clickButtonUpdate(QListWidgetItem*) {
    int ind = this->listWidget->currentIndex().row();
    if (ind > 0) this->upButton->setEnabled(true); else this->upButton->setEnabled(false);
    if (ind < (this->listWidget->count() - 1)) this->downButton->setEnabled(true); else this->downButton->setEnabled(false);
    this->delButton->setEnabled(true);
}

void LayerListWidget::doubleClicked(QModelIndex index) {
    int ind = index.row();
    if ((ind >= 0) && (ind < int(this->document->layers.size()))) {
        this->document->selectLayer(ind);
    }
}

void LayerListWidget::upClicked() {
    int ind = this->listWidget->currentIndex().row();
    if ((ind >= 0) && (ind < int(this->document->layers.size())))
        this->document->upLayer(ind);
}

void LayerListWidget::downClicked() {
    int ind = this->listWidget->currentIndex().row();
    if ((ind >= 0) && (ind < int(this->document->layers.size())))
        this->document->downLayer(ind);
}

void LayerListWidget::delClicked() {
    int ind = this->listWidget->currentIndex().row();
    if ((ind >= 0) && (ind < int(this->document->layers.size())))
        this->document->deleteLayer(ind);
}


