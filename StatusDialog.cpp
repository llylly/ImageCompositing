#include "StatusDialog.h"

StatusDialog::StatusDialog(QWidget *widget): QDialog(widget) {
    this->setWindowTitle("Quad Tree Compositing is Running...");
    this->mainLayout = new QVBoxLayout();
    this->listWidget = new QListWidget();
    this->mainLayout->addWidget(new QLabel("Status"));
    this->mainLayout->addWidget(this->listWidget);
    this->setLayout(this->mainLayout);
}

void StatusDialog::cleanStatus() {
    this->listWidget->clear();
}

void StatusDialog::addStatus(QString s) {
    this->listWidget->addItem(s);
    this->listWidget->verticalScrollBar()->setValue(this->listWidget->verticalScrollBar()->maximum());
}
