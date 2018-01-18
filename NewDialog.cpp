#include "NewDialog.h"

NewDialog::NewDialog(QWidget *parent): QDialog(parent) {
    this->width = 800;
    this->height = 600;
    this->available = false;

    this->mainLayout = new QVBoxLayout();
    this->widthLayout = new QHBoxLayout();
    this->widthLabel = new QLabel("Width");
    this->widthEdit = new QLineEdit(QString::number(this->width));
    this->widthLayout->addWidget(this->widthLabel);
    this->widthLayout->addWidget(this->widthEdit);
    this->mainLayout->addLayout(this->widthLayout);
    this->heightLayout = new QHBoxLayout();
    this->heightLabel = new QLabel("Height");
    this->heightEdit = new QLineEdit(QString::number(this->height));
    this->heightLayout->addWidget(this->heightLabel);
    this->heightLayout->addWidget(this->heightEdit);
    this->mainLayout->addLayout(this->heightLayout);
    this->buttonLayout = new QHBoxLayout();
    this->okButton = new QPushButton("Ok");
    this->okButton->setFocus();
    this->cancelButton = new QPushButton("Cancel");
    this->buttonLayout->addWidget(this->okButton);
    this->buttonLayout->addWidget(this->cancelButton);
    this->mainLayout->addLayout(this->buttonLayout);
    this->setLayout(this->mainLayout);

    connect(this->okButton, SIGNAL(clicked(bool)), this, SLOT(ok()));
    connect(this->cancelButton, SIGNAL(clicked(bool)), this, SLOT(close()));
}

void NewDialog::ok() {
    bool *legal = new bool(false);
    int h = this->heightEdit->text().toInt(legal);
    int w = this->widthEdit->text().toInt(legal);
    if ((h < 0) || (w < 0) || (h > 10000) || (w > 10000)) *legal = false;
    if (!(*legal)) {
        QMessageBox::information(this, "New Document", "The height and width should be integers between 1 and 10000.", QMessageBox::Ok);
    } else {
        this->height = h;
        this->width = w;
        this->available = true;
        this->close();
    }
    delete legal;
}
