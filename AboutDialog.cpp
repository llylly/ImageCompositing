#include "AboutDialog.h"

AboutDialog::AboutDialog(QWidget *widget): QDialog(widget) {

    this->setWindowTitle("About");
    QVBoxLayout *layout = new QVBoxLayout();
    layout->addWidget(new QLabel("Media Computing Course Project"));

    layout->addWidget(new QLabel("By Linyi Li @ Tsinghua"));

    QPushButton *closeButton = new QPushButton("Ok");
    layout->addWidget(closeButton);
    this->setLayout(layout);

    connect(closeButton, SIGNAL(clicked(bool)), this, SLOT(close()));
}
