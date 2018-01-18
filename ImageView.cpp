#include "ImageView.h"

ImageView::ImageView(QWidget *parent): QWidget(parent) {
    this->mainLayout = new QVBoxLayout(this);
    this->mainLayout->setSpacing(1);

    this->label = new ImageLabel(this);
    this->label->setAlignment(Qt::AlignHCenter | Qt::AlignVCenter);
    this->label->setSizePolicy(QSizePolicy::Ignored, QSizePolicy::Ignored);
    this->label->setScaledContents(true);
    this->label->setContentsMargins(0, 0, 0, 0);
    this->label->setStyleSheet("color: black");

    this->scrollArea = new QScrollArea(this);
    this->scrollArea->setBackgroundRole(QPalette::Dark);
    this->scrollArea->setWidget(this->label);
    this->scrollArea->setVisible(true);
    this->mainLayout->addWidget(this->scrollArea);

    this->buttonWidget = new QWidget(this);
    QHBoxLayout *buttonLayout = new QHBoxLayout(this);
    buttonLayout->setMargin(1);
    this->maskCropButton = new QPushButton("Crop by Mask");
    this->saveButton = new QPushButton("Crop");
    this->cancelButton = new QPushButton("√");
    buttonLayout->addWidget(this->maskCropButton);
    buttonLayout->addWidget(this->saveButton);
    buttonLayout->addSpacerItem(new QSpacerItem(10, 10, QSizePolicy::Expanding));
    buttonLayout->addWidget(this->cancelButton);
    this->buttonWidget->setLayout(buttonLayout);
    this->buttonWidget->hide();
    this->mainLayout->addWidget(this->buttonWidget);

    this->nowPic = NULL;

    this->updateView();
}

void ImageView::setImage(QPixmap *newPic) {
    if (newPic == NULL) return;
    if (this->nowPic != NULL)
        delete this->nowPic;
    this->nowPic = newPic;
}

void ImageView::updateView() {
    if (this->nowPic == NULL) {
        this->label->setText("No document for show.");
    } else {
        this->label->setText("");
        this->label->setPixmap(*(this->nowPic));
    }
    this->label->adjustSize();
}

void ImageView::showButton() {
    this->buttonWidget->show();
}

void ImageView::hideButton() {
    this->buttonWidget->hide();
}

ImageView::~ImageView() {

}
