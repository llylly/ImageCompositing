#include "ImageLabel.h"

ImageLabel::ImageLabel(QWidget *parent): QLabel(parent) {
    this->setAcceptDrops(true);
}

void ImageLabel::mousePressEvent(QMouseEvent *event) {
    if (event->button() == Qt::LeftButton) {
        startPos = event->pos();
        pressPos = event->pos();
    }
    QLabel::mousePressEvent(event);
}

void ImageLabel::mouseReleaseEvent(QMouseEvent *event) {
    int distance = (event->pos() - pressPos).manhattanLength();
    if (distance <= QApplication::startDragDistance()) {
        emit staticClick(event->pos());
    }
    QLabel::mouseReleaseEvent(event);
}

void ImageLabel::mouseMoveEvent(QMouseEvent *event) {
    if (event->buttons() & Qt::LeftButton) {
        int distance = (event->pos() - startPos).manhattanLength();
        if (distance >= QApplication::startDragDistance()) {
            emit dragged(event->pos() - startPos);
            startPos = event->pos();
        }
    }
    QLabel::mouseMoveEvent(event);
}


