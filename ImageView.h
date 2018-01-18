#ifndef IMAGEVIEW_H
#define IMAGEVIEW_H

#include <QObject>
#include <QWidget>
#include <QtWidgets>
#include "ImageLabel.h"

using namespace std;

class ImageView: public QWidget {
    Q_OBJECT

public:
    ImageView(QWidget *parent = Q_NULLPTR);
    virtual ~ImageView();

    void setImage(QPixmap *newPic);
    void updateView();

    void showButton();
    void hideButton();

    ImageLabel *label;
    QPushButton *maskCropButton, *saveButton, *cancelButton;

private:
    QPixmap *nowPic;
    QScrollArea *scrollArea;
    QVBoxLayout *mainLayout;
    QWidget *buttonWidget;
};

#endif // IMAGEVIEW_H
