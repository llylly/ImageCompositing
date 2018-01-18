#ifndef IMAGELABEL_H
#define IMAGELABEL_H

#include <iostream>
#include <QObject>
#include <QWidget>
#include <QtWidgets>

using namespace std;

class ImageLabel : public QLabel {

    Q_OBJECT

public:
    ImageLabel(QWidget *parent = 0);

signals:
    void dragged(QPoint delta);
    void staticClick(QPoint pos);

protected:
    void mousePressEvent(QMouseEvent *event);
    void mouseReleaseEvent(QMouseEvent *event);
    void mouseMoveEvent(QMouseEvent *event);

private:
    QPoint startPos;
    QPoint pressPos;
};

#endif // IMAGELABEL_H
