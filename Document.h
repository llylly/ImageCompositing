#ifndef DOCUMENT_H
#define DOCUMENT_H

#include <vector>
#include <iostream>
#include <QPixmap>
#include "DocumentLayer.h"

using namespace std;

class MainWindow;

class Document {

public:
    Document(MainWindow *theWindow, int width, int height);
    void addLayer(DocumentLayer *layer);
    void addLayerAtBottom(DocumentLayer *layer);
    QPixmap *genViewImage(bool forSave = false);

    MainWindow *theWindow;

    vector<DocumentLayer*> layers;
    int width, height;

    vector<QPoint> choosedPoints;

    int selectedLayer;

    void autoAdjustSize();
    void upLayer(int index);
    void downLayer(int index);
    void deleteLayer(int index);
    void selectLayer(int index);
    void crop();
    void cropByImageMask(QImage *maskimg);
    void unSelect();

    void save(string s);

    void addPoint(QPoint newPoint);
    void deletePoints();
    bool checkInside(int x, int y);
};

#endif // DOCUMENT_H
