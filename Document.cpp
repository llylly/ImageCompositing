#include "Document.h"
#include "MainWindow.h"

#define EPS 1e-6

double cross(double x1, double y1, double x2, double y2) {
    return x1 * y2 - x2 * y1;
}

int mark_x[9] = {-1,-1,-1,0,0,0,1,1,1};
int mark_y[9] = {-1,0,1,-1,0,1,-1,0,1};

Document::Document(MainWindow *theWindow, int width, int height) {
    this->theWindow = theWindow;

    this->layers.clear();
    this->width = width;
    this->height = height;
    this->selectedLayer = -1;

    this->choosedPoints.clear();
}

void Document::autoAdjustSize() {
    if (this->layers.size() == 0) return;
    int w = 0, h = 0;
    for (vector<DocumentLayer*>::iterator ite = this->layers.begin(); ite != this->layers.end(); ++ite) {
        DocumentLayer *nowLayer = *ite;
        if (nowLayer->height + nowLayer->hOffset > h)
            h = nowLayer->height + nowLayer->hOffset;
        if (nowLayer->width + nowLayer->wOffset > w)
            w = nowLayer->width + nowLayer->wOffset;
    }
    this->width = w;
    this->height = h;
}

void Document::addLayer(DocumentLayer *layer) {
    this->layers.push_back(layer);
}

QPixmap *Document::genViewImage(bool forSave) {
    int *arrR = new int[this->width * this->height];
    int *arrG = new int[this->width * this->height];
    int *arrB = new int[this->width * this->height];

    for (int i=0; i<this->height; ++i)
        for (int j=0; j<this->width; ++j) {
            if (forSave)
                arrR[I(i,j,this->width)] = arrG[I(i,j,this->width)] = arrB[I(i,j,this->width)] = 255;
            else {
                if (this->selectedLayer == -1)
                    if (((i/8)&1)^((j/8)&1))
                        arrR[I(i,j,this->width)] = arrG[I(i,j,this->width)] = arrB[I(i,j,this->width)] = 180;
                    else
                        arrR[I(i,j,this->width)] = arrG[I(i,j,this->width)] = arrB[I(i,j,this->width)] = 255;
                else
                    if (((i/8)&1)^((j/8)&1))
                        arrR[I(i,j,this->width)] = arrG[I(i,j,this->width)] = arrB[I(i,j,this->width)] = 90;
                    else
                        arrR[I(i,j,this->width)] = arrG[I(i,j,this->width)] = arrB[I(i,j,this->width)] = 127;
                }
        }

    for (int ii=(int)(this->layers.size() - 1); ii >=-1; --ii) {
        DocumentLayer *curLayer;
        if (ii == -1) {
            if ((this->selectedLayer != -1) && (!forSave))
                curLayer = this->layers[this->selectedLayer];
            else
                continue;
        } else {
            if ((this->selectedLayer != -1) && (ii == this->selectedLayer) && (!forSave))
                continue;
            else
                curLayer = this->layers[ii];
        }
        int ha = curLayer->hOffset, hb = curLayer->height + curLayer->hOffset,
            wa = curLayer->wOffset, wb = curLayer->width + curLayer->wOffset;
        if (ha < 0) ha = 0;
        if (wa < 0) wa = 0;
        if (hb > this->height) hb = this->height;
        if (wb > this->width) wb = this->width;
        for (int i = ha; i < hb; ++i)
            for (int j = wa; j < wb; ++j)
                if (!curLayer->D[I(i-curLayer->hOffset, j-curLayer->wOffset, curLayer->width)]) {
                    arrR[I(i,j,this->width)] = curLayer->R[I(i-curLayer->hOffset, j-curLayer->wOffset, curLayer->width)];
                    arrG[I(i,j,this->width)] = curLayer->G[I(i-curLayer->hOffset, j-curLayer->wOffset, curLayer->width)];
                    arrB[I(i,j,this->width)] = curLayer->B[I(i-curLayer->hOffset, j-curLayer->wOffset, curLayer->width)];
                    if ((ii != -1) && (this->selectedLayer != -1) && (!forSave)) {
                        arrR[I(i,j,this->width)] >>= 1, arrG[I(i,j,this->width)] >>= 1, arrB[I(i,j,this->width)] >>= 1;
                    }
                    if ((ii == -1) && (!forSave)) {
                        if (this->checkInside(i,j))
                            arrR[I(i,j,this->width)] >>= 1, arrG[I(i,j,this->width)] >>= 1, arrB[I(i,j,this->width)] >>= 1;
                    }
                }
    }

    if (!forSave) {
        for (vector<QPoint>::iterator ite = this->choosedPoints.begin();
             ite != this->choosedPoints.end();
             ++ite) {
            for (int j=0; j<9; ++j) {
                int now_x = ite->y() + mark_x[j];
                int now_y = ite->x() + mark_y[j];
                if ((now_x >= 0) && (now_x < this->width) && (now_y >= 0) && (now_y < this->height)) {
                    arrR[I(now_x,now_y,this->width)] = 255;
                    arrG[I(now_x,now_y,this->width)] = 255;
                    arrB[I(now_x,now_y,this->width)] = 255;
                }
            }
        }

        if (this->choosedPoints.size() > 1) {
            for (int i=0; i<this->choosedPoints.size(); ++i) {
                int nex = (i + 1) % this->choosedPoints.size();
                int x1 = this->choosedPoints[i].y(), y1 = this->choosedPoints[i].x();
                int x2 = this->choosedPoints[nex].y(), y2 = this->choosedPoints[nex].x();
                if (((x1 - x2) == 0) && ((y1 - y2) == 0)) continue;
                if (abs(x1-x2) >= abs(y1-y2)) {
                    int xmin = min(x1, x2);
                    int xmax = max(x1, x2);
                    for (int x = xmin; x <= xmax; ++x) {
                        int y = ((y1 - y2) * x + x1 * y2 - x2 * y1) / (x1 - x2);
                        if ((x >= 0) && (x < this->height) && (y >= 0) && (y < this->width)) {
                            arrR[I(x,y,this->width)] = 255;
                            arrG[I(x,y,this->width)] = 255;
                            arrB[I(x,y,this->width)] = 255;
                        }
                    }
                } else {
                    int ymin = min(y1, y2);
                    int ymax = max(y1, y2);
                    for (int y = ymin; y <= ymax; ++y) {
                        int x = ((x1 - x2) * y + x2 * y1 - x1 * y2) / (y1 - y2);
                        if ((x >= 0) && (x < this->height) && (y >= 0) && (y < this->width)) {
                            arrR[I(x,y,this->width)] = 255;
                            arrG[I(x,y,this->width)] = 255;
                            arrB[I(x,y,this->width)] = 255;
                        }
                    }
                }
            }
        }
    }

    QImage *qimg = new QImage(this->width, this->height, QImage::Format_RGB32);
    for (int i=0; i<this->height; ++i)
        for (int j=0; j<this->width; ++j) {
            qimg->setPixelColor(j, i, QColor(arrR[I(i,j,this->width)], arrG[I(i,j,this->width)], arrB[I(i,j,this->width)]));
        }

    QPixmap *qpixmap = new QPixmap(QPixmap::fromImage(*qimg));
    delete[] arrR;
    delete[] arrG;
    delete[] arrB;
    delete qimg;
    return qpixmap;
}

void Document::upLayer(int index) {
    if ((index > 0) && (index < (int)(this->layers.size()))) {
        DocumentLayer *tmp;
        tmp = layers[index - 1];
        layers[index - 1] = layers[index];
        layers[index] = tmp;
        if (selectedLayer == index) selectedLayer = index - 1; else
            if (selectedLayer == index - 1) selectedLayer = index;
    }
    this->theWindow->viewUpdate();
}

void Document::downLayer(int index) {
    if ((index >= 0) && (index < (int)(this->layers.size() - 1))) {
        DocumentLayer *tmp;
        tmp = layers[index + 1];
        layers[index + 1] = layers[index];
        layers[index] = tmp;
        if (selectedLayer == index) selectedLayer = index + 1; else
            if (selectedLayer == index + 1) selectedLayer = index;
    }
    this->theWindow->viewUpdate();
}

void Document::deleteLayer(int index) {
    if ((index >= 0) && (index < (int)(this->layers.size()))) {
        int len = this->layers.size();
        for (int i = index; i < len; ++i)
            this->layers[i] = this->layers[i + 1];
        this->layers.pop_back();
        if (selectedLayer == index) selectedLayer = -1; else
            if (selectedLayer > index) selectedLayer --;
    }
    this->theWindow->viewUpdate();
}

void Document::selectLayer(int index) {
    if ((index >= 0) && (index < (int)(this->layers.size())))
        this->selectedLayer = index;
    this->theWindow->viewUpdate();
}

void Document::crop() {
    if ((selectedLayer != -1) && (this->choosedPoints.size() >= 3)) {
        DocumentLayer *selectLayer = this->layers[this->selectedLayer];
        for (int i=0; i<selectLayer->height; ++i)
            for (int j=0; j<selectLayer->width; ++j)
                if (!selectLayer->D[I(i,j,selectLayer->width)])
                    if (this->checkInside(i+selectLayer->hOffset, j+selectLayer->wOffset))
                        selectLayer->D[I(i,j,selectLayer->width)] = true;
    }
    this->choosedPoints.clear();
}

void Document::cropByImageMask(QImage *maskimg) {
    if (selectedLayer != -1) {
        DocumentLayer *selectLayer = this->layers[this->selectedLayer];
        int lh = selectLayer->height, lw = selectLayer->width, mh = maskimg->height(), mw = maskimg->width();
        int rh = min(lh, mh), rw = min(lw, mw);
        for (int i=0; i<rh; ++i)
            for (int j=0; j<rw; ++j) {
                QRgb maskRgb = maskimg->pixel(j, i);
                if ((qRed(maskRgb) + qGreen(maskRgb) + qBlue(maskRgb)) / 3 > 240) {
                    selectLayer->D[I(i, j, lw)] = true;
                }
            }
    }
}

void Document::unSelect() {
    this->selectedLayer = -1;
    this->theWindow->viewUpdate();
}

void Document::save(string s) {
    QPixmap *img = this->genViewImage(true);
    img->save(QString::fromStdString(s));
    delete img;
}

void Document::addPoint(QPoint newPoint) {
    this->choosedPoints.push_back(newPoint);
}

void Document::deletePoints() {
    this->choosedPoints.clear();
}

bool Document::checkInside(int x, int y) {
    if (this->choosedPoints.size() < 3) return false;
    /*
    double xa = -1.0, ya = (double)y + EPS;
    double xb = (double)x + EPS, yb = (double)y + EPS;
    int crossCnt = 0;
    for (int i=0; i<this->choosedPoints.size(); ++i) {
        int nex = (i + 1) % this->choosedPoints.size();
        double xc = (double)(this->choosedPoints[i].y()), yc = (double)(this->choosedPoints[i].x());
        double xd = (double)(this->choosedPoints[nex].y()), yd = (double)(this->choosedPoints[nex].x());
        double a = cross(xb-xa, yb-ya, xc-xa, yc-ya) * cross(xb-xa, yb-ya, xd-xa, yd-ya);
        double b = cross(xd-xc, yd-yc, xa-xc, ya-yc) * cross(xd-xc, yd-yc, xb-xc, yb-yc);
        if (a*b < -EPS)
            crossCnt++;
    }
    return (crossCnt & 1);
    */

    int j=this->choosedPoints.size() - 1;
    bool oddNodes = false;

    for (int i=0; i<this->choosedPoints.size(); i++) {
        if ((this->choosedPoints[i].x() < y && this->choosedPoints[j].x() >= y
            || this->choosedPoints[j].x() < y && this->choosedPoints[i].x() >=y)
            && (this->choosedPoints[i].y() <= x || this->choosedPoints[j].y() <= x)) {
            oddNodes ^= ((float)this->choosedPoints[i].y()+(float)(y-this->choosedPoints[i].x())/(float)(this->choosedPoints[j].x()-this->choosedPoints[i].x())*(float)(this->choosedPoints[j].y()-this->choosedPoints[i].y())<(float)x);
        }
        j=i;
    }
    return oddNodes;
}

