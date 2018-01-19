#include "DocumentLayer.h"

DocumentLayer::DocumentLayer(QImage *qimg, string name, int wOffset, int hOffset) {
    int h, w;
    this->name = name;
    this->width = w = qimg->width();
    this->height = h = qimg->height();
    this->wOffset = wOffset;
    this->hOffset = hOffset;

    this->R = new int[this->width * this->height];
    this->G = new int[this->width * this->height];
    this->B = new int[this->width * this->height];
    for (int i=0; i<h; ++i)
        for (int j=0; j<w; ++j) {
            this->R[Q(i,j,w)] = qRed(qimg->pixel(j, i));
            this->G[Q(i,j,w)] = qGreen(qimg->pixel(j, i));
            this->B[Q(i,j,w)] = qBlue(qimg->pixel(j, i));
        }

    this->D = new bool[this->width * this->height];
    memset(D, 0, sizeof(bool) * this->width * this->height);
}

DocumentLayer::DocumentLayer(int *R, int *G, int *B, int width, int height, string name, int wOffset, int hOffset) {
    this->name = name;
    this->width = width;
    this->height = height;
    this->wOffset = wOffset;
    this->hOffset = hOffset;

    this->R = new int[width * height];
    this->G = new int[width * height];
    this->B = new int[width * height];
    for (int i=0; i<width*height; ++i) {
        this->R[i] = R[i];
        this->G[i] = G[i];
        this->B[i] = B[i];
    }
    this->D = new bool[this->width * this->height];
    memset(D, 0, sizeof(bool) * this->width * this->height);
}

void DocumentLayer::addOffset(int w_offset, int h_offset) {
    this->wOffset += w_offset;
    this->hOffset += h_offset;
}

DocumentLayer::~DocumentLayer() {
    delete[] this->R;
    delete[] this->G;
    delete[] this->B;
    delete[] this->D;
}
