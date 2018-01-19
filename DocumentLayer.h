#ifndef DOCUMENTLAYER_H
#define DOCUMENTLAYER_H

#include <string>
#include <QImage>

#define Q(i,j,h) ((i)*(h)+(j))

using namespace std;

class DocumentLayer
{
public:
    DocumentLayer(QImage *qimg, string name, int wOffset = 0, int hOffset = 0);
    DocumentLayer(int *R, int *G, int *B, int width, int height, string name, int wOffset = 0, int hOffset = 0);
    ~DocumentLayer();
    void addOffset(int w_offset, int h_offset);

    string name;
    int width, height, wOffset, hOffset;
    int *R, *G, *B;
    bool *D;
};

#endif // DOCUMENTLAYER_H
