#ifndef QUADTREENODE_H
#define QUADTREENODE_H

#include <cassert>
#include <iostream>

using namespace std;

class QuadTreeNode {
public:
    int xl, xr, yl, yr;
    int scope;
    QuadTreeNode *a, *b, *c, *d;
    bool isLeaf;

    /** left close, right open interval **/
    QuadTreeNode(int xl, int xr, int yl, int yr);
    ~QuadTreeNode();
    void split();

    QuadTreeNode *findNode(int x, int y);
    QuadTreeNode *findUpperBoundNode(int x, int y);

    static void splitFromRootByPixel(QuadTreeNode *root, int x, int y);

    static void splitFromRoot(QuadTreeNode *root, int x, int y, int scaleLim);
};

#endif // QUADTREENODE_H
