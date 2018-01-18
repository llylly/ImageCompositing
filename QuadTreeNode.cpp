#include "QuadTreeNode.h"

QuadTreeNode::QuadTreeNode(int xl, int xr, int yl, int yr) {
    this->xl = xl, this->xr = xr, this->yl = yl, this->yr = yr;
    this->a = this->b = this->c = this->d = NULL;
    this->scope = xr-xl;
    this->isLeaf = true;
}

QuadTreeNode::~QuadTreeNode() {
    if (!this->isLeaf) {
        delete a, b, c, d;
    }
}

void QuadTreeNode::split() {
    if (scope == 1) return;
    if (this->isLeaf) {
        this->isLeaf = false;
        int xmid = (xl + xr) >> 1, ymid = (yl + yr) >> 1;
        a = new QuadTreeNode(xl, xmid, yl, ymid);
        b = new QuadTreeNode(xl, xmid, ymid, yr);
        c = new QuadTreeNode(xmid, xr, yl, ymid);
        d = new QuadTreeNode(xmid, xr, ymid, yr);
    }
}

QuadTreeNode *QuadTreeNode::findNode(int x, int y) {
    if ((xl <= x) && (x < xr) && (yl <= y) && (y < yr)) {
        if (this->isLeaf)
            return this;
        else {
            int xmid = (this->xl + this->xr) >> 1, ymid = (this->yl + this->yr) >> 1;
            if (x < xmid)
                if (y < ymid)
                    return this->a->findNode(x, y);
                else
                    return this->b->findNode(x, y);
            else
                if (y < ymid)
                    return this->c->findNode(x, y);
                else
                    return this->d->findNode(x, y);
        }
    } else
        return NULL;
}

QuadTreeNode *QuadTreeNode::findUpperBoundNode(int x, int y) {
    if ((xl <= x) && (x <= xr) && (yl <= y) && (y <= yr)) {
        if (this->isLeaf)
            return this;
        else {
            int xmid = (this->xl + this->xr) >> 1, ymid = (this->yl + this->yr) >> 1;
            if (x <= xmid)
                if (y <= ymid)
                    return this->a->findUpperBoundNode(x, y);
                else
                    return this->b->findUpperBoundNode(x, y);
            else
                if (y <= ymid)
                    return this->c->findUpperBoundNode(x, y);
                else
                    return this->d->findUpperBoundNode(x, y);
        }
    } else
        return NULL;
}

void QuadTreeNode::splitFromRootByPixel(QuadTreeNode *root, int x, int y) {
    if ((x-1 >= root->xl) && (y-1 >= root->yl)) {
        QuadTreeNode *node = root->findNode(x-1, y-1);
        assert(node != NULL);
        if (node->scope > 1)
            QuadTreeNode::splitFromRoot(root, x-1, y-1, 1);
    }
    if ((x < root->xr) && (y < root->yr)) {
        QuadTreeNode *node = root->findNode(x, y);
        assert(node != NULL);
        if (node->scope > 1)
            QuadTreeNode::splitFromRoot(root, x, y, 1);
    }
}

void QuadTreeNode::splitFromRoot(QuadTreeNode *root, int x, int y, int scaleLim) {
    QuadTreeNode *p = root;
    while (p->scope > scaleLim) {
        if (p->isLeaf) {
            p->split();

            QuadTreeNode *neighbor;
            if (p->xl > root->xl) {
                neighbor = root->findNode(p->xl-1, p->yl);
                if (neighbor->scope > p->scope) {
                    QuadTreeNode::splitFromRoot(root, p->xl-1, p->yl, p->scope);
                }
            }
            if (p->xr < root->xr) {
                neighbor = root->findNode(p->xr, p->yl);
                if (neighbor->scope > p->scope) {
                    QuadTreeNode::splitFromRoot(root, p->xr, p->yl, p->scope);
                }
            }
            if (p->yl > root->yl) {
                neighbor = root->findNode(p->xl, p->yl-1);
                if (neighbor->scope > p->scope) {
                    QuadTreeNode::splitFromRoot(root, p->xl, p->yl-1, p->scope);
                }
            }
            if (p->yr < root->yr) {
                neighbor = root->findNode(p->xl, p->yr);
                if (neighbor->scope > p->scope) {
                    QuadTreeNode::splitFromRoot(root, p->xl, p->yr, p->scope);
                }
            }
        }

        int xmid = (p->xl + p->xr) >> 1, ymid = (p->yl + p->yr) >> 1;
        if (x < xmid)
            if (y < ymid)
                p = p->a;
            else
                p = p->b;
        else
            if (y < ymid)
                p = p->c;
            else
                p = p->d;
    }
}

