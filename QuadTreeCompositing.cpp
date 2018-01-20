#include "QuadTreeCompositing.h"

#define I(i,j,h) ((i)*(h)+(j))

int d_x[4] = {-1, 0, 0, 1};
int d_y[4] = {0, -1, 1, 0};

void traversePaint(QuadTreeNode *p, int *R, int *G, int *B, int n, int m) {
    if (p->isLeaf) {
        int c_R = rand() & 0xFF, c_G = rand() & 0xFF, c_B = rand() & 0xFF;
        for (int i=p->xl; i<p->xr; ++i)
            if (i < n)
                for (int j=p->yl; j<p->yr; ++j)
                    if (j < m)
                        R[I(i,j,m)] = c_R, G[I(i,j,m)] = c_G, B[I(i,j,m)] = c_B;
    } else {
        traversePaint(p->a, R, G, B, n, m);
        traversePaint(p->b, R, G, B, n, m);
        traversePaint(p->c, R, G, B, n, m);
        traversePaint(p->d, R, G, B, n, m);
    }
}

vector<pair<int, double>> *edgeInterpolate(int i, int j, QuadTreeNode *root, const bool* isKeyPoint, const int *keyPointNo, int, int m) {
    QuadTreeNode *nodes[2] = {root->findNode(i, j), root->findUpperBoundNode(i, j)};
    vector<pair<int, double>> *ans = new vector<pair<int, double>>();
    for (int k=0; k<2; ++k) {
        if ((i == nodes[k]->xl) || (i == nodes[k]->xr))
            if (isKeyPoint[I(i, nodes[k]->yl, m)] && isKeyPoint[I(i, nodes[k]->yr, m)]) {
                ans->push_back(make_pair(keyPointNo[I(i, nodes[k]->yl, m)], double(nodes[k]->yr - j)/double(nodes[k]->yr - nodes[k]->yl)));
                ans->push_back(make_pair(keyPointNo[I(i, nodes[k]->yr, m)], double(j - nodes[k]->yl)/double(nodes[k]->yr - nodes[k]->yl)));
                break;
            }
        if ((j == nodes[k]->yl) || (j == nodes[k]->yr))
            if (isKeyPoint[I(nodes[k]->xl, j, m)] && isKeyPoint[I(nodes[k]->xr, j, m)]) {
                ans->push_back(make_pair(keyPointNo[I(nodes[k]->xl, j, m)], double(nodes[k]->xr - i)/double(nodes[k]->xr - nodes[k]->xl)));
                ans->push_back(make_pair(keyPointNo[I(nodes[k]->xr, j, m)], double(i - nodes[k]->xl)/double(nodes[k]->xr - nodes[k]->xl)));
                break;
            }
    }
    return ans;
}

vector<pair<int, double>>* interpolate(int i, int j, QuadTreeNode *root, const bool* isKeyPoint, const int *keyPointNo, int n, int m) {
    if (isKeyPoint[I(i,j,m)]) {
        vector<pair<int, double>> *ans = new vector<pair<int, double>>();
        ans->push_back(make_pair(keyPointNo[I(i,j,m)], 1.0));
        return ans;
    } else {
        QuadTreeNode *node = root->findNode(i, j);
        map<int, double> points;
        points.clear();
        int x[4], y[4];
        double w[4];
        x[0] = node->xl, y[0] = node->yl, w[0] = double((node->xr - i) * (node->yr - j)) / double(node->scope * node->scope);
        x[1] = node->xl, y[1] = node->yr, w[1] = double((node->xr - i) * (j - node->yl)) / double(node->scope * node->scope);
        x[2] = node->xr, y[2] = node->yl, w[2] = double((i - node->xl) * (node->yr - j)) / double(node->scope * node->scope);
        x[3] = node->xr, y[3] = node->yr, w[3] = double((i - node->xl) * (j - node->yl)) / double(node->scope * node->scope);
        for (int i=0; i<4; ++i) {
            if (w[i] > 1e-6) {
                if (isKeyPoint[I(x[i], y[i], m)]) {
                    if (points.count(keyPointNo[I(x[i], y[i], m)]) == 0)
                        points[keyPointNo[I(x[i], y[i], m)]] = w[i];
                    else
                        points[keyPointNo[I(x[i], y[i], m)]] += w[i];
                } else {
                    vector<pair<int, double>> *edgePoints = edgeInterpolate(x[i], y[i], root, isKeyPoint, keyPointNo, n, m);
                    for (vector<pair<int, double>>::iterator ite = edgePoints->begin(); ite != edgePoints->end(); ++ite) {
                        if (points.count(ite->first) == 0)
                            points[ite->first] = w[i] * ite->second;
                        else
                            points[ite->first] += w[i] * ite->second;
                    }
                    delete edgePoints;
                }
            }
        }
        vector<pair<int, double>> *ans = new vector<pair<int, double>>();
        for (map<int, double>::iterator ite = points.begin(); ite != points.end(); ++ite)
            ans->push_back(make_pair(ite->first, ite->second));
        return ans;
    }
}

int getLatentColor(int channelNo, int x, int y, int discloseZ, Document *doc) {
    int layerN = doc->layers.size();
    for (int i=0; i<layerN; ++i)
        if (i != discloseZ) {
            DocumentLayer *nowLayer = doc->layers[i];
            if ((x - nowLayer->hOffset >= 0) && (x - nowLayer->hOffset < nowLayer->height) &&
                (y - nowLayer->wOffset >= 0) && (y - nowLayer->wOffset < nowLayer->width)) {
                if (!nowLayer->D[I(x - nowLayer->hOffset, y - nowLayer->wOffset, nowLayer->width)]) {
                    if (channelNo == 0)
                        return nowLayer->R[I(x - nowLayer->hOffset, y - nowLayer->wOffset, nowLayer->width)];
                    if (channelNo == 1)
                        return nowLayer->G[I(x - nowLayer->hOffset, y - nowLayer->wOffset, nowLayer->width)];
                    if (channelNo == 2)
                        return nowLayer->B[I(x - nowLayer->hOffset, y - nowLayer->wOffset, nowLayer->width)];
                }
            }
        }
    // white background by default
    return 255;
}

Eigen::SparseMatrix<double> *getMatrix(vector<vector<pair<int, double>>*> *A, int n) {
    int m = A->size();
    Eigen::SparseMatrix<double> *mat = new Eigen::SparseMatrix<double>(n, n);
    mat->setZero();
    int ti, tj;
    double tv;
    for (int i=0; i<m; ++i) {
        vector<pair<int, double>> *nowVec = (*A)[i];
        for (int j=0; j<int(nowVec->size()); ++j)
            for (int k=0; k<int(nowVec->size()); ++k) {
                ti = (*nowVec)[j].first;
                tj = (*nowVec)[k].first;
                tv = (*nowVec)[j].second * (*nowVec)[k].second;
                Eigen::SparseMatrix<double>::Scalar &scalar = mat->coeffRef(ti, tj);
                scalar += tv;
            }
    }
    mat->makeCompressed();
    return mat;
}

Eigen::SparseVector<double> *getBVector(vector<vector<pair<int, double>>*> *A, double *b, int n) {
    int m = A->size();
    Eigen::SparseVector<double> *vec = new Eigen::SparseVector<double>(n);
    vec->setZero();
    for (int i=0; i<m; ++i) {
        vector<pair<int, double>>* nowVec = (*A)[i];
        for (vector<pair<int, double>>::iterator ite = nowVec->begin(); ite != nowVec->end(); ++ite) {
            Eigen::SparseVector<double>::Scalar &scalar = vec->coeffRef(ite->first);
            scalar += ite->second * b[i];
        }
    }
    return vec;
}

void QuadTreeCompositing::run() {
    assert(doc != NULL);
    stringstream logger;
    logger << "Begin";
    this->startTime = time(0);
    this->printLog(logger);

    srand(time(0));

    int m = doc->width, n = doc->height, nm = doc->width * doc->height;
    int *R = new int[nm], *G = new int[nm], *B = new int[nm], *z = new int[nm];
    int *nR = new int[nm], *nG = new int[nm], *nB = new int[nm];
    int *qtreeR = new int[nm], *qtreeG = new int[nm], *qtreeB = new int[nm];

    vector<int> borderPointXs, borderPointYs;

    int scale = max(n, m) - 1;
    QuadTreeNode *root = NULL;

    bool *isKeyPoint = new bool[nm];
    int *keyPointNo = new int[nm];
    int totKeyPoint = 0;

    vector<pair<int, double>> **interpolates = new vector<pair<int, double>>*[nm];

    vector<vector<pair<int, double>>*> A;
    Eigen::SparseMatrix<double> *matA = NULL;

    Eigen::BiCGSTAB<Eigen::SparseMatrix<double>, Eigen::IdentityPreconditioner> solver;

    double *b = NULL;

    double *delta = new double[nm];
    double *x = NULL;

    /* STAGE 1: calculating x0 - simple covering compositing */
    {
        for (int i=0; i<nm; ++i)
            R[i] = G[i] = B[i] = 255, z[i] = -1;

        for (int k=0; k<int(doc->layers.size()); ++k) {
            DocumentLayer *curLayer = doc->layers[k];
            int ha = curLayer->hOffset, hb = curLayer->height + curLayer->hOffset;
            int wa = curLayer->wOffset, wb = curLayer->width + curLayer->wOffset;
            if (ha < 0) ha = 0;
            if (wa < 0) wa = 0;
            if (hb > n) hb = n;
            if (wb > m) wb = m;
            for (int i = ha; i<hb; ++i)
                for (int j=wa; j<wb; ++j)
                    if (z[I(i,j,m)] == -1)
                        if (!curLayer->D[I(i-curLayer->hOffset, j-curLayer->wOffset, curLayer->width)]) {
                            R[I(i,j,m)] = curLayer->R[I(i-curLayer->hOffset, j-curLayer->wOffset, curLayer->width)];
                            G[I(i,j,m)] = curLayer->G[I(i-curLayer->hOffset, j-curLayer->wOffset, curLayer->width)];
                            B[I(i,j,m)] = curLayer->B[I(i-curLayer->hOffset, j-curLayer->wOffset, curLayer->width)];
                            z[I(i,j,m)] = k;
                        }
        }
    }
    logger << "Initial image compositing finish";
    this->printLog(logger);

    /* STAGE 2: calculating border points */
    {
        borderPointXs.clear();
        borderPointYs.clear();
        for (int i=0; i<n; ++i)
            for (int j=0; j<m; ++j) {
                bool border = false;
                for (int k=0; k<4; ++k) {
                    int tx = i + d_x[k], ty = j + d_y[k];
                    if ((tx >= 0) && (tx < n) && (ty >= 0) && (ty < m))
                        if (z[I(i,j,m)] != z[I(tx,ty,m)]) {
                            border = true;
                            break;
                        }
                }
                if (border) {
                    borderPointXs.push_back(i);
                    borderPointYs.push_back(j);
                }
            }
    }
    logger << "Border pixels figure out";
    this->printLog(logger);

    /* STAGE 3: building quad tree */
    {
        int tmp = 1;
        while (tmp <= scale) tmp <<= 1;
        scale = tmp;

        root = new QuadTreeNode(0, scale, 0, scale);

        for (int i=0; i<int(borderPointXs.size()); ++i)
            QuadTreeNode::splitFromRootByPixel(root, borderPointXs[i], borderPointYs[i]);

        // constrain image border pixels to be border pixels in quadtree
        for (int i=0; i<m; ++i)
            QuadTreeNode::splitFromRootByPixel(root, n-1, i);
        for (int i=0; i<n; ++i)
            QuadTreeNode::splitFromRootByPixel(root, i, m-1);
    }
    logger << "Quad tree building finish";
    this->printLog(logger);

    /* STAGE 4: detect key points */
    {
        memset(isKeyPoint, 0, sizeof(bool) * nm);
        memset(keyPointNo, 0, sizeof(int) * nm);
        for (int i=0; i<n; ++i)
            for (int j=0; j<m; ++j) {
                QuadTreeNode *theNode = root->findNode(i,j);
                assert(theNode != NULL);
                if ((theNode->xl == i) && (theNode->yl == j)) {
                    QuadTreeNode *neighborNode = root->findUpperBoundNode(i, j);
                    assert(neighborNode != NULL);
                    if (((neighborNode->xr == i) && (neighborNode->yr == j)) || (theNode->xl == 0) || (theNode->yl == 0))
                        isKeyPoint[I(i,j,m)] = true, keyPointNo[I(i,j,m)] = totKeyPoint, ++totKeyPoint;
                }
            }
    }

    /*// B DEBUG
    for (int i=0; i<n; ++i)
        for (int j=0; j<m; ++j)
            if (isKeyPoint[I(i,j,m)])
                cout << "NO " << keyPointNo[I(i,j,m)] << " X" << i << " Y" << j << endl;
    // E DEBUG*/

    /* STAGE 5: interpolate */
    for (int i=0; i<n; ++i)
        for (int j=0; j<m; ++j)
            interpolates[I(i,j,m)] = interpolate(i, j, root, isKeyPoint, keyPointNo, n, m);
    logger << "Interpolation finish";
    this->printLog(logger);

    /* STAGE 6: construct sparse matrix A */
    {
        vector<pair<int, double>> *intera, *interb, *interc;
        for (int i=0; i<n; ++i)
            for (int j=0; j<m; ++j)
                for (int tt = 0; tt < 2; ++tt) {
                    intera = interpolates[I(i,j,m)];
                    if (tt == 0) {
                        if (i == 0) continue;
                        interb = interpolates[I(i-1,j,m)];
                    }
                    if (tt == 1) {
                        if (j == 0) continue;
                        interb = interpolates[I(i,j-1,m)];
                    }
                    interc = new vector<pair<int, double>>();
                    for (int k=0; k<int(intera->size()); ++k)
                        interc->push_back((*intera)[k]);
                    for (int k=0; k<int(interb->size()); ++k) {
                        int index = -1;
                        for (int l=0; l<interc->size(); ++l)
                            if ((*interc)[l].first == (*interb)[k].first) {
                                index = l;
                                (*interc)[l].second -= (*interb)[k].second;
                                break;
                            }
                        if (index == -1)
                            interc->push_back(make_pair((*interb)[k].first, -(*interb)[k].second));
                    }
                    A.push_back(interc);
                }
        // constrain last pixel to be 0, to prevent incomplete-constrained
        vector<pair<int, double>> *lastRow = new vector<pair<int, double>>();
        vector<pair<int, double>> *proto = interpolates[I(n-1,m-1,m)];
        for (int i=0; i<int(proto->size()); ++i)
            lastRow->push_back((*proto)[i]);
        A.push_back(lastRow);
        // calculate A^T * A
        matA = getMatrix(&A, totKeyPoint);
    }
    logger << "Matrix A construction finish";
    this->printLog(logger);

    logger << "Constraints: " << A.size();
    this->printLog(logger);
    logger << "Variables: " << totKeyPoint;
    this->printLog(logger);

    /* STAGE 7: matrix A, vector b, vector x initialize */
    {
        solver.compute(*matA);
        b = new double[A.size()];
        x = new double[totKeyPoint];
    }

    /* STAGE 8: solve Ax=b for the 3 channels respectively */
    for (int channel = 0; channel < 3; ++channel) {
        if (channel == 0)
            logger << "Calculating channel R";
        if (channel == 1)
            logger << "Calculating channel G";
        if (channel == 2)
            logger << "Calculating channel B";
        this->printLog(logger);
        int *color = NULL, *nColor = NULL;
        if (channel == 0) color = R, nColor = nR;
        if (channel == 1) color = G, nColor = nG;
        if (channel == 2) color = B, nColor = nB;
        int ni, nj, bp = 0;
        int s, minz, t;
        for (int i=0; i<n; ++i)
            for (int j=0; j<m; ++j)
                for (int tt=0; tt<2; ++tt) {
                    if (tt == 0) {
                        if (i == 0) continue;
                        ni = i-1, nj = j;
                    }
                    if (tt == 1) {
                        if (j == 0) continue;
                        ni = i, nj = j-1;
                    }
                    if (z[I(i,j,m)] == z[I(ni,nj,m)])
                        b[bp++] = 0.0;
                    else {

                        if (z[I(i,j,m)] == -1)
                            minz = z[I(ni,nj,m)];
                        else if (z[I(ni,nj,m)] == -1)
                            minz = z[I(i,j,m)];
                        else
                            minz = min(z[I(i,j,m)], z[I(ni,nj,m)]);
                        s = color[I(i,j,m)] - color[I(ni,nj,m)];
                        t = getLatentColor(channel, i, j, minz, doc) - getLatentColor(channel, ni, nj, minz, doc);
                        b[bp++] = t-s;
                    }
                }
        b[bp++] = 0.0;
        Eigen::SparseVector<double> *vecB = getBVector(&A, b, totKeyPoint);

        /*// B DEBUG
        cout << endl << "A:" << endl;
        for (int i=0; i<A.size(); ++i) {
            vector<pair<int, double>> *tmp = A[i];
            cout << "line " << i << ": ";
            for (int j=0; j<tmp->size(); ++j)
                cout << "[" << (*tmp)[j].first << "]=" << (*tmp)[j].second << " ";
            cout << endl;
        }
        cout << endl << "b:" << endl;
        for (int i=0; i<b.size(); ++i) {
            cout << "[" << i << "]=" << b[i] << endl;
        }
        cout << endl << "A^T*A:" << endl << (*matA) << endl;
        cout << endl << "A^T*b:" << endl << (*vecB) << endl;
        // E DEBUG*/

        Eigen::SparseVector<double> vecAns = solver.solve(*vecB);
        delete vecB;
        logger << "'x' solved out";
        this->printLog(logger);

        for (int i=0; i<totKeyPoint; ++i)
            x[i] = 0.0;
        for (Eigen::SparseVector<double>::InnerIterator it(vecAns); it; ++it)
            x[it.index()] = it.value();

        /*// B DEBUG
        cout << endl << "x:" << endl;
        for (int i=0; i<totKeyPoint; ++i) {
            cout << "[" << i << "]=" << keyV[i] << endl;
        }
        cout << endl;
        // E DEBUG*/

        double deltaAvg = 0.0;
        for (int i=0; i<n; ++i)
            for (int j=0; j<m; ++j) {
                delta[I(i,j,m)] = 0.0;
                vector<pair<int, double>> *inter = interpolates[I(i,j,m)];
                for (vector<pair<int, double>>::iterator ite = inter->begin(); ite != inter->end(); ++ite)
                    delta[I(i,j,m)] += x[ite->first] * ite->second;
                deltaAvg += delta[I(i,j,m)];
            }
        deltaAvg /= double(nm);
        logger << "'x' average: " << deltaAvg;
        this->printLog(logger);
        int itmp;
        for (int i=0; i<n; ++i)
            for (int j=0; j<m; ++j) {
                itmp = (int)((double)color[I(i,j,m)] + delta[I(i,j,m)] - deltaAvg + 0.5);
                if (itmp > 255) itmp = 255;
                if (itmp < 0) itmp = 0;
                nColor[I(i,j,m)] = itmp;
            }
    }

    /* STAGE 9: render the answer */
    {
        traversePaint(root, qtreeR, qtreeG, qtreeB, n, m);
        doc->addLayerAtBottom(new DocumentLayer(qtreeR, qtreeG, qtreeB, m, n, "QuadTree"));
        doc->addLayerAtBottom(new DocumentLayer(nR, nG, nB, m, n, "Result"));
    }

    /* Clean */
    delete[] R;
    delete[] G;
    delete[] B;
    delete[] z;
    delete[] nR;
    delete[] nG;
    delete[] nB;
    delete[] qtreeR;
    delete[] qtreeG;
    delete[] qtreeB;
    delete root;
    delete[] isKeyPoint;
    delete[] keyPointNo;
    for (int i=0; i<nm; ++i) {
        vector<pair<int, double>> *now = interpolates[i];
        delete now;
    }
    delete[] interpolates;
    for (int i=0; i<int(A.size()); ++i)
        delete A[i];
    delete matA;
    delete[] b;
    delete[] delta;
    delete[] x;

    ans = doc;

    logger << "Finish";
    this->printLog(logger);
}

void QuadTreeCompositing::printLog(stringstream &logger) {
    logger.flush();
    emit this->updateStatus(QString("[") + QString::number((int)(time(0) - this->startTime)) + QString(" s] ") + QString::fromStdString(logger.str()));
    logger.str("");
}
