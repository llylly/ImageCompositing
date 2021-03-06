#include "MVCCompositingThread.h"

void drawCDT(CDT &cdt, int *R, int *G, int *B, bool *D, int n, int m) {
    for (CDT::Finite_edges_iterator edge_ite = cdt.finite_edges_begin(); edge_ite != cdt.finite_edges_end(); ++edge_ite) {
//    for (CDT::Finite_faces_iterator face_ite = cdt.finite_faces_begin(); face_ite != cdt.finite_faces_end(); ++face_ite) {
//        for (int i=0; i<3; ++i) {
        Point source = cdt.segment(edge_ite).source();
        Point dest = cdt.segment(edge_ite).target();
//        Point source = face_ite->vertex(i)->point();
//        Point dest = face_ite->vertex((i+1)%3)->point();
        int x1 = int(source.x() / 2. + 0.5), y1 = int(source.y() / 2. + 0.5), x2 = int(dest.x() / 2. + 0.5), y2 = int(dest.y() / 2. + 0.5);
        x1 = max(min(x1, n-1), 0);
        x2 = max(min(x2, n-1), 0);
        y1 = max(min(y1, m-1), 0);
        y2 = max(min(y2, m-1), 0);
        if (abs(x1-x2) >= abs(y1-y2)) {
            int xmin = min(x1, x2);
            int xmax = max(x1, x2);
            if (xmin == xmax) continue;
            for (int x = xmin; x <= xmax; ++x) {
                int y = ((y1 - y2) * x + x1 * y2 - x2 * y1) / (x1 - x2);
                if ((x >= 0) && (x < n) && (y >= 0) && (y < m)) {
                    R[Q(x,y,m)] = 255;
                    G[Q(x,y,m)] = 255;
                    B[Q(x,y,m)] = 255;
                    D[Q(x,y,m)] = false;
                }
            }
        } else {
            int ymin = min(y1, y2);
            int ymax = max(y1, y2);
            if (ymin == ymax) continue;
            for (int y = ymin; y <= ymax; ++y) {
                int x = ((x1 - x2) * y + x2 * y1 - x1 * y2) / (y1 - y2);
                if ((x >= 0) && (x < n) && (y >= 0) && (y < m)) {
                    R[Q(x,y,m)] = 255;
                    G[Q(x,y,m)] = 255;
                    B[Q(x,y,m)] = 255;
                    D[Q(x,y,m)] = false;
                }
            }
        }
//        }
    }
}

int traverse(int *z, int *traverseVisit, int &traverseCnt,
              int s_x, int s_y, int z_lim,
              vector<Point> &borderPoints, vector<Point> &vitalPoints,
              int n, int m) {
    ++traverseCnt;
    int x = s_x, y = s_y;
    int maxColor = -1;
    while (true) {
        traverseVisit[Q(x, y, 2*m+1)] = traverseCnt;
        borderPoints.push_back(Point(x, y));
        int x1, y1, x2, y2;
        if (x & 1)
            x1 = x2 = (x - 1) >> 1;
        else
            x1 = (x >> 1) - 1, x2 = x >> 1;
        if (y & 1)
            y1 = y2 = (y - 1) >> 1;
        else
            y1 = (y >> 1) - 1, y2 = y >> 1;
        if ((x1 >= 0) && (x1 < n) && (y1 >= 0) && (y1 < m))
            if ((x2 >= 0) && (x2 < n) && (y2 >= 0) && (y2 < m))
                if (((z[Q(x1,y1,m)]<=z_lim) && (z[Q(x2,y2,m)]>z_lim)) ||
                    ((z[Q(x1,y1,m)]>z_lim) && (z[Q(x2,y2,m)]<=z_lim))) {
                    vitalPoints.push_back(Point(x, y));
                    maxColor = max(maxColor, max(z[Q(x1,y1,m)], z[Q(x2,y2,m)]));
                }

        int nex_x = -1, nex_y = -1;
        vector<Point> series[2]; bool valid[2] = {true, true};
        series[0].clear(), series[1].clear();
        if (x & 1) {
            series[0].push_back(Point(x-1, y-1)), series[0].push_back(Point(x-2, y)), series[0].push_back(Point(x-1, y+1));
            series[1].push_back(Point(x+1, y-1)), series[1].push_back(Point(x+2, y)), series[1].push_back(Point(x+1, y+1));
        }
        if (y & 1) {
            series[0].push_back(Point(x, y-2)), series[0].push_back(Point(x-1, y-1)), series[0].push_back(Point(x+1, y-1));
            series[1].push_back(Point(x, y+2)), series[1].push_back(Point(x-1, y+1)), series[1].push_back(Point(x+1, y+1));
        }
        for (int i=0; i<2; ++i)
            for (int j=0; j<3; ++j)
                if ((series[i][j].x() >= 0) && (series[i][j].x() <= 2*n) &&
                    (series[i][j].y() >= 0) && (series[i][j].y() <= 2*m)) {
                    int t_x = series[i][j].x(), t_y = series[i][j].y();
                    if (traverseVisit[Q(t_x, t_y, 2*m+1)] == traverseCnt)
                        valid[i] = false;
                }
        for (int i=0; i<2; ++i) {
            if (valid[i])
                for (int j=0; j<3; ++j)
                    if ((series[i][j].x() >= 0) && (series[i][j].x() <= 2*n) &&
                        (series[i][j].y() >= 0) && (series[i][j].y() <= 2*m)) {
                        int n_x = series[i][j].x(), n_y = series[i][j].y();
                        int x1, y1, x2, y2;
                        if (n_x & 1)
                            x1 = x2 = (n_x - 1) >> 1;
                        else
                            x1 = (n_x >> 1) - 1, x2 = n_x >> 1;
                        if (n_y & 1)
                            y1 = y2 = (n_y - 1) >> 1;
                        else
                            y1 = (n_y >> 1) - 1, y2 = n_y >> 1;
                        bool has_out = false, has_in = false;
                        for (int k=0; k<2; ++k) {
                            int nn_x, nn_y;
                            if (k == 0) nn_x = x1, nn_y = y1; else nn_x = x2, nn_y = y2;
                            if ((nn_x >= 0) && (nn_x < n) && (nn_y >= 0) && (nn_y < m)) {
                                if ((z[Q(nn_x, nn_y, m)] >= 0) && (z[Q(nn_x, nn_y, m)] <= z_lim))
                                    has_in = true;
                                else {
                                    has_out = true;
                                    maxColor = max(maxColor, z[Q(nn_x, nn_y, m)]);
                                }
                            } else {
                                has_out = true;
                            }
                        }
                        if (has_in && has_out) {
                            nex_x = n_x, nex_y = n_y;
                            break;
                        }
                        if ((nex_x >= 0) && (nex_y >= 0)) break;
                    }
            if ((nex_x >= 0) && (nex_y >= 0)) break;
        }
        if ((nex_x >= 0) && (nex_y >= 0)) {
            x = nex_x, y = nex_y;
        } else {
            break;
        }
    }
    return maxColor;
}

void floodfill(int *z, int *floodfillVisit, int &floodfillCnt,
               int s_x, int s_y, int z_lim, int *que_x, int *que_y,
               int n, int m) {
    ++floodfillCnt;
    int l=0, r=1;
    que_x[0] = s_x, que_y[0] = s_y;
    floodfillVisit[Q(s_x, s_y, m)] = floodfillCnt;

    int d_x[4] = {-1, 0, 0, 1};
    int d_y[4] = {0, -1, 1, 0};

    while (l < r) {
        for (int k=0; k<4; ++k) {
            int nx = que_x[l] + d_x[k], ny = que_y[l] + d_y[k];
            if ((nx >= 0) && (nx < n) && (ny >= 0) && (ny < m))
                if (floodfillVisit[Q(nx, ny, m)] != floodfillCnt)
                    if ((z[Q(nx, ny, m)] >= 0) && (z[Q(nx, ny, m)] <= z_lim)) {
                        que_x[r] = nx, que_y[r] = ny;
                        floodfillVisit[Q(nx, ny, m)] = floodfillCnt;
                        ++r;
                    }
        }
        ++l;
    }
}

double distance(const Point &a, const Point &b) {
    return (sqrt((a.x() - b.x()) * (a.x() - b.x()) + (a.y() - b.y()) * (a.y() - b.y())));
}

double calcArea(double x1, double y1, double x2, double y2, double x3, double y3) {
    double ax = x2-x1, ay = y2-y1;
    double bx = x3-x1, by = y3-y1;
    return ax * by - ay * bx;
}

double triangleInterpolate(double px, double py, Face_handle &face, double *mesh_delta) {
    double x[3], y[3], d[3], w[3];
    for (int i=0; i<3; ++i) {
        Vertex_handle v = face->vertex(i);
        x[i] = v->point().x(), y[i] = v->point().y();
        d[i] = mesh_delta[v->info() - 1];
    }
    for (int i=0; i<3; ++i) {
        w[i] = calcArea(px, py, x[(i+1)%3], y[(i+1)%3], x[(i+2)%3], y[(i+2)%3]);
    }
    double w_sum = w[0] + w[1] + w[2];
    double ans = 0.0;
    for (int i=0; i<3; ++i)
        ans += w[i] * d[i];
    ans /= w_sum;
    return ans;
}

void MVCCompositingThread::run() {
    assert(doc != NULL);
    stringstream logger;
    logger << "Begin";
    this->startTime = time(0);
    this->printLog(logger);

    int m = doc->width, n = doc->height, nm = doc->width * doc->height;

    int *oR = new int[nm], *oG = new int[nm], *oB = new int[nm], *z = new int[nm];
    double *dR = new double[nm], *dG = new double[nm], *dB = new double[nm];
    int *newR = new int[nm], *newG = new int[nm], *newB = new int[nm];
    int *zPixelCount = new int[doc->layers.size()];
    int *traverseVisit = new int[(2*n+1)*(2*m+1)];
    int *floodfillVisit = new int[nm];
    int *meshR = new int[nm], *meshG = new int[nm], *meshB = new int[nm];
    bool *meshD = new bool[nm];
    int *que_x = new int[nm], *que_y = new int[nm];

    int traverseCnt = 0;
    int floodfillCnt = 0;

    for (int i=0; i<(2*n+1)*(2*m+1); ++i)
        traverseVisit[i] = traverseCnt;
    for (int i=0; i<nm; ++i)
        floodfillVisit[i] = floodfillCnt;

    for (int i=0; i<nm; ++i) {
        oR[i] = oG[i] = oB[i] = 255;
        z[i] = -1;
        meshR[i] = meshG[i] = meshB[i] = 0;
        meshD[i] = true;
    }

    for (int ii=0; ii < (int)(doc->layers.size()); ++ii)
        zPixelCount[ii] = 0;

    for (int ii = (int)(doc->layers.size() - 1); ii >= 0; --ii) {
        DocumentLayer *curLayer = doc->layers[ii];
        int ha = curLayer->hOffset, hb = curLayer->height + curLayer->hOffset,
            wa = curLayer->wOffset, wb = curLayer->width + curLayer->wOffset;
        if (ha < 0) ha = 0;
        if (wa < 0) wa = 0;
        if (hb > n) hb = n;
        if (wb > m) wb = m;
        for (int i=ha; i<hb; ++i)
            for (int j=wa; j<wb; ++j)
                if (!curLayer->D[Q(i-curLayer->hOffset, j-curLayer->wOffset, curLayer->width)]) {
                    oR[Q(i,j,m)] = curLayer->R[Q(i-curLayer->hOffset, j-curLayer->wOffset, curLayer->width)];
                    oG[Q(i,j,m)] = curLayer->G[Q(i-curLayer->hOffset, j-curLayer->wOffset, curLayer->width)];
                    oB[Q(i,j,m)] = curLayer->B[Q(i-curLayer->hOffset, j-curLayer->wOffset, curLayer->width)];
                    if (z[Q(i,j,m)] != -1) --zPixelCount[z[Q(i,j,m)]];
                    z[Q(i,j,m)] = ii;
                    ++zPixelCount[ii];
                }
    }
    for (int i=0; i<nm; ++i)
        dR[i] = oR[i], dG[i] = oG[i], dB[i] = oB[i];

    logger << "Initial image render finish";
    this->printLog(logger);

    int d_x[] = {-1, 0, 0, 1};
    int d_y[] = {0, -1, 1, 0};

    for (int ii=0; ii < (int)(doc->layers.size() - 1); ++ii) {

        int tryTimes = 0;
        int s_x, s_y, inner_x, inner_y;

        while ((zPixelCount[ii] > 0) && (tryTimes < 100)) {
            logger << "Iteration for layer #" << ii << " block #" << tryTimes;
            this->printLog(logger);

            s_x = -1, s_y = -1, inner_x = -1, inner_y = -1;
            for (int i=0; i<n; ++i) {
                for (int j=0; j<m; ++j)
                    if ((z[Q(i,j,m)] != -1) && (z[Q(i,j,m)] <= ii)) {
                        for (int k=0; k<4; ++k) {
                            int ni = i + d_x[k], nj = j + d_y[k];
                            if ((ni >= 0) && (ni < n) && (nj >= 0) && (nj < m)) {
                                if (z[Q(ni,nj,m)] > ii) {
                                    s_x = i + ni + 1, s_y = j + nj + 1;
                                    inner_x = i, inner_y = j;
                                    break;
                                }
                            }
                        }
                        if ((s_x > 0) && (s_y > 0)) break;
                    }
                if ((s_x > 0) && (s_y > 0)) break;
            }
            if ((s_x == -1) && (s_y == -1)) break;

            vector<Point> borderPoints;
            vector<Point> vitalPoints;
            vector<Point> meshPoints;
            borderPoints.clear();
            vitalPoints.clear();
            meshPoints.clear();

            int zToFill = traverse(z, traverseVisit, traverseCnt, s_x, s_y, ii, borderPoints, vitalPoints, n, m);
            floodfill(z, floodfillVisit, floodfillCnt, inner_x, inner_y, ii, que_x, que_y, n, m);

            logger << "  Block and boundary traverse finish";
            this->printLog(logger);

            CDT cdt;
            for (int i=0; i<int(borderPoints.size()-1); ++i)
                cdt.insert_constraint(borderPoints[i], borderPoints[i+1]);
            try {
                cdt.insert_constraint(borderPoints[borderPoints.size()-1], borderPoints[0]);
            } catch (exception e) {
                break;
            }

            Mesher mesher(cdt);
            mesher.set_criteria(Criteria(0.125));
            try {
                mesher.refine_mesh();
            } catch (exception e) {
                break;
            }

            drawCDT(cdt, meshR, meshG, meshB, meshD, n, m);

            logger << "  Adaptive mesh generation finish";
            this->printLog(logger);

            int nv = vitalPoints.size(), nmesh = cdt.number_of_vertices();
            double *vitalDelta = new double[nv];
            double *meshDelta = new double[nmesh];

            logger << "  # of boundary pixel: " << borderPoints.size();
            this->printLog(logger);
            logger << "  # of constrained pixel: " << nv;
            this->printLog(logger);
            logger << "  # of mesh vertices: " << nmesh;
            this->printLog(logger);

            MVCHierarchyList *hierarchyList = new MVCHierarchyList(vitalPoints);

            logger << "  Hierarchical list generation finish";
            this->printLog(logger);

            for (int channel = 0; channel < 3; ++channel) {
                double *curChannel;
                if (channel == 0) curChannel = dR;
                if (channel == 1) curChannel = dG;
                if (channel == 2) curChannel = dB;

                if (channel == 0) logger << "    Channel R";
                if (channel == 1) logger << "    Channel G";
                if (channel == 2) logger << "    Channel B";
                this->printLog(logger);

                for (CDT::Finite_vertices_iterator vertex_ite = cdt.finite_vertices_begin(); vertex_ite != cdt.finite_vertices_end(); ++vertex_ite) {
                     vertex_ite->info() = 0;
                     meshPoints.push_back(vertex_ite->point());
                }
                for (int i=0; i<nv; ++i) {
                    Face_handle face = cdt.locate(vitalPoints[i]);
                    for (int j=0; j<3; ++j) {
                        Vertex_handle vertex = face->vertex(j);
                        if (distance(vertex->point(), vitalPoints[i]) < 1e-6) {
                            vertex->info() = i + 1;
                        }
                    }
                }

                for (int i=0; i<nv; ++i) {
                    int ns_x = vitalPoints[i].x(), ns_y = vitalPoints[i].y();
                    int x1, y1, x2, y2;
                    if (ns_x & 1)
                        x1 = x2 = (ns_x - 1) >> 1;
                    else
                        x1 = (ns_x >> 1) - 1, x2 = ns_x >> 1;
                    if (ns_y & 1)
                        y1 = y2 = (ns_y - 1) >> 1;
                    else
                        y1 = (ns_y >> 1) - 1, y2 = ns_y >> 1;
                    int x3, y3;
                    if ((z[Q(x1,y1,m)] >= 0) && (z[Q(x1,y1,m)] <= ii)) x3 = x1, y3 = y1;
                    if ((z[Q(x2,y2,m)] >= 0) && (z[Q(x2,y2,m)] <= ii)) x3 = x2, y3 = y2;
                    double now_color = curChannel[Q(x3,y3,m)], latentColor = 255.0;
                    for (int j=ii + 1; j < int(doc->layers.size()); ++j) {
                        DocumentLayer *curLayer = doc->layers[j];
                        int lx = x3 - curLayer->hOffset, ly = y3 - curLayer->wOffset;
                        if ((lx >= 0) && (lx < curLayer->height) && (ly >= 0) && (ly < curLayer->width))
                            if (!curLayer->D[Q(lx,ly,curLayer->width)]) {
                                if (channel == 0) latentColor = curLayer->R[Q(lx,ly,curLayer->width)];
                                if (channel == 1) latentColor = curLayer->G[Q(lx,ly,curLayer->width)];
                                if (channel == 2) latentColor = curLayer->B[Q(lx,ly,curLayer->width)];
                                break;
                            }
                    }
                    vitalDelta[i] = latentColor - now_color;
                }
                logger << "      Boundary difference calculation finish";
                this->printLog(logger);

                hierarchyList->loadValue(vitalDelta, nv);

                for (int i=0; i<nmesh; ++i) {
                    Face_handle face = cdt.locate(meshPoints[i]);
                    for (int j=0; j<3; ++j) {
                        Vertex_handle vertex = face->vertex(j);
                        if (distance(vertex->point(), meshPoints[i]) < 1e-6) {
                            if (vertex->info() > 0) {
                                meshDelta[i] = vitalDelta[vertex->info() - 1];
                            } else {
                                meshDelta[i] = hierarchyList->calc(meshPoints[i].x(), meshPoints[i].y());
                            }
                        }
                    }
                }
                logger << "      Mesh vertices difference calculation finish";
                this->printLog(logger);

                for (int i=0; i<nmesh; ++i) {
                    Face_handle face = cdt.locate(meshPoints[i]);
                    for (int j=0; j<3; ++j) {
                        Vertex_handle vertex = face->vertex(j);
                        if (distance(vertex->point(), meshPoints[i]) < 1e-6)
                            vertex->info() = i + 1;
                    }
                }

                Face_handle face;
                bool initial = true;
                for (int i=0; i<n; ++i)
                    for (int j=0; j<m; ++j)
                        if (floodfillVisit[Q(i,j,m)] == floodfillCnt) {
                            if (initial)
                                face = cdt.locate(Point(2*i+1, 2*j+1)), initial = false;
                            else
                                face = cdt.locate(Point(2*i+1, 2*j+1), face);
                            double delta = triangleInterpolate(2*i+1, 2*j+1, face, meshDelta);
                            curChannel[Q(i,j,m)] += delta;
                        }
                logger << "      Inner pixels difference calculation finish";
                this->printLog(logger);
            }

            for (int i=0; i<n; ++i)
                for (int j=0; j<m; ++j)
                    if (floodfillVisit[Q(i,j,m)] == floodfillCnt) {
                        zPixelCount[z[Q(i,j,m)]]--;
                        z[Q(i,j,m)] = zToFill;
                        zPixelCount[zToFill]++;
                    }

            ++tryTimes;

            delete[] vitalDelta;
            delete[] meshDelta;
            delete hierarchyList;
        }
    }

    for (int i=0; i<nm; ++i) {
        newR[i] = int(dR[i] + 0.5);
        newR[i] = max(min(newR[i], 255), 0);
        newG[i] = int(dG[i] + 0.5);
        newG[i] = max(min(newG[i], 255), 0);
        newB[i] = int(dB[i] + 0.5);
        newB[i] = max(min(newB[i], 255), 0);
    }

    doc->addLayerAtBottom(new DocumentLayer(meshR, meshG, meshB, m, n, "AdaptiveMesh", meshD));
    doc->addLayerAtBottom(new DocumentLayer(newR, newG, newB, m, n, "Result"));

    delete[] oR;
    delete[] oG;
    delete[] oB;
    delete[] z;
    delete[] dR;
    delete[] dG;
    delete[] dB;
    delete[] newR;
    delete[] newG;
    delete[] newB;
    delete[] zPixelCount;
    delete[] traverseVisit;
    delete[] meshR;
    delete[] meshG;
    delete[] meshB;
    delete[] meshD;
    delete[] que_x;
    delete[] que_y;

    ans = doc;

    logger << "Finish";
    this->printLog(logger);
}

void MVCCompositingThread::printLog(stringstream &logger) {
    logger.flush();
    emit this->updateStatus(QString("[") + QString::number((int)(time(0) - this->startTime)) + QString(" s] ") + QString::fromStdString(logger.str()));
    logger.str("");
}
