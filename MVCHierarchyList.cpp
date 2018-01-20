#include "MVCHierarchyList.h"

Element::Element(double x, double y) {
    this->x = x, this->y = y, this->v = 0.0;
    this->froms[0] = this->froms[1] = this->froms[2] = -1;
    this->fromWs[0] = this->fromWs[1] = this->fromWs[2] = 0.0;
    this->fromCnt = 0;
}

double dist2(const Element &a, const Element &b) {
    return (a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y);
}

double dist(const Element &a, const Element &b) {
    return sqrt(double(dist2(a, b)));
}

double getAng(const Element &a, const Element &b, const Element &c) {
    double x1 = a.x - b.x, y1 = a.y - b.y;
    double x2 = c.x - b.x, y2 = c.y - b.y;
    double e1 = x1 * x2 + y1 * y2;
    double e2 = sqrt((x1 * x1 + y1 * y1) * (x2 * x2 + y2 * y2));
    if (abs(e2) - abs(e1) < 1e-6) {
        if (e1 * e2 > 0.0)
            return 0.0;
        else
            return M_PI - 1e-6;
    }
    if (abs(e2) < 1e-6) return M_PI - 1e-6;
    if (!isfinite(acos(e1/e2))) {
        cerr << "x1: " << x1 << " x2: " << x2 << " y1: " << y1 << " y2: " << y2 << " e1: " << e1 << " e2: " << e2 << " ERROR IN ACOS!" << endl;
    }
    return acos(e1/e2);
}

MVCHierarchyList::MVCHierarchyList(vector<Point> &pointVec) {
    vector<Element> *bottom = new vector<Element>();
    int n = pointVec.size();
    this->tot_n = n;
    for (int i=0; i<n; ++i) {
        Element now(pointVec[i].x(), pointVec[i].y());
        bottom->push_back(now);
    }
    layers.push_back(bottom);
    int kpow = 1;
    while (layers[layers.size() - 1]->size()  > 16) {
        vector<Element> *source = layers[layers.size() - 1];
        vector<Element> *dest = new vector<Element>();
        kpow *= 64;

        int n = source->size() - 1, p = 1;

        Element start((*source)[0].x, (*source)[0].y);
        start.fromCnt = 1, start.froms[0] = 0, start.fromWs[0] = 1.;
        dest->push_back(start);

        while (p < n) {
            int t = p;
            while ((t + 1 < n) && (t - p < 2) && (dist2((*source)[t], (*source)[t+1]) <= kpow)) {
                ++t;
            }
            int center_t = (p + t) >> 1;
            Element now((*source)[center_t].x, (*source)[center_t].y);
            for (int i=p; i<=t; ++i) {
                now.froms[now.fromCnt] = i;
                now.fromWs[now.fromCnt] = exp(-sqrt(double(dist2((*source)[i], (*source)[center_t]))/double(kpow)));
                now.fromCnt++;
            }
            dest->push_back(now);
            p = t + 1;
        }

        Element end((*source)[n].x, (*source)[n].y);
        end.fromCnt = 1, end.froms[0] = n, end.fromWs[0] = 1.;
        dest->push_back(end);
        layers.push_back(dest);
    }
}

void MVCHierarchyList::loadValue(double *varr, int n) {
    assert(n == int(layers[0]->size()));
    vector<Element> *bottom = layers[0];
    for (int i=0; i<n; ++i)
        (*bottom)[i].v = varr[i];
    for (int i=1; i<int(layers.size()); ++i) {
        vector<Element> *nowl = layers[i];
        vector<Element> *lowl = layers[i-1];
        for (int j=0; j<int(nowl->size()); ++j) {
            Element &ele = (*nowl)[j];
            double sum = 0.;
            for (int k=0; k<ele.fromCnt; ++k)
                sum += ele.fromWs[k];
            ele.v = 0.;
            for (int k=0; k<ele.fromCnt; ++k)
                ele.v += (*lowl)[ele.froms[k]].v * ele.fromWs[k];
            ele.v /= sum;
        }
    }
}

double MVCHierarchyList::calc(double x, double y) {

    vector<Element> *top = layers[layers.size() - 1];
    int top_n = top->size();
    double ans = 0.; double wtot = 0.;
    for (int i=1; i<top_n-1; ++i) {
        pair<double, double> ret = this->calc(x, y, i, 0);
        ans += ret.first; wtot += ret.second;
    }

    return ans / double(wtot);
}

MVCHierarchyList::~MVCHierarchyList() {
    for (int i=0; i<int(this->layers.size()); ++i) {
        delete this->layers[i];
    }
}

pair<double, double> MVCHierarchyList::calc(double x, double y, int p, int level) {
    vector<Element> *curLayer = layers[layers.size() - 1 - level];
    Element &pb = (*curLayer)[p], &pa = (*curLayer)[p-1], &pc = (*curLayer)[p+1];
    Element now(x, y);
    double d = dist(now, pb);
    double ang1 = getAng(pa, now, pb), ang2 = getAng(pb, now, pc);
    double dist_cond = (double)tot_n / ((double)16 * pow(2.5, level));
    double ang_cond = 0.75 * pow(0.8, level);

    double w;
    if (abs(d) < 1e-6)
        w = 1e+20;
    else
        w = (tan(ang1 / 2) + tan(ang2 / 2)) / d;

    if ((level == int(layers.size() - 1)) || ((d > dist_cond) && (ang1 < ang_cond) && (ang2 < ang_cond))) {
        return make_pair(pb.v * w, w);
    } else {
        double sum_v = 0.0;
        double sum_w = 0.0;
        for (int i=0; i<pb.fromCnt; ++i) {
            pair<double, double> ret = this->calc(x, y, pb.froms[i], level+1);
            sum_v += ret.first;
            sum_w += ret.second;
        }
        return make_pair(sum_v, sum_w);
    }
}
