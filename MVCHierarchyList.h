#ifndef MVCHIERARCHYLIST_H
#define MVCHIERARCHYLIST_H

#include <cmath>

#include <vector>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_conformer_2.h>
#include <CGAL/Delaunay_mesher_2.h>
#include <CGAL/Delaunay_mesh_face_base_2.h>
#include <CGAL/Delaunay_mesh_size_criteria_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Point_set_2.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Triangulation_vertex_base_with_info_2<unsigned int, K> Vb;
typedef CGAL::Delaunay_mesh_face_base_2<K> Fb;
typedef CGAL::Triangulation_data_structure_2<Vb, Fb> Tds;
typedef CGAL::Constrained_Delaunay_triangulation_2<K, Tds> CDT;
typedef CDT::Point Point;

using namespace std;

struct Element {
    double x, y;
    double v;

    int froms[3];
    double fromWs[3];
    int fromCnt;

    Element(double x, double y);
};

double dist2(const Element &a, const Element &b);

double dist(const Element &a, const Element &b);

double getAng(const Element &a, const Element &b, const Element &c);

class MVCHierarchyList
{
public:
    MVCHierarchyList(vector<Point> &pointVec);
    void loadValue(double *varr, int n);
    double calc(double x, double y);
    ~MVCHierarchyList();

    int tot_n;

private:
    pair<double, double> calc(double x, double y, int p, int level);

    vector<vector<Element>*> layers;
};

#endif // MVCHIERARCHYLIST_H
