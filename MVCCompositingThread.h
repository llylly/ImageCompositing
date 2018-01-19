#ifndef MVCCOMPOSITINGTHREAD_H
#define MVCCOMPOSITINGTHREAD_H

#include <cassert>
#include <cmath>
#include <sstream>
#include <QObject>
#include <QThread>
#include "Document.h"
#include "MVCHierarchyList.h"

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_conformer_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Delaunay_mesher_2.h>
#include <CGAL/Delaunay_mesh_face_base_2.h>
#include <CGAL/Delaunay_mesh_size_criteria_2.h>
#include <CGAL/Polygon_2.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Triangulation_vertex_base_with_info_2<unsigned int, K> Vb;
typedef CGAL::Delaunay_mesh_face_base_2<K> Fb;
typedef CGAL::Triangulation_data_structure_2<Vb, Fb> Tds;
typedef CGAL::Constrained_Delaunay_triangulation_2<K, Tds> CDT;
typedef CGAL::Delaunay_mesh_criteria_2<CDT> Criteria;
typedef CGAL::Delaunay_mesher_2<CDT, Criteria> Mesher;
typedef CGAL::Polygon_2<K> Polygon_2;
typedef CDT::Vertex_handle Vertex_handle;
typedef CDT::Face_handle Face_handle;
typedef CDT::Point Point;

using namespace std;

#define Q(i,j,h) ((i)*(h)+(j))

class MVCCompositingThread : public QThread {
    Q_OBJECT

public:
    Document *doc;
    Document *ans;

protected:
    void run();

private:
    time_t startTime;
    void printLog(stringstream &logger);

signals:
    void updateStatus(QString s);
};

#endif // MVCCOMPOSITINGTHREAD_H
