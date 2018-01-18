#ifndef QUADTREECOMPOSITING_H
#define QUADTREECOMPOSITING_H

#include <ctime>
#include <QThread>
#include <QObject>
#include <eigen3/Eigen/Sparse>
#include "Document.h"
#include "QuadTreeNode.h"

class QuadTreeCompositing: public QThread {
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

#endif // QUADTREECOMPOSITING_H
