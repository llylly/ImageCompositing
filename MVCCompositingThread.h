#ifndef MVCCOMPOSITINGTHREAD_H
#define MVCCOMPOSITINGTHREAD_H

#include <cassert>
#include <sstream>
#include <QObject>
#include <QThread>
#include "Document.h"

using namespace std;

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
