#include "MVCCompositingThread.h"

void MVCCompositingThread::run() {
    assert(doc != NULL);
    stringstream logger;
    logger << "Begin";
    this->startTime = time(0);
    this->printLog(logger);

    Document *docAns = new Document(doc->theWindow, doc->width, doc->height);

    /* TODO */

    ans = docAns;
    logger << "Finish";
    this->printLog(logger);
}

void MVCCompositingThread::printLog(stringstream &logger) {
    logger.flush();
    emit this->updateStatus(QString("[") + QString::number((int)(time(0) - this->startTime)) + QString(" s] ") + QString::fromStdString(logger.str()));
    logger.str("");
}
