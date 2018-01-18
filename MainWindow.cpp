#include "MainWindow.h"

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent) {
    if (this->objectName().isEmpty())
        this->setObjectName(Context::projectName.c_str());
    this->resize(Context::defaultWindowWidth, Context::defaultWindowHeight);
    this->move((QApplication::desktop()->width() - this->width()) / 2, (QApplication::desktop()->height() - this->height()) / 2);
    this->setWindowTitle(Context::projectName.c_str());

    this->centralWidget = new QWidget(this);
    this->centralWidget->setEnabled(true);
    this->mainLayout = new QHBoxLayout(this->centralWidget);
    this->mainLayout->setSpacing(6);
    this->mainLayout->setContentsMargins(11, 11, 11, 11);
        this->imageView = new ImageView(this);
        this->mainLayout->addWidget(this->imageView, 3);

        this->rightLayout = new QGridLayout();
            this->newButton = new QPushButton(this);
            this->newButton->setText("New Document");
            this->rightLayout->addWidget(this->newButton, 0,0,1,2);

            QFrame *line;
            line = new QFrame();
            line->setFrameShape(QFrame::HLine);
            line->setFrameShadow(QFrame::Sunken);
            this->rightLayout->addWidget(line, 1,0,1,2);

            this->rightLayout->addWidget(new QLabel("Size"), 2,0,1,2);

            this->adjustSizeButton = new QPushButton(this);
            this->adjustSizeButton->setText("Adjust");
            this->rightLayout->addWidget(this->adjustSizeButton, 3,0,1,1);

            this->autoAdjustSizeButton = new QPushButton(this);
            this->autoAdjustSizeButton->setText("Auto Adjust");
            this->rightLayout->addWidget(this->autoAdjustSizeButton, 3,1,1,1);

            line = new QFrame();
            line->setFrameShape(QFrame::HLine);
            line->setFrameShadow(QFrame::Sunken);
            this->rightLayout->addWidget(line, 4,0,1,2);

            this->addLayerButton = new QPushButton(this);
            this->addLayerButton->setText("Add New Layer");
            this->rightLayout->addWidget(this->addLayerButton, 5,0,1,2);

            line = new QFrame();
            line->setFrameShape(QFrame::HLine);
            line->setFrameShadow(QFrame::Sunken);
            this->rightLayout->addWidget(line, 6,0,1,2);

            this->rightLayout->addWidget(new QLabel("Composite"), 7,0,1,2);

            this->compositeButton = new QPushButton(this);
            this->compositeButton->setText("Quadtree Gradient Domain Composite");
            this->rightLayout->addWidget(this->compositeButton, 8,0,1,2);

            this->mvcCompositeButton = new QPushButton(this);
            this->mvcCompositeButton->setText("Mean Value Coordinates Composite");
            this->rightLayout->addWidget(this->mvcCompositeButton, 9,0,1,2);

            line = new QFrame();
            line->setFrameShape(QFrame::HLine);
            line->setFrameShadow(QFrame::Sunken);
            this->rightLayout->addWidget(line, 10,0,1,2);

            this->saveButton = new QPushButton(this);
            this->saveButton->setText("Save");
            this->rightLayout->addWidget(this->saveButton, 11,0,1,2);

            this->layerView = new LayerListWidget(this);
            this->rightLayout->addWidget(this->layerView, 12,0,1,2);

            this->aboutButton = new QPushButton(this);
            this->aboutButton->setText("About");
            this->rightLayout->addWidget(this->aboutButton, 13,1,1,1);
        this->mainLayout->addLayout(this->rightLayout, 1);

    this->setCentralWidget(this->centralWidget);

    this->newDialog = new NewDialog(this);
    this->adjustSizeDialog = new AdjustSizeDialog(this);
    this->statusDialog = new StatusDialog(this);
    this->aboutDialog = new AboutDialog(this);

    this->curDocment = NULL;

    this->quadTreeCompositingThread = new QuadTreeCompositing();
    this->mvcCompositingThread = new MVCCompositingThread();

    this->newButton->setFocus();

    connect(this->newButton, SIGNAL(clicked(bool)), this, SLOT(newClick()));
    connect(this->adjustSizeButton, SIGNAL(clicked(bool)), this, SLOT(adjustSizeClick()));
    connect(this->autoAdjustSizeButton, SIGNAL(clicked(bool)), this, SLOT(autoAdjustClick()));
    connect(this->addLayerButton, SIGNAL(clicked(bool)), this, SLOT(addLayerClick()));
    connect(this->imageView->cancelButton, SIGNAL(clicked(bool)), this, SLOT(viewCancelClick()));
    connect(this->imageView->saveButton, SIGNAL(clicked(bool)), this, SLOT(viewSaveClick()));
    connect(this->imageView->maskCropButton, SIGNAL(clicked(bool)), this, SLOT(viewMaskCropClick()));
    connect(this->imageView->label, SIGNAL(dragged(QPoint)), this, SLOT(dragHandler(QPoint)));
    connect(this->imageView->label, SIGNAL(staticClick(QPoint)), this, SLOT(clickhandler(QPoint)));
    connect(this->compositeButton, SIGNAL(clicked(bool)), this, SLOT(quadTreeCompositingToggle()));
    connect(this->mvcCompositeButton, SIGNAL(clicked(bool)), this, SLOT(mvcCompositingToggle()));
    connect(this->saveButton, SIGNAL(clicked(bool)), this, SLOT(saveClick()));
    connect(this->aboutButton, SIGNAL(clicked(bool)), this, SLOT(aboutClick()));

    connect(this->quadTreeCompositingThread, SIGNAL(updateStatus(QString)), this->statusDialog, SLOT(addStatus(QString)), Qt::QueuedConnection);
    connect(this->quadTreeCompositingThread, SIGNAL(finished()), this->statusDialog, SLOT(close()), Qt::QueuedConnection);
    connect(this->quadTreeCompositingThread, SIGNAL(finished()), this, SLOT(quadTreeCompositingFinish()), Qt::QueuedConnection);

    connect(this->mvcCompositingThread, SIGNAL(updateStatus(QString)), this->statusDialog, SLOT(addStatus(QString)), Qt::QueuedConnection);
    connect(this->mvcCompositingThread, SIGNAL(finished()), this->statusDialog, SLOT(close()), Qt::QueuedConnection);
    connect(this->mvcCompositingThread, SIGNAL(finished()), this, SLOT(mvcCompositingFinish()), Qt::QueuedConnection);
}

void MainWindow::newClick() {
    this->newDialog->available = false;
    this->newDialog->exec();
    if (this->newDialog->available) {
        if (this->curDocment != NULL)
            delete this->curDocment;
        this->curDocment = new Document(this, this->newDialog->width, this->newDialog->height);
        this->layerView->setNewDocument(this->curDocment);
        this->viewUpdate();
    }
}

void MainWindow::adjustSizeClick() {
    if (this->curDocment == NULL) {
        QMessageBox::information(this, "New Document Required", "Please create new document at first.", QMessageBox::Ok);
        return;
    }
    this->adjustSizeDialog->initRefresh(this->curDocment->width, this->curDocment->height);
    this->adjustSizeDialog->exec();
    if (this->adjustSizeDialog->available) {
        this->curDocment->width = this->adjustSizeDialog->width;
        this->curDocment->height = this->adjustSizeDialog->height;
        this->viewUpdate();
    }
}

void MainWindow::autoAdjustClick() {
    if (this->curDocment == NULL) {
        QMessageBox::information(this, "New Document Required", "Please create new document at first.", QMessageBox::Ok);
        return;
    }
    this->curDocment->autoAdjustSize();
    this->viewUpdate();
}

void MainWindow::addLayerClick() {
    if (this->curDocment == NULL) {
        QMessageBox::information(this, "New Document Required", "Please create new document at first.", QMessageBox::Ok);
        return;
    }

    QString fileName = QFileDialog::getOpenFileName(this, "Add Layer From Image File", "",
                                                    "Images (*.bmp *.png *.jpg *.jpeg *.tif *.GIF)");
    if (!fileName.isEmpty()) {
        QImage *qimg = new QImage();
        if (qimg->load(fileName)) {
            string layerName = string("Layer ") + char('A' + this->curDocment->layers.size());
            DocumentLayer *newLayer = new DocumentLayer(qimg, layerName);
            this->curDocment->addLayer(newLayer);
            this->viewUpdate();
        }
        delete qimg;
    }
}

void MainWindow::viewSaveClick() {
    if (this->curDocment == NULL) return;
    this->curDocment->crop();
    this->viewUpdate();
}

void MainWindow::viewCancelClick() {
    if (this->curDocment == NULL) return;
    this->curDocment->unSelect();
    this->curDocment->deletePoints();
    this->viewUpdate();
}

void MainWindow::viewMaskCropClick() {
    if (this->curDocment == NULL) return;
    if (this->curDocment->selectedLayer == -1) return;
    QString fileName = QFileDialog::getOpenFileName(this, "Add Layer From Image File", "",
                                                    "Images (*.bmp *.png *.jpg *.jpeg *.tif *.GIF)");
    if (!fileName.isEmpty()) {
        QImage *qimg = new QImage();
        if (qimg->load(fileName)) {
            this->curDocment->cropByImageMask(qimg);
            this->viewUpdate();
        }
        delete qimg;
    }
}

void MainWindow::saveClick() {
    if (this->curDocment == NULL) {
        QMessageBox::information(this, "New Document Required", "Please create new document at first.", QMessageBox::Ok);
        return;
    }

    QString fileName = QFileDialog::getSaveFileName(
                this,
                QString("Save image to"),
                "",
                QString("BMP Images (*.bmp);;JPG Images (*.jpg);;JPEG Images (*.jpeg);;PNG Images (*.png)"));
    if (!fileName.isEmpty()) {
        this->curDocment->save(fileName.toStdString());
    }
}

void MainWindow::aboutClick() {
    this->aboutDialog->exec();
}

void MainWindow::dragHandler(QPoint delta) {
    if ((this->curDocment == NULL) || (this->curDocment->selectedLayer == -1)) return;
    this->curDocment->layers[this->curDocment->selectedLayer]->addOffset(delta.x(), delta.y());
    this->viewUpdate();
}

void MainWindow::clickhandler(QPoint pos) {
    if ((this->curDocment == NULL) || (this->curDocment->selectedLayer == -1)) return;
    this->curDocment->addPoint(pos);
    this->viewUpdate();
}

void MainWindow::quadTreeCompositingToggle() {
    if (this->curDocment == NULL) {
        QMessageBox::information(this, "New Document Required", "Please create new document at first.", QMessageBox::Ok);
        return;
    }
    this->quadTreeCompositingThread->doc = this->curDocment;
    this->quadTreeCompositingThread->start();
    this->statusDialog->cleanStatus();
    this->statusDialog->exec();
}

void MainWindow::mvcCompositingToggle() {
    if (this->curDocment == NULL) {
        QMessageBox::information(this, "New Document Required", "Please create new document at first.", QMessageBox::Ok);
        return;
    }
    this->mvcCompositingThread->doc = this->curDocment;
    this->mvcCompositingThread->start();
    this->statusDialog->cleanStatus();
    this->statusDialog->exec();
}

void MainWindow::quadTreeCompositingFinish() {
    delete this->curDocment;
    this->curDocment = this->quadTreeCompositingThread->ans;
    this->layerView->setNewDocument(this->curDocment);
    this->viewUpdate();
}

void MainWindow::mvcCompositingFinish() {
    delete this->curDocment;
    this->curDocment = this->mvcCompositingThread->ans;
    this->layerView->setNewDocument(this->curDocment);
    this->viewUpdate();
}

void MainWindow::viewUpdate() {
    if (this->curDocment == NULL) return;
    this->imageView->setImage(this->curDocment->genViewImage());
    this->imageView->updateView();
    this->layerView->updateView();
    if (this->curDocment->selectedLayer != -1)
        this->imageView->showButton();
    else
        this->imageView->hideButton();
}

MainWindow::~MainWindow() {

}
