#pragma once

#include "bsp.hpp"
#include "toy.hpp"
#include <QtCore>
#include <QImage>
#include <qwt_plot.h>
#include <qwt_plot_layout.h>
#include <qwt_plot_canvas.h>
#include <QPainter>
#include <QKeyEvent>
#include <qwt_painter.h>
#include <qevent.h>

namespace ToyBSP {
  class ToyPlotter: public QWidget {
    Q_OBJECT
  public:
    ToyPlotter(QWidget * = NULL);
  public slots:
    void plot(DblVec&, ToyBSPProblemHelperPtr);
  signals:
    void inited();
    void finished();
    void proceed_signal();
  protected:
    virtual void closeEvent(QCloseEvent*);
    virtual void paintEvent(QPaintEvent*);
    virtual void keyPressEvent(QKeyEvent* event);
    virtual void showEvent(QShowEvent*);
    vector<VectorXd> ellipse_params;
    double old_alpha, cur_alpha;
    QImage distmap;
    ToyBSPProblemHelperPtr helper;
  };
}
