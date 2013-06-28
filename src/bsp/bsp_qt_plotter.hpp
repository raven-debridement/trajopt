#pragma once

#include "common.hpp"
#include "bsp_problem_helper_base.hpp"
#include <QtCore>
#include <QWidget>
#include <qevent.h>
#include <QKeyEvent>
#include <QPainter>

namespace BSP {
  class BSPQtPlotter : public QWidget {
    Q_OBJECT
  public:
    BSPQtPlotter(double x_min, double x_max, double y_min, double y_max, BSPProblemHelperBasePtr helper, QWidget* parent=NULL);
    double scale_x_length(double x) const;
    double scale_y_length(double y) const;
    double scale_x(double x) const;
    double scale_y(double y) const;
    double unscale_x(double x) const;
    double unscale_y(double y) const;
    double x_min;
    double x_max;
    double y_min;
    double y_max;
    BSPProblemHelperBasePtr helper;
    virtual void update_plot_data(void* d1) {}
    virtual void update_plot_data(void* d1, void* d2) {}
    virtual void update_plot_data(void* d1, void* d2, void* d3) {}
  signals:
    void finished_signal();
    void proceed_signal();
  protected:
    virtual void closeEvent(QCloseEvent*);
    virtual void keyPressEvent(QKeyEvent* event);
    void draw_ellipse(const VectorXd& mean, const MatrixXd& cov, QPainter& painter, double scale_factor=1);
    void draw_ellipse(double cx, double cy, double mx, double my, double theta, QPainter& painter);
    void draw_line(double x1, double y1, double x2, double y2, QPainter& painter);
    void draw_point(double x, double y, QPainter& painter);
    vector<VectorXd> states;
    vector<MatrixXd> sigmas;
  };

  typedef boost::shared_ptr<BSPQtPlotter> BSPQtPlotterPtr;
}
