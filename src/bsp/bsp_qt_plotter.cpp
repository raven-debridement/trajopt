#include "bsp_qt_plotter.hpp"

namespace BSP {

  BSPQtPlotter::BSPQtPlotter(double x_min, double x_max, double y_min, double y_max, BSPProblemHelperBasePtr helper, QWidget* parent) :
   QWidget(parent), x_min(x_min), x_max(x_max), y_min(y_min), y_max(y_max), helper(helper) {
    QPalette pal = palette();
    pal.setColor(backgroundRole(), Qt::white);
    setPalette(pal);
  }

  double BSPQtPlotter::scale_x_length(double x) const {
    return x / (x_max - x_min) * width();
  }

  double BSPQtPlotter::scale_y_length(double y) const {
    return y / (y_max - y_min) * height();
  }

  double BSPQtPlotter::scale_x(double x) const {
    return scale_x_length(x - x_min);
  }

  double BSPQtPlotter::scale_y(double y) const {
    return height() - scale_y_length(y - y_min);
  }

  double BSPQtPlotter::unscale_x(double x) const {
    return (x / width()) * (x_max - x_min) + x_min;
  }

  double BSPQtPlotter::unscale_y(double y) const {
    return ((height() - y) / height()) * (y_max - y_min) + y_min;
  }

  void BSPQtPlotter::draw_ellipse(const Vector2d& mean, const Matrix2d& cov, QPainter& painter, double scale_factor) {
    assert(mean.size() == 2);
    assert(cov.rows() == 2 && cov.cols() == 2);
    double cx = mean(0), cy = mean(1);
    double sx = cov(0, 0), sy = cov(1, 1), sxy = cov(0, 1);
    double sxy_diff = sx - sy;
    double ratio = height() * 1.0 / width();
    double theta = 0.5 * atan2(2*sxy, sxy_diff*ratio);
    EigenSolver<MatrixXd> es(cov);
    VectorXd eigenvals = es.eigenvalues().real();
    double max_val = eigenvals(0), min_val = eigenvals(1);
    if (min_val > max_val) {
      double tmp = max_val; max_val = min_val; min_val = tmp;
    }
    double mx, my; // semi major/minor axis lengths
    if (sx >= sy) {
      mx = sqrt(max_val)*scale_factor;
      my = sqrt(min_val)*scale_factor;
    } else {
      mx = sqrt(min_val)*scale_factor;
      my = sqrt(max_val)*scale_factor;
    }
    return draw_ellipse(cx, cy, mx, my, theta, painter);
  }

  void BSPQtPlotter::draw_ellipse(double cx, double cy, double mx, double my, double theta, QPainter& painter) {
    cx = scale_x(cx);
    cy = scale_y(cy);
    mx = scale_x_length(mx);
    my = scale_y_length(my);
    painter.save();
    painter.translate(cx-mx, cy-my);
    painter.rotate(theta);
    painter.drawEllipse(0, 0, 2*mx, 2*my);
    painter.restore();
  }

  void BSPQtPlotter::draw_line(double x1, double y1, double x2, double y2, QPainter& painter) {
    x1 = scale_x(x1);
    y1 = scale_y(y1);
    x2 = scale_x(x2);
    y2 = scale_y(y2);
    painter.drawLine(x1, y1, x2, y2);
  }

  void BSPQtPlotter::draw_point(double x, double y, QPainter& painter) {
    x = scale_x(x);
    y = scale_y(y);
    painter.drawPoint(x, y);
  }

  void BSPQtPlotter::keyPressEvent(QKeyEvent *event) {
    if (event->key() == Qt::Key_P) {
      emit(proceed_signal());
    }
  }

  void BSPQtPlotter::closeEvent(QCloseEvent*) {
    emit finished_signal();
  }

}
